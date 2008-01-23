//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldTest2D.cxx,v 1.3 2008-01-23 23:27:21 adye Exp $
//
// Description:
//      2D test of RooUnfold package using toy MC generated according to PDFs
//      defined in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Fergus Wilson <F.F.Wilson@rl.ac.uk>
//      Tim Adye <T.J.Adye@rl.ac.uk>
//
// Copyright Information:
//      Copyleft (C) 2005-6     Rutherford Appleton Laboratory
//
//==============================================================================

#ifndef NOROOFIT
#define USE_ROOFIT
#endif

#include <cfloat>
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVectorD.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"
#endif

//==============================================================================
// MC generation routine: RooUnfoldTestPdf()
// This routine is included inline so it does not have to be part of the library
// and does not need to be loaded explicitly from the ROOT prompt.
//==============================================================================

#ifdef USE_ROOFIT
#include "RooUnfoldTestPdfRooFit.icc"
#else
#include "RooUnfoldTestPdf.icc"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
Bool_t nosmear= false;
TCanvas* canvas= 0;
TH2D *hTrain2= 0, *hTrainTrue2= 0, *hTrue2= 0, *hMeas2= 0, *hReco2= 0, *hRes2= 0, *hTrue02=   0;
TH1D *hTrainX= 0, *hTrainTrueX= 0, *hTrueX= 0, *hMeasX= 0, *hRecoX= 0, *hPDFx= 0, *hTestPDFx= 0;
TH1D *hTrainY= 0, *hTrainTrueY= 0, *hTrueY= 0, *hMeasY= 0, *hRecoY= 0, *hPDFy= 0, *hTestPDFy= 0;

RooUnfoldResponse* response= 0;
RooUnfold*         unfold=   0;
Int_t regparm= 0, ntoys= 0;

//==============================================================================
// Utility routines
//==============================================================================

//==============================================================================
// Set histogram Y-axis display range
//==============================================================================

void setmax2 (TH1* h,
              const TH1* h1= 0, const TH1* h2= 0, const TH1* h3= 0,
              const TH1* h4= 0, const TH1* h5= 0, const TH1* h6= 0)
{
  Double_t maxval= h1 ? h1->GetMaximum() : -DBL_MAX;
  if (h2 && h2->GetMaximum() > maxval) maxval= h2->GetMaximum();
  if (h3 && h3->GetMaximum() > maxval) maxval= h3->GetMaximum();
  if (h4 && h4->GetMaximum() > maxval) maxval= h4->GetMaximum();
  if (h5 && h5->GetMaximum() > maxval) maxval= h5->GetMaximum();
  if (h6 && h6->GetMaximum() > maxval) maxval= h6->GetMaximum();
  if (maxval > h->GetMaximum()) h->SetMaximum (1.1*maxval);
}

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.7 + (1.0-0.7)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff)  return cutdummy;
  if (nosmear) return xt;   // bin-by-bin correction can't handle bias and smearing
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

Bool_t smear (Double_t& x, Double_t& y)
{
  Double_t xs= smear (x);
  if (xs==cutdummy) return false;
  Double_t ys= smear (y);
  if (ys==cutdummy) return false;
  if (nosmear) {
    x= xs;
    y= ys;
    return true;
  }
  Double_t a= gRandom->Gaus (0.25, 0.05);
  if      (a<0.0) a= 0.0;
  else if (a>1.0) a= 1.0;
  Double_t b= 1.0-a;
  x= a*xs + b*ys;
  y= b*xs + a*xs;
  return true;
}

//==============================================================================
// Train unfolding algorithm
//==============================================================================

Int_t Train2 (Int_t fxtrain, Int_t fytrain, Int_t nb, Int_t ntrain,
              Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi)
{
  const Int_t nbPDF= 500;
  TVectorD xtrue(ntrain), ytrue(ntrain);

  hPDFx= new TH1D ("trainpdfx", "Training PDF X", nbPDF, xlo, xhi);
  hPDFx->SetLineColor(3);
  hPDFx->SetLineWidth(2);
  if (!RooUnfoldTestPdf (fxtrain, ntrain, xlo, xhi, xtrue, hPDFx)) return 0;
  hPDFx->Scale (nbPDF/Double_t(nb));

  hPDFy= new TH1D ("trainpdfy", "Training PDF Y", nbPDF, ylo, yhi);
  hPDFy->SetLineColor(3);
  hPDFy->SetLineWidth(2);
  if (!RooUnfoldTestPdf (fytrain, ntrain, ylo, yhi, ytrue, hPDFy)) return 0;
  hPDFy->Scale (nbPDF/Double_t(nb));

  hTrainTrue2= new TH2D ("traintrue", "Training Truth", nb, xlo, xhi, nb, ylo, yhi);
  hTrainTrue2->SetLineColor(4);
  hTrain2= new TH2D ("train", "Training Measured", nb, xlo, xhi, nb, ylo, yhi);
  hTrain2->SetLineColor(2);

  response->Setup (hTrain2, hTrainTrue2);

  for (Int_t i= 0; i<ntrain; i++) {
    Double_t xt= xtrue[i], yt= ytrue[i];
    hTrainTrue2->Fill (xt, yt);

    Double_t x= xt, y= yt;
    if (smear (x, y)) {
      hTrain2 ->Fill (x, y);
      response->Fill (x, y, xt, yt);
    } else {
      response->Miss (xt, yt);
    }
  }

  hTrainTrueX= hTrainTrue2->ProjectionX("Training Truth X");
  hTrainTrueY= hTrainTrue2->ProjectionY("Training Truth Y");
  hTrainX=     hTrain2    ->ProjectionX("Training Measured X");
  hTrainY=     hTrain2    ->ProjectionY("Training Measured Y");

  setmax2 (hTrainTrueX, hPDFx, hTrainX);
  setmax2 (hTrainTrueY, hPDFy, hTrainY);

  canvas->cd(1);
  hTrainTrueX->Draw();
  hPDFx      ->Draw("LSAME");
  hTrainX    ->Draw("SAME");

  canvas->cd(2);
  hTrainTrueY->Draw();
  hPDFy      ->Draw("LSAME");
  hTrainY    ->Draw("SAME");

  canvas->Update();

  return 1;
}

//==============================================================================
// Test unfolding algorithm
//==============================================================================

Int_t Test2 (Int_t fxtest, Int_t fytest, Int_t nb, Int_t ntest,
             Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi)
{
  const Int_t nbPDF= 500;
  TVectorD xtest(ntest), ytest(ntest);

  hTestPDFx= new TH1D ("pdfx", "PDF X", nbPDF, xlo, xhi);
  hTestPDFx->SetLineColor(3);
  hTestPDFx->SetLineWidth(2);
  if (!RooUnfoldTestPdf (fxtest, ntest, xlo, xhi, xtest, hTestPDFx, 1.0, 2.5)) return 0;
  hTestPDFx->Scale (nbPDF/Double_t(nb));

  hTestPDFy= new TH1D ("pdfy", "PDF Y", nbPDF, ylo, yhi);
  hTestPDFy->SetLineColor(3);
  hTestPDFy->SetLineWidth(2);
  if (!RooUnfoldTestPdf (fytest, ntest, ylo, yhi, ytest, hTestPDFy, -0.5, 2.5)) return 0;
  hTestPDFy->Scale (nbPDF/Double_t(nb));

  hTrue2= new TH2D ("true", "Test Truth", nb, xlo, xhi, nb, ylo, yhi);
  hTrue2->SetLineColor(4);
  hMeas2= new TH2D ("meas", "Test Measured", nb, xlo, xhi, nb, ylo, yhi);
  hMeas2->SetLineColor(2);

  for (Int_t i=0; i<ntest ; i++) {
    Double_t xt= xtest[i], yt= ytest[i];
    hTrue2->Fill (xt, yt);
    Double_t x= xt, y= yt;
    if (smear (x, y))
      hMeas2->Fill (x, y);
  }

  hTrueX= hTrue2->ProjectionX("Test Truth X");
  hTrueY= hTrue2->ProjectionY("Test Truth Y");
  hMeasX= hMeas2->ProjectionX("Test Measured X");
  hMeasY= hMeas2->ProjectionY("Test Measured Y");

  return 1;
}

//==============================================================================
// Unfold
//==============================================================================

void Unfold2 (Int_t method, Int_t nb, Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi)
{
  switch (method) {
    case 1:  unfold= new RooUnfoldBayes    (response, hMeas2, regparm);
             break;
    case 2:  unfold= new RooUnfoldSvd      (response, hMeas2, regparm, ntoys);
             break;
    case 3:  unfold= new RooUnfoldBinByBin (response, hMeas2);
             break;
    default: cerr << "Unknown RooUnfold method " << method << endl;
             return;
  }
  hReco2= (TH2D*) unfold->Hreco();

  hRecoX= hReco2->ProjectionX("Reconstructed X", -1, -1, "E");
  hRecoY= hReco2->ProjectionY("Reconstructed Y", -1, -1, "E");
  hRecoX->SetMarkerStyle(8);
  hRecoY->SetMarkerStyle(8);

  setmax2 (hTrueX, hTestPDFx, hMeasX, hRecoX);
  setmax2 (hTrueY, hTestPDFy, hMeasY, hRecoY);

  canvas->cd(3);
  hTrueX   ->Draw();
  hTestPDFx->Draw("LSAME");
  hMeasX   ->Draw("SAME");
  hRecoX   ->Draw("SAME");

  canvas->cd(4);
  hTrueY   ->Draw();
  hTestPDFy->Draw("LSAME");
  hMeasY   ->Draw("SAME");
  hRecoY   ->Draw("SAME");

  canvas->cd(5);
  hTrue2->Draw();

  canvas->cd(6);
  hMeas2->Draw();

  canvas->cd(7);
  hReco2->Draw();

  // I think hReco2 already includes the statistical error on hTrue2, so
  // don't include that twice when calculating residuals.
  hRes2= new TH2D ("res", "Pull (reco-true)/error", nb, xlo, xhi, nb, ylo, yhi);
  for (Int_t i= 0; i <= nb+1; i++) {
    for (Int_t j= 0; j <= nb+1; j++) {
      Double_t rec= hReco2->GetBinContent (i, j);
      Double_t tru= hTrue2->GetBinContent (i, j);
      Double_t err=  hReco2->GetBinError   (i, j);
      if      (err > 0)
        hRes2->SetBinContent (i, j, (rec-tru)/err);
      else if (rec > 0)
        hRes2->SetBinContent (i, j, (rec-tru)/sqrt(rec));
    }
  }

  canvas->cd(8);
  hRes2->Draw("COLZ");
  canvas->Update();
}

//==============================================================================
// Controlling routine
//==============================================================================

void RooUnfoldTest2D (
                      Int_t    method=      1,
                      Int_t    stage=       0,
                      Int_t    ftrainx=     1,
                      Int_t    ftrainy=     1,
                      Int_t    ftestx=      3,
                      Int_t    ftesty=      5,
                      Int_t    nb=         40,
                      Int_t    ntest=   10000,
                      Int_t    ntrain= 100000,
                      Double_t xlo=     -12.5,
                      Double_t xhi=      10.0,
                      Double_t ylo=     -12.5,
                      Double_t yhi=      10.0,
                      Int_t    regparm_in= -999,  // Bayes niter=4, SVD kterm=20
                      Int_t    ntoys_in= 1000   // SVD only
                     )
{
#ifdef __CINT__
// If run interactively, remove canvas and all histograms that might have been
// created with a previous invocation.
  delete response; response= 0;
  delete unfold;   unfold= 0;
  delete canvas;   canvas= 0;
  gDirectory->Clear();
  hTrain2= hTrainTrue2= hTrue2= hMeas2= hReco2= hRes2= hTrue02=   0;
  hTrainX= hTrainTrueX= hTrueX= hMeasX= hRecoX= hPDFx= hTestPDFx= 0;
  hTrainY= hTrainTrueY= hTrueY= hMeasY= hRecoY= hPDFy= hTestPDFy= 0;
#endif
  regparm= regparm_in!=-999 ? regparm_in : (method==1 ?  4 :
                                            method==2 ? 20 : 0);
  if (method == 3) nosmear= true;  // bin-by-bin can't handle smearing or bias
  ntoys= ntoys_in;
  cout << "RooUnfoldTest2D"
       << " (method="  << method  // RooUnfold method: 1=Bayes, 2=SVD, 3=bin-by-bin
       << ", stage="   << stage   // 1=train (writes RooUnfoldTest2D.root), 2=test (reads), 0=both (default)
       << ", ftrainx=" << ftrainx // selected X training function
       << ", ftestx="  << ftestx  // selected X test function
       << ", ftrainy=" << ftrainy // selected Y training function
       << ", ftesty="  << ftesty  // selected Y test function
       << ", nb="      << nb      // #bins
       << ", ntest="   << ntest   // # events to use for unsmearing
       << ", ntrain="  << ntrain  // #events to use for training
       << ", xlo="     << xlo     // range minimum
       << ", xhi="     << xhi     // range maximum
       << ", ylo="     << ylo     // range minimum
       << ", yhi="     << yhi     // range maximum
       << ", regparm=" << regparm // regularisation parameter (eg. number of iterations)
       << ", ntoys="   << ntoys   // number of toys used to obtain SVD covariances
       << ")" << endl;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  TPostScript ps("RooUnfoldTest2D.ps", 111);
  canvas= new TCanvas("RooUnfoldTest2D","RooUnfoldTest2D",1);
  canvas->SetGrid();
  canvas->Divide(2,4);

  if (stage != 2) {
    response= new RooUnfoldResponse ("response", "test 2-D Unfolding");
    cout   << "==================================== TRAIN ========================" << endl;
    if (Train2 (ftrainx, ftrainy, nb, ntrain, xlo, xhi, ylo, yhi)) {
      TFile f ("RooUnfoldTest2D.root", "recreate");
      f.WriteTObject (response, "response");
      f.Close();
    } else
      return;   // training failed - skip testing
  }

  if (stage != 1) {
    if (!response) {
      TFile f ("RooUnfoldTest2D.root");
      f.GetObject ("response", response);
      if (!response) {
        cerr << "could not read 'response' object from file RooUnfoldTest2D.root" << endl;
        return;
      }
      f.Close();
    }
    cout   << "==================================== TEST =========================" << endl;
    if (Test2 (ftestx, ftesty, nb, ntest,  xlo, xhi, ylo, yhi)) {
      cout << "==================================== UNFOLD =======================" << endl;
      Unfold2 (method,         nb,         xlo, xhi, ylo, yhi);
    }
  }

  ps.Close();
}

//==============================================================================
// Main program when run stand-alone
//==============================================================================

#ifndef __CINT__
int main (int argc, char *argv[])
{
  switch (argc) {
    case  1:  RooUnfoldTest2D(); break;
    case  2:  RooUnfoldTest2D(atoi(argv[1])); break;
    case  3:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2])); break;
    case  4:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3])); break;
    case  5:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])); break;
    case  6:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5])); break;
    case  7:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6])); break;
    case  8:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7])); break;
    case  9:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8])); break;
    case 10:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9])); break;
    case 11:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atof(argv[10])); break;
    case 12:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atof(argv[10]), atof(argv[11])); break;
    case 13:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12])); break;
    case 14:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13])); break;
    case 15:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]), atoi(argv[14])); break;
    case 16:  RooUnfoldTest2D(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atof(argv[10]), atof(argv[11]), atof(argv[12]), atof(argv[13]), atoi(argv[14]), atoi(argv[15])); break;
    default: cerr << argv[0] << ": too many arguments (" << argc-1 << ")" << endl;
             return 1;
  }
  return 0;
}
#endif
