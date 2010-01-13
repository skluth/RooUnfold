//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldTest3D.cxx,v 1.2 2010-01-13 00:18:19 adye Exp $
//
// Description:
//      3D test of RooUnfold package using toy MC generated according to PDFs
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
#include "TH3D.h"
#include "TFile.h"
#include "TVectorD.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"
#endif

#include "RooUnfoldTestArgs.icc"

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

// Does TH3D support ProjectX and ProjectY methods?
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,24,0)
#define PROJECT3D 1
#endif

const Double_t cutdummy= -99999.0;
Bool_t nosmear= false;
TCanvas* canvas= 0;
TH3D *hTrain3= 0, *hTrainTrue3= 0, *hTrue3= 0, *hMeas3= 0, *hReco3= 0, *hTrue03=   0;
TH1D *hTrainX= 0, *hTrainTrueX= 0, *hTrueX= 0, *hMeasX= 0, *hRecoX= 0, *hPDFx= 0, *hTestPDFx= 0;
TH1D *hTrainY= 0, *hTrainTrueY= 0, *hTrueY= 0, *hMeasY= 0, *hRecoY= 0, *hPDFy= 0, *hTestPDFy= 0;
TH1D *hTrainZ= 0, *hTrainTrueZ= 0, *hTrueZ= 0, *hMeasZ= 0, *hRecoZ= 0, *hPDFz= 0, *hTestPDFz= 0;

RooUnfoldResponse* response= 0;
RooUnfold*         unfold=   0;
Int_t regparm= 0, ntoys= 0;

//==============================================================================
// Utility routines
//==============================================================================

//==============================================================================
// Set histogram Y-axis display range
//==============================================================================

void setmax (TH1* h,
             const TH1* h1= 0, const TH1* h2= 0, const TH1* h3= 0,
             const TH1* h4= 0, const TH1* h5= 0, const TH1* h6= 0)
{
  // Get the maximum y value of up to 7 histograms
  // Add 10% to match behaviour of ROOT's automatic scaling
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

Double_t smear (Double_t xt, Int_t nb, Double_t xlo, Double_t xhi)
{
  // Apply a gaussian smearing, systematic translation, and an efficiency
  // function to the truth.
  // Efficiency: 70% at x=xlo, 100% at x=xhi.
  // Shift = -10% of the range.
  // Smear = half a bin width.

  const Double_t ylo= 0.7, yhi= 1.0, relshift= -0.1, binsmear= 0.5;
  Double_t xwidth =  (xhi-xlo);

  Double_t slope = (yhi-ylo) / xwidth;
  Double_t yeff= ylo + slope * (xt-xlo);  // efficiency

  // MC test: if random number > eff then reject
  if (gRandom->Rndm() > yeff)  return cutdummy;
  if (nosmear) return xt;   // bin-by-bin correction can't handle bias and smearing

  Double_t xshift = xwidth*relshift;                 // shift
  Double_t xsigma = xwidth*binsmear / Double_t(nb);  // smear sigma

  Double_t xsmear= gRandom->Gaus(xshift, xsigma);     // bias and smear
  //cout << "SMEAR " << xt << " " << xsmear << " " << xwidth << " " << xsigma << endl;
  return xt+xsmear;
}

void rot (Double_t& x, Double_t& y)
{
  Double_t a= gRandom->Gaus (0.25, 0.05);
  if      (a<0.0) a= 0.0;
  else if (a>1.0) a= 1.0;
  Double_t b= 1.0-a;
  Double_t xr= a*x + b*y;
  Double_t yr= b*x + a*y;
  x= xr;
  y= yr;
}

Bool_t smear (Double_t& x, Double_t& y, Double_t& z,
              Int_t nx, Double_t xlo, Double_t xhi,
              Int_t ny, Double_t ylo, Double_t yhi,
              Int_t nz, Double_t zlo, Double_t zhi)
{
  Double_t xs= smear (x, nx, xlo, xhi);
  if (xs==cutdummy) return false;
  Double_t ys= smear (y, ny, ylo, yhi);
  if (ys==cutdummy) return false;
  Double_t zs= smear (z, nz, zlo, zhi);
  if (ys==cutdummy) return false;
  x= xs;
  y= ys;
  z= zs;
  if (nosmear) return true;
  rot (x, y);
  rot (x, z);
  rot (y, z);
  return true;
}

//==============================================================================
// Train unfolding algorithm
//==============================================================================

Int_t Train3 (Int_t fxtrain, Int_t fytrain, Int_t fztrain, Int_t nb, Int_t ntrain,
              Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi, Double_t zlo, Double_t zhi)
{
  const Int_t nbPDF= 500;
  TVectorD xtrue(ntrain), ytrue(ntrain), ztrue(ntrain);

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

  hPDFz= new TH1D ("trainpdfz", "Training PDF Z", nbPDF, zlo, zhi);
  hPDFz->SetLineColor(3);
  hPDFz->SetLineWidth(2);
  if (!RooUnfoldTestPdf (fztrain, ntrain, zlo, zhi, ztrue, hPDFz)) return 0;
  hPDFz->Scale (nbPDF/Double_t(nb));

  hTrainTrue3= new TH3D ("traintrue", "Training Truth", nb, xlo, xhi, nb, ylo, yhi, nb, zlo, zhi);
  hTrainTrue3->SetLineColor(4);
  hTrain3= new TH3D ("train", "Training Measured", nb, xlo, xhi, nb, ylo, yhi, nb, zlo, zhi);
  hTrain3->SetLineColor(2);

  response->Setup (hTrain3, hTrainTrue3);

  for (Int_t i= 0; i<ntrain; i++) {
    Double_t xt= xtrue[i], yt= ytrue[i], zt= ztrue[i];
    hTrainTrue3->Fill (xt, yt, zt);

    Double_t x= xt, y= yt, z= zt;
    if (smear (x, y, z, nb, xlo, xhi, nb, ylo, yhi, nb, zlo, zhi)) {
      hTrain3 ->Fill (x, y, z);
      response->Fill (x, y, z, xt, yt, zt);
    } else {
      response->Miss (xt, yt, zt);
    }
  }

#if PROJECT3D
  hTrainTrueX= hTrainTrue3->ProjectionX("Training Truth X");
  hTrainX=     hTrain3    ->ProjectionX("Training Measured X");
  setmax (hTrainTrueX, hPDFx, hTrainX);
  canvas->cd(1);
  hTrainTrueX->Draw();
  hPDFx      ->Draw("LSAME");
  hTrainX    ->Draw("SAME");

  hTrainTrueY= hTrainTrue3->ProjectionY("Training Truth Y");
  hTrainY=     hTrain3    ->ProjectionY("Training Measured Y");
  setmax (hTrainTrueY, hPDFy, hTrainY);
  canvas->cd(2);
  hTrainTrueY->Draw();
  hPDFy      ->Draw("LSAME");
  hTrainY    ->Draw("SAME");
#endif

  hTrainTrueZ= hTrainTrue3->ProjectionZ("Training Truth Z");
  hTrainZ=     hTrain3    ->ProjectionZ("Training Measured Z");
  setmax (hTrainTrueZ, hPDFz, hTrainZ);
  canvas->cd(3);
  hTrainTrueZ->Draw();
  hPDFz      ->Draw("LSAME");
  hTrainZ    ->Draw("SAME");

  canvas->Update();

  return 1;
}

//==============================================================================
// Test unfolding algorithm
//==============================================================================

Int_t Test3 (Int_t fxtest, Int_t fytest, Int_t fztest, Int_t nb, Int_t ntest,
             Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi, Double_t zlo, Double_t zhi)
{
  const Int_t nbPDF= 500;
  TVectorD xtest(ntest), ytest(ntest), ztest(ntest);

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

  hTestPDFz= new TH1D ("pdfz", "PDF Z", nbPDF, zlo, zhi);
  hTestPDFz->SetLineColor(3);
  hTestPDFz->SetLineWidth(2);
  if (!RooUnfoldTestPdf (fztest, ntest, zlo, zhi, ztest, hTestPDFz, -0.5, 2.5)) return 0;
  hTestPDFz->Scale (nbPDF/Double_t(nb));

  hTrue3= new TH3D ("true", "Test Truth", nb, xlo, xhi, nb, ylo, yhi, nb, zlo, zhi);
  hTrue3->SetLineColor(4);
  hMeas3= new TH3D ("meas", "Test Measured", nb, xlo, xhi, nb, ylo, yhi, nb, zlo, zhi);
  hMeas3->SetLineColor(2);

  for (Int_t i=0; i<ntest ; i++) {
    Double_t xt= xtest[i], yt= ytest[i], zt= ztest[i];
    hTrue3->Fill (xt, yt, zt);
    Double_t x= xt, y= yt, z= zt;
    if (smear (x, y, z, nb, xlo, xhi, nb, ylo, yhi, nb, zlo, zhi))
      hMeas3->Fill (x, y, z);
  }

#if PROJECT3D
  hTrueX= hTrue3->ProjectionX("Test Truth X");
  hMeasX= hMeas3->ProjectionX("Test Measured X");
  hTrueY= hTrue3->ProjectionY("Test Truth Y");
  hMeasY= hMeas3->ProjectionY("Test Measured Y");
#endif
  hTrueZ= hTrue3->ProjectionZ("Test Truth Z");
  hMeasZ= hMeas3->ProjectionZ("Test Measured Z");

  return 1;
}

//==============================================================================
// Unfold
//==============================================================================

void Unfold3 (Int_t method, Int_t nb, Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi, Double_t zlo, Double_t zhi)
{
  cout << "Create RooUnfold object for method " << method << endl;
  switch (method) {
    case 1:  unfold= new RooUnfoldBayes    (response, hMeas3, regparm);
             break;
    case 2:  unfold= new RooUnfoldSvd      (response, hMeas3, regparm, ntoys);
             break;
    case 3:  unfold= new RooUnfoldBinByBin (response, hMeas3);
             break;
    default: cerr << "Unknown RooUnfold method " << method << endl;
             return;
  }
  cout << "Created "; unfold->Print();
  bool withError= !(dynamic_cast<RooUnfoldBayes*>(unfold) && nb > 20);   // nb^4<160000 (causes*effects, where each is nb*nb)
  if (!withError) {
    cerr << "Don't calculate errors - this would take too long ("
         << nb*nb << " bins - skip if more than 400 for Bayes algorithm)" << endl;
  }
  hReco3= (TH3D*) unfold->Hreco(withError);
  unfold->PrintTable (cout, hTrue3, withError);

#if PROJECT3D
  hRecoX= hReco3->ProjectionX("Reconstructed X", 0, 0, 0, 0, "E");
  hRecoX->SetMarkerStyle(8);
  setmax (hTrueX, hTestPDFx, hMeasX, hRecoX);
  canvas->cd(4);
  hTrueX   ->Draw();
  hTestPDFx->Draw("LSAME");
  hMeasX   ->Draw("SAME");
  hRecoX   ->Draw("SAME");

  hRecoY= hReco3->ProjectionY("Reconstructed Y", 0, 0, 0, 0, "E");
  hRecoY->SetMarkerStyle(8);
  setmax (hTrueY, hTestPDFy, hMeasY, hRecoY);
  canvas->cd(5);
  hTrueY   ->Draw();
  hTestPDFy->Draw("LSAME");
  hMeasY   ->Draw("SAME");
  hRecoY   ->Draw("SAME");
#endif

  hRecoZ= hReco3->ProjectionZ("Reconstructed Z", 0, 0, 0, 0, "E");
  hRecoZ->SetMarkerStyle(8);
  setmax (hTrueZ, hTestPDFz, hMeasZ, hRecoZ);
  canvas->cd(6);
  hTrueZ   ->Draw();
  hTestPDFz->Draw("LSAME");
  hMeasZ   ->Draw("SAME");
  hRecoZ   ->Draw("SAME");

  canvas->cd(7);
  hTrue3->Draw();

  canvas->cd(8);
  hMeas3->Draw();

  canvas->cd(9);
  hReco3->Draw();

  canvas->Update();
}

//==============================================================================
// Controlling routine
//==============================================================================

void RooUnfoldTest3D (
                      Int_t    method=      1,
                      Int_t    stage=       0,
                      Int_t    ftrainx=     1,
                      Int_t    ftrainy=     1,
                      Int_t    ftrainz=     1,
                      Int_t    ftestx=      3,
                      Int_t    ftesty=      5,
                      Int_t    ftestz=      5,
                      Int_t    nb=          6,
                      Int_t    ntest=   10000,
                      Int_t    ntrain= 100000,
                      Double_t xlo=     -12.5,
                      Double_t xhi=      10.0,
                      Double_t ylo=     -12.5,
                      Double_t yhi=      10.0,
                      Double_t zlo=     -12.5,
                      Double_t zhi=      10.0,
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
  hTrain3= hTrainTrue3= hTrue3= hMeas3= hReco3= hTrue03=   0;
  hTrainX= hTrainTrueX= hTrueX= hMeasX= hRecoX= hPDFx= hTestPDFx= 0;
  hTrainY= hTrainTrueY= hTrueY= hMeasY= hRecoY= hPDFy= hTestPDFy= 0;
#endif
  regparm= regparm_in!=-999 ? regparm_in : (method==1 ?  4 :
                                            method==2 ? 20 : 0);
  if (method == 3) nosmear= true;  // bin-by-bin can't handle smearing or bias
  ntoys= ntoys_in;
  cout << "RooUnfoldTest3D"
       << " (method="  << method  // RooUnfold method: 1=Bayes, 2=SVD, 3=bin-by-bin
       << ", stage="   << stage   // 1=train (writes RooUnfoldTest3D.root), 2=test (reads), 0=both (default)
       << ", ftrainx=" << ftrainx // selected X training function
       << ", ftestx="  << ftestx  // selected X test function
       << ", ftrainy=" << ftrainy // selected Y training function
       << ", ftesty="  << ftesty  // selected Y test function
       << ", ftrainz=" << ftrainz // selected Z training function
       << ", ftestz="  << ftestz  // selected Z test function
       << ", nb="      << nb      // #bins
       << ", ntest="   << ntest   // # events to use for unsmearing
       << ", ntrain="  << ntrain  // #events to use for training
       << ", xlo="     << xlo     // range minimum
       << ", xhi="     << xhi     // range maximum
       << ", ylo="     << ylo     // range minimum
       << ", yhi="     << yhi     // range maximum
       << ", zlo="     << zlo     // range minimum
       << ", zhi="     << zhi     // range maximum
       << ", regparm=" << regparm // regularisation parameter (eg. number of iterations)
       << ", ntoys="   << ntoys   // number of toys used to obtain SVD covariances
       << ")" << endl;
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  TPostScript ps("RooUnfoldTest3D.ps", 111);
  canvas= new TCanvas("RooUnfoldTest3D","RooUnfoldTest3D",1);
  canvas->SetGrid();
  canvas->Divide(2,5);

  if (stage != 2) {
    response= new RooUnfoldResponse ("response", "test 3-D Unfolding");
    cout   << "==================================== TRAIN ========================" << endl;
    if (Train3 (ftrainx, ftrainy, ftrainz, nb, ntrain, xlo, xhi, ylo, yhi, zlo, zhi)) {
      TFile f ("RooUnfoldTest3D.root", "recreate");
      f.WriteTObject (response, "response");
      f.Close();
    } else
      return;   // training failed - skip testing
  }

  if (stage != 1) {
    if (!response) {
      TFile f ("RooUnfoldTest3D.root");
      f.GetObject ("response", response);
      if (!response) {
        cerr << "could not read 'response' object from file RooUnfoldTest3D.root" << endl;
        return;
      }
      f.Close();
    }
    cout   << "==================================== TEST =========================" << endl;
    if (Test3 (ftestx, ftesty, ftestz, nb, ntest,  xlo, xhi, ylo, yhi, zlo, zhi)) {
      cout << "==================================== UNFOLD =======================" << endl;
      Unfold3 (method,         nb,         xlo, xhi, ylo, yhi, zlo, zhi);
    }
  }

  ps.Close();
}

//==============================================================================
// Parse arguments and set defaults (uses RooUnfoldTestArgs.icc).
// Defaults here should probably match those in RooUnfoldTest3D(method,...) above.
//==============================================================================

int RooUnfoldTest3D (int argc, const char* const* argv, Bool_t split= false)
{
  Int_t    method=      1;
  Int_t    stage=       0;
  Int_t    ftrainx=     1;
  Int_t    ftrainy=     1;
  Int_t    ftrainz=     1;
  Int_t    ftestx=      3;
  Int_t    ftesty=      5;
  Int_t    ftestz=      5;
  Int_t    nb=          6;
  Int_t    ntest=   10000;
  Int_t    ntrain= 100000;
  Double_t xlo=     -12.5;
  Double_t xhi=      10.0;
  Double_t ylo=     -12.5;
  Double_t yhi=      10.0;
  Double_t zlo=     -12.5;
  Double_t zhi=      10.0;
  Int_t    regparm=  -999;
  Int_t    ntoys=    1000;
  const setargs_t args[]= {
    { "method",  &method,  0 },
    { "stage",   &stage,   0 },
    { "ftrainx", &ftrainx, 0 },
    { "ftrainy", &ftrainy, 0 },
    { "ftrainz", &ftrainz, 0 },
    { "ftestx",  &ftestx,  0 },
    { "ftesty",  &ftesty,  0 },
    { "ftestz",  &ftestz,  0 },
    { "nb",      &nb,      0 },
    { "ntest",   &ntest,   0 },
    { "ntrain",  &ntrain,  0 },
    { "xlo",     0,     &xlo },
    { "xhi",     0,     &xhi },
    { "ylo",     0,     &ylo },
    { "yhi",     0,     &yhi },
    { "zlo",     0,     &zlo },
    { "zhi",     0,     &zhi },
    { "regparm", &regparm, 0 },
    { "ntoys",   &ntoys,   0 },
  };
  if (!setargs (args, (sizeof(args)/sizeof(setargs_t)), argc, argv, split)) return 1;
  RooUnfoldTest3D (method, stage, ftrainx, ftrainy, ftrainz, ftestx, ftesty, ftestz, nb, ntest, ntrain, xlo, xhi, ylo, yhi, zlo, zhi, regparm, ntoys);
  return 0;
}

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

int RooUnfoldTest3D (const char* args)
{
  const char* const argv[]= { "RooUnfoldTest3D", args };
  return RooUnfoldTest3D (2, argv, true);
}

//==============================================================================
// Main program when run stand-alone
//==============================================================================

#ifndef __CINT__
int main (int argc, char** argv) { return RooUnfoldTest3D (argc, argv); }
#endif
