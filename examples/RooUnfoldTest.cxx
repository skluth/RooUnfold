// File and Version Information:
//      $Id: RooUnfoldTest.cxx,v 1.8 2009-06-12 00:44:40 adye Exp $
//
// Description:
//      Tests RooUnfold package using toy MC generated according to PDFs defined
//      in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
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

#define ONEPAGE 4

#ifndef NOROOFIT
#define USE_ROOFIT
#endif

#include <cfloat>
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <cmath>
#include <iostream>

#include "TROOT.h"
#include "TObject.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TLine.h"

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

#if !defined(__CINT__) || defined(__MAKECINT__)
using std::fabs;
using std::cout;
using std::cerr;
using std::endl;
#endif

//==============================================================================
// Global definitions
//==============================================================================

const char * const methodnames[]= {"bayes", "svd", "binbybin"};
const Double_t cutdummy= -99999.0;
Bool_t nosmear= false;
#ifdef ONEPAGE
Int_t onepage= ONEPAGE;
#else
Int_t onepage= 0;
#endif
Int_t ipad= 0;
TCanvas* canvas= 0;
TH1D *hPDF= 0, *hTrain= 0, *hTrainTrue= 0, *hTestPDF= 0, *hTrue= 0;
TH1D  *hMeas= 0, *hReco= 0, *hTrue0= 0, *hRes= 0, *hPulls= 0;
TH2D* hResmat= 0;
Int_t regparm= 0, ntoys= 0;

RooUnfoldResponse* response= 0;
RooUnfold*         unfold=   0;

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
  // Efficiency: 30% at x=xlo, 100% at x=xhi.
  // Shift = -10% of the range.
  // Smear = half a bin width.

  const Double_t ylo= 0.3, yhi= 1.0, relshift= -0.1, binsmear= 0.5;
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

//==============================================================================
// Train unfolding algorithm
//==============================================================================

Int_t Train (Int_t ftrain, Int_t nm, Int_t nt, Int_t ntrain, Double_t xlo, Double_t xhi)
{
  const Int_t nbPDF= 500;
  TVectorD xtrue(ntrain);
  hPDF= new TH1D ("trainpdf", "Training PDF", nbPDF, xlo, xhi);
  hPDF->SetLineColor(kGreen);
  hPDF->SetLineWidth(2);
  if (!RooUnfoldTestPdf (ftrain, ntrain, xlo, xhi, xtrue, hPDF)) return 0;
  hPDF->Scale (nbPDF/Double_t(nt));


  hTrainTrue= new TH1D ("traintrue", "Training Truth", nt, xlo, xhi);
  hTrainTrue->SetLineColor(kBlue);
  hTrain= new TH1D ("train", "Training Measured", nm, xlo, xhi);
  hTrain->SetLineColor(kRed);
  hResmat= new TH2D ("resmat", "Response Matrix", nm, xlo, xhi, nt, xlo, xhi);

  response->Setup (nm, xlo, xhi, nt, xlo, xhi);
  //  response->Setup (hTrain, hTrainTrue);

  for (Int_t i= 0; i<ntrain; i++) {
    Double_t xt= xtrue[i];
    hTrainTrue->Fill (xt);
    Double_t x= smear (xt, nt, xlo, xhi);
    if (x!=cutdummy) {
      hTrain ->Fill (x);
      hResmat->Fill  (x, xt);
      response->Fill (x, xt);
    } else
      response->Miss (xt);
  }

  //  response->Setup (hTrain, hTrainTrue, hResmat);

  setmax (hTrainTrue, hPDF, hTrain);

  if (onepage>=4) canvas->cd(++ipad);
  if (!onepage || onepage >= 4) {
    hTrainTrue->Draw();
    hPDF      ->Draw("LSAME");
    hTrain    ->Draw("SAME");
    canvas->Update();
  }

  return 1;
}

//==============================================================================
// Test distribution
//==============================================================================

Int_t Test (Int_t ftest, Int_t nm, Int_t nt, Int_t ntest, Double_t xlo, Double_t xhi)
{
  const Int_t nbPDF= 500;
  TVectorD xtest(ntest);
  hTestPDF= new TH1D ("pdf", "PDF", nbPDF, xlo, xhi);
  hTestPDF->SetLineColor(kGreen); // green
  hTestPDF->SetLineWidth(2);
  if (!RooUnfoldTestPdf (ftest, ntest, xlo, xhi, xtest, hTestPDF, 1.0, 2.0)) return 0;
  hTestPDF->Scale (nbPDF/Double_t(nt));

  hTrue= new TH1D ("true", "Test Truth", nt, xlo, xhi);
  hTrue->SetLineColor(kBlue); // Blue
  hMeas= new TH1D ("meas", "Test Measured", nm, xlo, xhi);
  hMeas->SetLineColor(kRed); // Red
  for (Int_t i=0; i<ntest ; i++) {
    Double_t xt= xtest[i];
    hTrue->Fill(xt);
    Double_t x = smear (xt, nt, xlo, xhi);
    if (x!=cutdummy)
      hMeas->Fill(x);
  }

  return 1;
}

//==============================================================================
// Unfold
//==============================================================================

void Unfold (Int_t method, Int_t nm, Int_t nt, Double_t xlo, Double_t xhi)
{
  cout << "Create RooUnfold object for method " << method << endl;
  switch (method) {
    case 1:  unfold= new RooUnfoldBayes    (response, hMeas, regparm);
             break;
    case 2:  unfold= new RooUnfoldSvd      (response, hMeas, regparm, ntoys);
             break;
    case 3:  unfold= new RooUnfoldBinByBin (response, hMeas);
             break;
    default: cerr << "Unknown RooUnfold method " << method << endl;
             return;
  }
  cout << "Created "; unfold->Print();

  // Unfolded distribution
  hReco= (TH1D*) unfold->Hreco();
  unfold->PrintTable (cout, hTrue);

  hReco->SetLineColor(kBlack);
  hReco->SetMarkerStyle(kFullDotLarge);
  setmax (hTrue, hTestPDF, hMeas, hReco);

  if (onepage) canvas->cd(++ipad);
  hTrue   ->Draw();
  hTestPDF->Draw("LSAME");
  hMeas   ->Draw("SAME");
  hReco   ->Draw("SAME");
  // canvas->Update();

  // I think hReco already includes the statistical error on hTrue, so
  // don't include that twice when calculating residuals.
  hTrue0= (TH1D*) hTrue->Clone();
  hTrue0->SetNameTitle("true0", "Truth with zero errors");
  for (Int_t i = 0 ; i <= nt+1 ; i++)
    hTrue0->SetBinError (i, 0.0);

  hRes= new TH1D ("reco-true", "Residuals", nt, xlo, xhi);
  hRes->SetMarkerStyle(kFullDotLarge);
  hRes->Sumw2();
  hRes->Add (hTrue0, hReco, -1, 1);

  if (onepage>=2) canvas->cd(++ipad);
  if (!onepage || onepage >= 2) {
    hRes->Draw();
    TLine l; // draw a line at y=0;
    l.DrawLine(hRes->GetBinLowEdge(1), 0.0, hRes->GetBinLowEdge(nt+1), 0.0);
    //   canvas->Update();
  }

  // how about doing some pulls
  hPulls = new TH1D ("pulls", "Pulls", nt, xlo, xhi);
  const Double_t MAXPULL = 5.0;

  for (Int_t i = 0 ; i <= nt+1; i++) {
    
    Double_t ydiff    = hRes->GetBinContent(i);
    Double_t ydiffErr = hRes->GetBinError(i);

    Double_t ypull;
    if (ydiffErr==0) {
      ypull= ydiff> 0 ? MAXPULL : -MAXPULL;
    } else {
      ypull = ydiff/ydiffErr;
    }

    //cout << ypull << " " << ydiffErr << endl;

    if (fabs(ypull)<=MAXPULL) {
      hPulls->SetBinContent(i,ypull);
      hPulls->SetBinError(i,1.0);
    }
  }

  if (onepage>=2) canvas->cd(++ipad);
  if (!onepage || onepage >= 2) {
    hPulls->SetMarkerStyle(kFullDotLarge);
    hPulls->Draw("P");
    TLine l; // draw a line at y=0;
    l.DrawLine(hPulls->GetBinLowEdge(1), 0.0, hPulls->GetBinLowEdge(nt+1), 0.0);
    //    canvas->Update();
  }
 canvas->Update();

}

//==============================================================================
// Controlling routine
//==============================================================================

void RooUnfoldTest (
                    Int_t    method=      1,
                    Int_t    stage=       0,
                    Int_t    ftrain=      2,
                    Int_t    ftest=       5,
                    Int_t    nt=         40,
                    Int_t    ntest=   10000,
                    Int_t    ntrain= 100000,
                    Double_t xlo=     -12.5,
                    Double_t xhi=      10.0,
                    Int_t    regparm_in= -999,  // Bayes niter=4, SVD kterm=20
                    Int_t    ntoys_in= 1000,   // SVD only
                    Int_t    nm=         -1
                   )
{
#ifdef __CINT__
// If run interactively, remove canvas and all histograms that might have been
// created with a previous invocation.
  delete response; response= 0;
  delete unfold;   unfold=   0;
  delete canvas;   canvas=   0;
  gDirectory->Clear();
  hPDF= hTestPDF= hTrain= hTrainTrue= hResmat= hTrue= hMeas= hReco= hTrue0= hRes= 0;
#endif
  if (nm==-1) nm= nt;
  regparm= regparm_in!=-999 ? regparm_in : (method==1 ?  4 :
                                            method==2 ? (nt>5 ? nt/2 : 2) : 0);
  if (method == 3) nosmear= true;  // bin-by-bin can't handle smearing or bias

  ntoys= ntoys_in;
  cout << "RooUnfoldTest"
       << " (method=" << method  // RooUnfold method: 1=Bayes, 2=SVD, 3=bin-by-bin
       << ", stage="  << stage   // 1=train (writes RooUnfoldTest.root), 2=test (reads), 0=both (default)
       << ", ftrain=" << ftrain  // selected training function
       << ", ftest="  << ftest   // selected test function
       << ", nt="     << nt      // #truth bins
       << ", ntest="  << ntest   // # events to use for unsmearing
       << ", ntrain=" << ntrain  // #events to use for training
       << ", xlo="    << xlo     // range minimum
       << ", xhi="    << xhi     // range maximum
       << ", regparm="<< regparm // regularisation parameter (eg. number of iterations)
       << ", ntoys="  << ntoys   // number of toys used to obtain SVD covariances
       << ", nm="     << nm      // #measured bins
       << ")" << endl;

  // some stupid tests
  Bool_t statusOk = true;
  if (xlo >= xhi) { cout << "Error: xlo (" << xlo << ") >= xhi(" << xhi << ")" << endl; statusOk = false;}
  if (ntest<=0) {cout << "Error: ntest (" << ntest << ") <= 0" << endl; statusOk = false;}
  if (ntrain<=0) {cout << "Error: ntrain (" << ntrain << ") <= 0" << endl; statusOk = false;}
  if (nm<=0) {cout << "Error: nm (" << nm << ") <= 0" << endl; statusOk = false;}
  if (nt<=0) {cout << "Error: nt (" << nt << ") <= 0" << endl; statusOk = false;}
  if (ftest<1) {cout << "Error: ftest (" << ftest << ") < 1" << endl; statusOk = false;}

  if (!statusOk) {exit(1);}

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  TPostScript ps("RooUnfoldTest.ps", 112);
  canvas= new TCanvas("RooUnfoldTest","RooUnfoldTest",1);
  //canvas= new TCanvas("RooUnfoldTest","RooUnfoldTest",10,10,900,500);
  canvas->Clear();

  if (onepage==2) {
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadBottomMargin(0.08);
    gStyle->SetPadLeftMargin(0.06);
    canvas->Divide(1,2);
    TPad* pad1= (TPad*) canvas->GetPad(1);
    TPad* pad2= (TPad*) canvas->GetPad(2);
    pad1->SetPad (pad1->GetXlowNDC(), .3, pad1->GetXlowNDC()+pad1->GetWNDC(), pad1->GetYlowNDC()+pad1->GetHNDC());
    pad2->SetPad (pad2->GetXlowNDC(), pad2->GetYlowNDC(), pad2->GetXlowNDC()+pad2->GetWNDC(), .3);
    pad1->SetGrid(1);
    pad2->SetGrid(1);
  } else if (onepage>=3)
    canvas->Divide(1,onepage);
  else
    canvas->SetGrid(1);
  ipad= 0;

  if (stage != 2) {
    response= new RooUnfoldResponse ("response", "test 1-D Unfolding");
    if (!response) return;
    cout   << "==================================== TRAIN ====================================" << endl;
    if (Train (ftrain, nm, nt, ntrain, xlo, xhi)) {
      TFile f ("RooUnfoldTest.root", "recreate");
      f.WriteTObject (response, "response");
      f.Close();
    } else
      return;   // training failed - skip testing
  }

  if (stage != 1) {
    if (!response) {
      TFile f ("RooUnfoldTest.root");
      f.GetObject ("response", response);
      f.Close();
      if (!response) {
        cerr << "could not read 'response' object from file RooUnfoldTest.root" << endl;
        return;
      }
    }
    cout   << "==================================== TEST =====================================" << endl;
    if (Test (ftest,  nm, nt, ntest,  xlo, xhi)) {
      cout << "==================================== UNFOLD ===================================" << endl;
      Unfold (method, nm, nt,         xlo, xhi);
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
    case  1:  RooUnfoldTest(); break;
    case  2:  RooUnfoldTest(atoi(argv[1])); break;
    case  3:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2])); break;
    case  4:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3])); break;
    case  5:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4])); break;
    case  6:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5])); break;
    case  7:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6])); break;
    case  8:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7])); break;
    case  9:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atof(argv[8])); break;
    case 10:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atof(argv[8]), atof(argv[9])); break;
    case 11:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atof(argv[8]), atof(argv[9]), atoi(argv[10])); break;
    case 12:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atof(argv[8]), atof(argv[9]), atoi(argv[10]), atoi(argv[11])); break;
    case 13:  RooUnfoldTest(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atof(argv[8]), atof(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12])); break;
    default: cerr << argv[0] << ": too many arguments (" << argc-1 << ")" << endl;
             return 1;
  }
  return 0;
}
#endif
