//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTestHarness2D.h,v 1.6 2010-01-19 23:30:57 adye Exp $
//
// Description:
//      Harness class to test the RooUnfold package using 2D toy MC generated
//      according to PDFs defined in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDTESTHARNESS2D_HH
#define ROOUNFOLDTESTHARNESS2D_HH

#include "RooUnfoldTestHarness.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1D.h"
#include "TH2.h"
#endif

class RooUnfoldTestHarness2D : public RooUnfoldTestHarness {
public:
  // Parameters
  Int_t    ftrainy, ftesty, nty, nmy;
  Double_t ylo, yhi, bkgty, bkgey;

  TH1D *hTrainX, *hTrainTrueX, *hTrueX, *hMeasX, *hRecoX;
  TH1D *hTrainY, *hTrainTrueY, *hTrueY, *hMeasY, *hRecoY, *hPDFy, *hTestPDFy;

  // Constructors
  RooUnfoldTestHarness2D (const char* name= "RooUnfoldTest2D");
  RooUnfoldTestHarness2D (const char* name, const char* args);
  RooUnfoldTestHarness2D (const char* name, int argc, const char* const* argv);
  virtual ~RooUnfoldTestHarness2D() {}

  virtual void  Reset();
  virtual void  Init();
  virtual Int_t Train();
  virtual Int_t Test();
  virtual void  Results();
  virtual Int_t Check();
  virtual void  Parms (ArgVars& args);

  void rot (Double_t& x, Double_t& y);
  Double_t smear (Double_t xt, Int_t nt, Double_t xlo, Double_t xhi) { return RooUnfoldTestHarness::smear(xt,nt,xlo,xhi); }
  Bool_t smear (Double_t& x, Double_t& y, Int_t nx, Double_t xlo, Double_t xhi, Int_t ny, Double_t ylo, Double_t yhi);

  static TH1D* ProjectionX (const TH1* h, const char* name="_px", const char* title=0, Option_t* opt="")
    {const TH2* h2=dynamic_cast<const TH2*>(h); TH1D* h1= h2->ProjectionX(name,0,-1,opt); if (title) h1->SetTitle(title); return h1;}
  static TH1D* ProjectionY (const TH1* h, const char* name="_py", const char* title=0, Option_t* opt="")
    {const TH2* h2=dynamic_cast<const TH2*>(h); TH1D* h1= h2->ProjectionY(name,0,-1,opt); if (title) h1->SetTitle(title); return h1;}
};

#endif
