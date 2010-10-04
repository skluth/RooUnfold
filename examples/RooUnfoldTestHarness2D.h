//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Test Harness class for the RooUnfold package using 2D toy MC.
//      Inherits from RooUnfoldTestHarness.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDTESTHARNESS2D_HH
#define ROOUNFOLDTESTHARNESS2D_HH

#include "RooUnfoldTestHarness.h"
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH2.h"
#endif

class TH1D;

class RooUnfoldTestHarness2D : public RooUnfoldTestHarness {
public:
  // Parameters
  Int_t    ftrainy, ftesty, nty, nmy;
  Double_t ylo, yhi, mtrainy, wtrainy, btrainy, mtesty, wtesty, btesty;
  Double_t effylo, effyhi, rotxy, ybias, ysmear;

  // Data
  TH1D *hTrainX, *hTrainTrueX, *hTrueX, *hMeasX, *hRecoX, *hPullsX;
  TH1D *hTrainY, *hTrainTrueY, *hTrueY, *hMeasY, *hRecoY, *hPullsY, *hPDFy, *hTestPDFy;

  // Constructors
  RooUnfoldTestHarness2D (const char* name= "RooUnfoldTest2D");
  RooUnfoldTestHarness2D (const char* name, const char* args);
  RooUnfoldTestHarness2D (const char* name, int argc, const char* const* argv);
  virtual ~RooUnfoldTestHarness2D() {}

  // Methods and functions
  virtual void  Parms (ArgVars& args);
  virtual Int_t Train();
  virtual void  TrainResults();
  virtual Int_t Test();
  virtual void  Results();
  static  void  Rot     (Double_t& x, Double_t& y, Double_t angle, Double_t x0, Double_t y0);
  virtual bool  Eff2D   (Double_t  x, Double_t  y) const;
  virtual void  Smear2D (Double_t& x, Double_t& y) const;
  virtual void  Reset();
  virtual void  Init();
  virtual Int_t CheckParms();

  static TH1D* ProjectionX (const TH1* h, const char* name="_px", const char* title=0, Option_t* opt="")
    {const TH2* h2=dynamic_cast<const TH2*>(h); TH1D* h1= h2->ProjectionX(name,0,-1,opt); if (title) h1->SetTitle(title); return h1;}
  static TH1D* ProjectionY (const TH1* h, const char* name="_py", const char* title=0, Option_t* opt="")
    {const TH2* h2=dynamic_cast<const TH2*>(h); TH1D* h1= h2->ProjectionY(name,0,-1,opt); if (title) h1->SetTitle(title); return h1;}
};

#ifndef NOINLINE
#include "RooUnfoldTestHarness2D.icc"
#endif

#endif
