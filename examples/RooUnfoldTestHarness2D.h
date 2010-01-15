// File and Version Information:
//      $Id: RooUnfoldTestHarness2D.h,v 1.2 2010-01-15 21:02:07 adye Exp $
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

#ifndef ROOUNFOLDTESTHARNESS2D_HH
#define ROOUNFOLDTESTHARNESS2D_HH

#include "RooUnfoldTestHarness.h"
#include "TH1.h"
#include "TH2.h"

class RooUnfoldTestHarness2D : public RooUnfoldTestHarness {
public:
  // Parameters
  Int_t    ftrainx, ftrainy, ftestx, ftesty;
  Double_t ylo, yhi;

  TH1D *hTrainX, *hTrainTrueX, *hTrueX, *hMeasX, *hRecoX, *hPDFx, *hTestPDFx;
  TH1D *hTrainY, *hTrainTrueY, *hTrueY, *hMeasY, *hRecoY, *hPDFy, *hTestPDFy;

  // Constructors
  RooUnfoldTestHarness2D (const char* name= "RooUnfoldTest2D")                 : RooUnfoldTestHarness(name) {}
  RooUnfoldTestHarness2D (const char* name, const char* args)                  : RooUnfoldTestHarness(name,args) {}
  RooUnfoldTestHarness2D (const char* name, int argc, const char* const* argv) : RooUnfoldTestHarness(name,argc,argv) {}
  virtual ~RooUnfoldTestHarness2D() {}

  virtual void  Reset();
  virtual void  Defaults();
  virtual Int_t Train();
  virtual Int_t Test();
  virtual void  Unfold();
  virtual int SetArgs (int argc, const char* const* argv, bool split= false);

  void rot (Double_t& x, Double_t& y);
  Double_t smear (Double_t xt, Int_t nt, Double_t xlo, Double_t xhi) { return RooUnfoldTestHarness::smear(xt,nt,xlo,xhi); }
  Bool_t smear (Double_t& x, Double_t& y, Int_t nx, Double_t xlo, Double_t xhi, Int_t ny, Double_t ylo, Double_t yhi);

  static TH1D* ProjectionX (const TH1* h, const char* name="_px", Int_t firstybin= 0, Int_t lastybin= -1, Option_t* option= "")
    {return dynamic_cast<const TH2*>(h)->ProjectionX(name,firstybin,lastybin,option);}
  static TH1D* ProjectionY (const TH1* h, const char* name="_py", Int_t firstxbin= 0, Int_t lastxbin= -1, Option_t* option= "")
    {return dynamic_cast<const TH2*>(h)->ProjectionY(name,firstxbin,lastxbin,option);}
};

#endif
