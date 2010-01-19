// File and Version Information:
//      $Id: RooUnfoldTestHarness3D.h,v 1.1 2010-01-19 00:05:55 adye Exp $
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

#ifndef ROOUNFOLDTESTHARNESS3D_HH
#define ROOUNFOLDTESTHARNESS3D_HH

#include "RooUnfoldTestHarness2D.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1.h"
#include "TH3.h"
#endif

class RooUnfoldTestHarness3D : public RooUnfoldTestHarness2D {
public:
  // Parameters
  Int_t    ftrainz, ftestz, ntz, nmz;
  Double_t zlo, zhi;

  TH1D *hTrainZ, *hTrainTrueZ, *hTrueZ, *hMeasZ, *hRecoZ, *hPDFz, *hTestPDFz;

  // Constructors
  RooUnfoldTestHarness3D (const char* name= "RooUnfoldTest3D");
  RooUnfoldTestHarness3D (const char* name, const char* args);
  RooUnfoldTestHarness3D (const char* name, int argc, const char* const* argv);
  virtual ~RooUnfoldTestHarness3D() {}

  virtual void  Reset();
  virtual void  Defaults();
  virtual void  Init();
  virtual Int_t Train();
  virtual Int_t Test();
  virtual void  Results();
  virtual Int_t Check();
  virtual void  Parms (ArgVars& args);

  Double_t smear (Double_t xt, Int_t nt, Double_t xlo, Double_t xhi) { return RooUnfoldTestHarness::smear(xt,nt,xlo,xhi); }
  Bool_t smear (Double_t& x, Double_t& y, Double_t& z,
                Int_t nx, Double_t xlo, Double_t xhi,
                Int_t ny, Double_t ylo, Double_t yhi,
                Int_t nz, Double_t zlo, Double_t zhi);

  Int_t Fill (TH1* h, Double_t x, Double_t y, Double_t z) {return dynamic_cast<TH3*>(h)->Fill (x, y, z);}
  static TH1D* ProjectionX (const TH1* h, const char* name=0, Option_t* opt="")
    {TH1D* g= dynamic_cast<TH1D*>(dynamic_cast<const TH3*>(h)->Project3D(TString("x")+opt)); if (g&&name) g->SetName(name); return g;}
  static TH1D* ProjectionY (const TH1* h, const char* name=0, Option_t* opt="")
    {TH1D* g= dynamic_cast<TH1D*>(dynamic_cast<const TH3*>(h)->Project3D(TString("y")+opt)); if (g&&name) g->SetName(name); return g;}
  static TH1D* ProjectionZ (const TH1* h, const char* name=0, Option_t* opt="")
    {TH1D* g= dynamic_cast<TH1D*>(dynamic_cast<const TH3*>(h)->Project3D(TString("z")+opt)); if (g&&name) g->SetName(name); return g;}
};

#endif
