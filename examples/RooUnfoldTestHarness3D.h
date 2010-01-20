//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTestHarness3D.h,v 1.4 2010-01-20 15:41:36 adye Exp $
//
// Description:
//      Harness class to test the RooUnfold package using 3D toy MC generated
//      according to PDFs defined in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
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
  Double_t zlo, zhi, bkgtz, bkgez;

  TH1D *hTrainZ, *hTrainTrueZ, *hTrueZ, *hMeasZ, *hRecoZ, *hPDFz, *hTestPDFz;

  // Constructors
  RooUnfoldTestHarness3D (const char* name= "RooUnfoldTest3D");
  RooUnfoldTestHarness3D (const char* name, const char* args);
  RooUnfoldTestHarness3D (const char* name, int argc, const char* const* argv);
  virtual ~RooUnfoldTestHarness3D() {}

  virtual void  Reset();
  virtual void  Init();
  virtual Int_t Train();
  virtual Int_t Test();
  virtual void  Results();
  virtual Int_t CheckParms();
  virtual void  Parms (ArgVars& args);

  Double_t smear (Double_t xt, Int_t nt, Double_t xlo, Double_t xhi) { return RooUnfoldTestHarness::smear(xt,nt,xlo,xhi); }
  Bool_t smear (Double_t& x, Double_t& y, Double_t& z,
                Int_t nx, Double_t xlo, Double_t xhi,
                Int_t ny, Double_t ylo, Double_t yhi,
                Int_t nz, Double_t zlo, Double_t zhi);

  Int_t Fill (TH1* h, Double_t x, Double_t y, Double_t z) {TH3* h3= dynamic_cast<TH3*>(h); return h3->Fill (x, y, z);}
  static TH1D* Projection3D (const TH1* h, Option_t* xyz, const char* name, const char* title, Option_t* opt);
  static TH1D* ProjectionX (const TH1* h, const char* name=0, const char* title=0, Option_t* opt="") {return Projection3D(h,"x",name,title,opt);}
  static TH1D* ProjectionY (const TH1* h, const char* name=0, const char* title=0, Option_t* opt="") {return Projection3D(h,"y",name,title,opt);}
  static TH1D* ProjectionZ (const TH1* h, const char* name=0, const char* title=0, Option_t* opt="") {return Projection3D(h,"z",name,title,opt);}
};

#endif
