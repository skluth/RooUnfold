//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTestHarness.h,v 1.9 2010-01-20 20:36:25 adye Exp $
//
// Description:
//      Harness class to test the RooUnfold package using toy MC generated
//      according to PDFs defined in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#ifndef ROOUNFOLDTESTHARNESS_HH
#define ROOUNFOLDTESTHARNESS_HH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TNamed.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
#include "TVectorDfwd.h"
#else
class TVectorD;
#endif
#endif

#ifdef __CINT__
#include "ArgVars.h"
#else
class ArgVars;
#endif

class TCanvas;
class TH1;
class TH1D;
class TH2D;
class RooUnfoldResponse;
class RooUnfold;

class RooUnfoldTestHarness : public TNamed {
public:
  // Parameters
  Int_t    method, stage, ftrainx, ftestx, ntx, ntest, ntrain;
  Double_t xlo, xhi, bkgtx, bkgex, efflo, effhi, bias, smear;
  Int_t    regparm, ntoys, nmx, onepage, doerror, dim;

  Bool_t             nosmear;
  Int_t              error, ipad, ntbins, nmbins;
  TCanvas*           canvas;
  TH1                *hTrain, *hTrainTrue, *hTrue, *hMeas, *hReco, *hRes, *hPulls;
  TH2D*              hResmat;
  RooUnfoldResponse* response;
  RooUnfold*         unfold;

  static const Int_t    nbPDF=         500;
  static const Double_t cutdummy= -99999.0;

  TH1D               *hPDFx, *hTestPDFx;

  // Constructors
  RooUnfoldTestHarness (const char* name= "RooUnfoldTest");
  RooUnfoldTestHarness (const char* name, const char* args);
  RooUnfoldTestHarness (const char* name, int argc, const char* const* argv);
  virtual ~RooUnfoldTestHarness();

  TH1D* Generate (TVectorD& x, const char* name, const char* title, Int_t nt, Int_t fpdf, Int_t nx, Double_t xlo, Double_t xhi,
                  Double_t bkg= 0.0, Double_t mean= 0.0, Double_t width= 2.5);
  virtual void  Reset();
  virtual void  Init();
  virtual void  Init2();
  virtual Int_t Train();
  virtual Int_t Test();
  virtual Int_t Unfold();
  virtual void  Results();
  virtual Int_t Run();
  virtual Int_t RunStuff();
  virtual void  Print      (std::ostream& o)                       const;
  virtual void  PrintParms (std::ostream& o, const char* sep= " ") const;
  virtual Int_t CheckParms();
  virtual void  Parms (ArgVars& args);
  virtual int   SetArgs (int argc, const char* const* argv, bool split= false);
  virtual void  SetDefaults();

  void setmax (TH1* h, const TH1* h1= 0, const TH1* h2= 0, const TH1* h3= 0,
                       const TH1* h4= 0, const TH1* h5= 0, const TH1* h6= 0);
  Double_t Smear (Double_t xt, Int_t nt, Double_t xlo, Double_t xhi);
};

#ifndef NOINLINE
#include "RooUnfoldTestHarness.icc"
#endif

#endif
