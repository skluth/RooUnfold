// File and Version Information:
//      $Id: RooUnfoldTestHarness.h,v 1.1 2010-01-14 01:42:57 adye Exp $
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

#ifndef ROOUNFOLDTESTHARNESS_HH
#define ROOUNFOLDTESTHARNESS_HH

#include "TNamed.h"

class TCanvas;
class TH1;
class TH1D;
class TH2D;
class RooUnfoldResponse;
class RooUnfold;

class RooUnfoldTestHarness : public TNamed {
public:
  // Parameters
  Int_t    method, stage, ftrain, ftest, nt, ntest, ntrain;
  Double_t xlo, xhi;
  Int_t    regparm, ntoys, nm, onepage;

  Bool_t             nosmear;
  Int_t              error, ipad;
  TCanvas*           canvas;
  TH1D               *hPDF, *hTrain, *hTrainTrue, *hTestPDF, *hTrue, *hMeas, *hReco, *hTrue0, *hRes, *hPulls;
  TH2D*              hResmat;
  RooUnfoldResponse* response;
  RooUnfold*         unfold;

  static const Double_t cutdummy= -99999.0;

  // Constructors
  RooUnfoldTestHarness (const char* name= "RooUnfoldTest");
  RooUnfoldTestHarness (const char* name, const char* args);
  RooUnfoldTestHarness (const char* name, int argc, const char* const* argv);
  virtual ~RooUnfoldTestHarness();

  virtual void  Reset();
  virtual void  Defaults();
  virtual Int_t Train();
  virtual Int_t Test();
  virtual void  Unfold();
  virtual Int_t Run();
  virtual int SetArgs (int argc, const char* const* argv, bool split= false);

  void setmax (TH1* h, const TH1* h1= 0, const TH1* h2= 0, const TH1* h3= 0,
                       const TH1* h4= 0, const TH1* h5= 0, const TH1* h6= 0);
  Double_t smear (Double_t xt);
};

#endif
