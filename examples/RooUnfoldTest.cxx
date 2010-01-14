// File and Version Information:
//      $Id: RooUnfoldTest.cxx,v 1.10 2010-01-14 01:42:57 adye Exp $
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

#include "RooUnfoldTestHarness.icc"

#ifdef __CINT__

RooUnfoldTestHarness* test= 0;
bool RooUnfoldLoaded= false;

void RooUnfoldTestReset()
{
  if (!RooUnfoldLoaded) {
    gSystem->Load("libRooUnfold");
    RooUnfoldLoaded= true;
  }
// If run interactively, remove canvas and all histograms that might have been
// created with a previous invocation.
  delete test; test= 0;
  gDirectory->Clear();
}

RooUnfoldTestError()
{
  Int_t error= test->error;
  if (error) {
    delete test; test= 0;
  }
  return error;
}

//==============================================================================
// Routine to run with specified arguments.
// These defaults should probably match those in RooUnfoldTestHarness::Defaults.
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
                    Int_t    regparm=  -999,  // Bayes niter=4, SVD kterm=20
                    Int_t    ntoys=    1000   // SVD only
                   )
{
  RooUnfoldTestReset();
  test= new RooUnfoldTestHarness ("RooUnfoldTest");
  if (RooUnfoldTestError()) return;

  test->method=  method;
  test->stage=   stage;
  test->ftrain=  ftrain;
  test->ftest=   ftest;
  test->nt=      nt;
  test->ntest=   ntest;
  test->ntrain=  ntrain;
  test->xlo=     xlo;
  test->xhi=     xhi;
  test->regparm= regparm;
  test->ntoys=   ntoys;

  test->Run();
}

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

RooUnfoldTest (const char* args)
{
  RooUnfoldTestReset();
  test= new RooUnfoldTestHarness ("RooUnfoldTest", args);
  if (RooUnfoldTestError()) return;
  test->Run();
}

#else   // __CINT__

//==============================================================================
// Main program when run stand-alone
//==============================================================================

int main (int argc, char** argv) {
  RooUnfoldTestHarness test ("RooUnfoldTest", argc, argv);
  return test.Run();
}

#endif
