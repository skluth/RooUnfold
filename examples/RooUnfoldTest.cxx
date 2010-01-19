//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTest.cxx,v 1.13 2010-01-19 23:30:57 adye Exp $
//
// Description:
//      Tests RooUnfold package using toy MC generated according to PDFs defined
//      in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//      This is the main program. The actual tests are performed using the
//      RooUnfoldTestHarness class.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
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

//==============================================================================
// Routine to run with specified arguments.
// These defaults should probably match those in RooUnfoldTestHarness::Parms.
//==============================================================================

void RooUnfoldTest (
                    Int_t    method=      1,
                    Int_t    stage=       0,
                    Int_t    ftrainx=     2,
                    Int_t    ftestx=      5,
                    Int_t    ntx=        40,
                    Int_t    nmx=        -1,
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

  test->method=  method;
  test->stage=   stage;
  test->ftrainx= ftrainx;
  test->ftestx=  ftestx;
  test->ntx=     ntx;
  test->nmx=     nmx;
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

void RooUnfoldTest (const char* args)
{
  RooUnfoldTestReset();
  test= new RooUnfoldTestHarness ("RooUnfoldTest", args);
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
