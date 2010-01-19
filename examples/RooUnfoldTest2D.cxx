//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTest2D.cxx,v 1.13 2010-01-19 23:30:57 adye Exp $
//
// Description:
//      2D test of RooUnfold package using toy MC generated according to PDFs
//      defined in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//      This is the main program. The actual tests are performed using the
//      RooUnfoldTestHarness2D class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldTestHarness2D.icc"

#ifdef __CINT__

RooUnfoldTestHarness2D* test2d= 0;
bool RooUnfoldLoaded2D= false;

void RooUnfoldTestReset2D()
{
  if (!RooUnfoldLoaded2D) {
    gSystem->Load("libRooUnfold");
    RooUnfoldLoaded2D= true;
  }
// If run interactively, remove canvas and all histograms that might have been
// created with a previous invocation.
  delete test2d; test2d= 0;
  gDirectory->Clear();
}

//==============================================================================
// Routine to run with specified arguments.
// These defaults should probably match those in RooUnfoldTestHarness::Parms.
//==============================================================================

void RooUnfoldTest2D (
                      Int_t    method=      1,
                      Int_t    stage=       0,
                      Int_t    ftrainx=     1,
                      Int_t    ftestx=      3,
                      Int_t    ntx=        40,
                      Int_t    nmx=        -1,
                      Int_t    ntest=   10000,
                      Int_t    ntrain= 100000,
                      Double_t xlo=     -12.5,
                      Double_t xhi=      10.0,
                      Int_t    regparm=  -999,  // Bayes niter=4, SVD kterm=20
                      Int_t    ntoys=    1000     // SVD only
                   )
{
  RooUnfoldTestReset2D();
  test2d= new RooUnfoldTestHarness2D ("RooUnfoldTest2D");

  test2d->method=  method;
  test2d->stage=   stage;
  test2d->ftrainx= ftrainx;
  test2d->ftestx=  ftestx;
  test2d->ntx=     ntx;
  test2d->nmx=     nmx;
  test2d->ntest=   ntest;
  test2d->ntrain=  ntrain;
  test2d->xlo=     xlo;
  test2d->xhi=     xhi;
  test2d->regparm= regparm;
  test2d->ntoys=   ntoys;

  test2d->Run();
}

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

void RooUnfoldTest2D (const char* args)
{
  RooUnfoldTestReset2D();
  test2d= new RooUnfoldTestHarness2D ("RooUnfoldTest2D", args);
  test2d->Run();
}

#else   // __CINT__

//==============================================================================
// Main program when run stand-alone
//==============================================================================

int main (int argc, char** argv) {
  RooUnfoldTestHarness2D test2d ("RooUnfoldTest2D", argc, argv);
  return test2d.Run();
}

#endif
