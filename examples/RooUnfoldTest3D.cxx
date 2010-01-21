//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTest3D.cxx,v 1.8 2010-01-21 20:05:14 adye Exp $
//
// Description:
//      3D test of RooUnfold package using toy MC generated according to PDFs
//      defined in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//      This is the main program. The actual tests are performed using the
//      RooUnfoldTestHarness3D class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldTestHarness3D.h"

#ifdef __CINT__

RooUnfoldTestHarness3D* test3d= 0;
bool RooUnfoldLoaded3D= false;

void RooUnfoldTestReset3D()
{
  if (!RooUnfoldLoaded3D) {
    gSystem->Load("libRooUnfold");
    RooUnfoldLoaded3D= true;
  }
// If run interactively, remove canvas and all histograms that might have been
// created with a previous invocation.
  delete test3d; test3d= 0;
  gDirectory->Clear();
}

//==============================================================================
// Routine to run with specified arguments.
// These defaults should probably match those in RooUnfoldTestHarness3D::Parms.
//==============================================================================

void RooUnfoldTest3D (
                      Int_t    method=      1,
                      Int_t    stage=       0
                     )
{
  RooUnfoldTestReset3D();
  test3d= new RooUnfoldTestHarness3D ("RooUnfoldTest3D");

  test3d->method=  method;
  test3d->stage=   stage;

  test3d->Run();
}

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

void RooUnfoldTest3D (const char* args)
{
  RooUnfoldTestReset3D();
  test3d= new RooUnfoldTestHarness3D ("RooUnfoldTest3D", args);
  test3d->Run();
}

#else   // __CINT__

//==============================================================================
// Main program when run stand-alone
//==============================================================================

int main (int argc, char** argv) {
  RooUnfoldTestHarness3D test3d ("RooUnfoldTest3D", argc, argv);
  return test3d.Run();
}

#endif
