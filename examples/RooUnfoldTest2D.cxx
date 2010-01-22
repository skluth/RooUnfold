//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTest2D.cxx,v 1.17 2010-01-22 15:46:04 adye Exp $
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

#include "RooUnfoldTestHarness2D.h"

#ifdef __CINT__

RooUnfoldTestHarness2D* test2d= 0;
bool RooUnfoldLoaded= false;

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

void RooUnfoldTest2D (const char* args= "")
{
  if (!(RooUnfoldLoaded++)) gSystem->Load("libRooUnfold");
// If run interactively, remove canvas and all histograms that might have been
// created with a previous invocation.
  delete test2d; test2d= 0;
  gDirectory->Clear();

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
