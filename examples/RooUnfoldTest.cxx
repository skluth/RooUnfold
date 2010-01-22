//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTest.cxx,v 1.17 2010-01-22 15:46:04 adye Exp $
//
// Description:
//      Tests RooUnfold package using toy MC generated according to PDFs defined
//      in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//      This is the main program. The actual tests are performed using the
//      RooUnfoldTestHarness class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldTestHarness.h"

#ifdef __CINT__

RooUnfoldTestHarness* test= 0;
bool RooUnfoldLoaded= false;

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

void RooUnfoldTest (const char* args= "")
{
  if (!(RooUnfoldLoaded++)) gSystem->Load("libRooUnfold");
// If run interactively, remove canvas and all histograms that might have been
// created with a previous invocation.
  delete test; test= 0;
  gDirectory->Clear();

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
