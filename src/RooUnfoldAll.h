//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldAll.h,v 1.7 2010-08-10 14:19:09 fwx38934 Exp $
//
// Description:
//      Graph Drawing Class for use with RooUnfold.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Richard Claridge <richard.claridge@stfc.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDALL_H_
#define ROOUNFOLDALL_H_
#include "TNamed.h"
class TH1;
class RooUnfold;
class TNtuple;
#include "TMatrixD.h"

class RooUnfoldAll : public TNamed {

public:
  
  int iterations; // Number of iterations 
  RooUnfold* unfold; // Input unfolding object
  const TH1* hTrue;
  RooUnfoldAll (int iterations,RooUnfold* unfold,const TH1* Truth);
  virtual ~RooUnfoldAll();
  TNtuple* Chi2(); 

  TH1* Spread();
  TH1* Unf_err();

protected:
  void Plotting();
  TH1* h_err; // Output plot
  TH1* h_err_res; // Output plot
  TNtuple* hchi2;  // Output plot 
  void All_hMeas(); //
  double xlo; // Minimum x-axis value 
  double xhi; // Maximum x-axis value
  int ntx; // Number of bins in true distribution
  const TH1* hMeas_const; // Measured Distribution
  
  
public:

  ClassDef (RooUnfoldAll, 0)
};
#endif
