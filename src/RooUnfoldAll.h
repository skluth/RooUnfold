//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldAll.h,v 1.5 2010-08-04 16:24:35 fwx38934 Exp $
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
  const RooUnfold* unfold; // Input unfolding object
  RooUnfoldAll (int iterations,const RooUnfold* unfold);
  virtual ~RooUnfoldAll();
  TNtuple* Chi2(const TH1* hTrue=0,Int_t doerror=1); 
  void Plotting();
  TH1* Spread();
  TH1* Unf_err();
  TMatrixD True_err();

protected:
  TH1* h_err; // Output plot
  TH1* h_err_res; // Output plot
  TMatrixD error_matrix; //Matrix of covariances
  TNtuple* chi2;  // Output plot
  TH1* Add_Random(TH1* hMeas);  
  void All_hMeas();
  double xlo; // Minimum x-axis value 
  double xhi; // Maximum x-axis value
  int ntx; // Number of bins in true distribution
  const TH1* hMeas_const; // Measured Distribution
  TH1* hReco; // Reconstructed Distribution
  
  
public:

  ClassDef (RooUnfoldAll, 0)
};
#endif
