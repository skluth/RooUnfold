//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldAll.h,v 1.2 2010-07-22 10:30:25 fwx38934 Exp $
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
class TLegend;
class TH1;
class RooUnfold;
class TNtuple;

class RooUnfoldAll : public TNamed {

public:
  //Data
  double xlo, xhi;
  int ntx, iterations, nmx;
  const TH1* hTrue;
  const TH1* hMeas_const;
  TH1* hReco;
  
  TH1* h_err;
  TH1* h_err_res;
  TNtuple* chi2;  

  
  const RooUnfold* unfold;

  // Constructor
  RooUnfoldAll (int iterations,const TH1* hTrue,const RooUnfold* unfold);
  virtual ~RooUnfoldAll();

  // Methods and functions
  void All_hTrue();
  void All_hMeas();
  TNtuple* Chi2();
  void Plotting();
  TH1* Spread();
  TH1* Add_Random(TH1* hMeas);
  TH1* Unf_err();


  
  
public:

  ClassDef (RooUnfoldAll, 0)
};
#endif
