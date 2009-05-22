//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldBinByBin.h,v 1.2 2009-05-22 17:10:20 adye Exp $
//
// Description:
//      Unfold
//
// Author List:
//      Tim Adye <T.J.Adye@rl.ac.uk>
//
// Copyright Information:
//      Copyleft (C) 2006 Rutherford Appleton Laboratory
//
//==============================================================================

#ifndef ROOUNFOLDBINBYBIN_HH
#define ROOUNFOLDBINBYBIN_HH

#include "RooUnfoldBayes.h"

class RooUnfoldResponse;
class TH1;

class RooUnfoldBinByBin : public RooUnfoldBayes {

public:

  // Special constructors
  RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit= false,
                     const char* name= 0, const char* title= 0);

  // Set up an existing object
  virtual RooUnfoldBinByBin& Setup (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit= false);

protected:

  virtual Int_t unfold (vector<Double_t>& causes);
  virtual Int_t train();
  virtual Int_t getCovariance() const;

  // instance variables

public:

  ClassDef (RooUnfoldBinByBin, 0) // Bin-by-bin Unfolding
};

#endif
