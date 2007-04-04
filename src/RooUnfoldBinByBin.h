//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldBinByBin.h,v 1.1.1.1 2007-04-04 21:27:02 adye Exp $
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

  // instance variables

public:

  ClassDef (RooUnfoldBinByBin, 0) // Bin-by-bin Unfolding
};

#endif
