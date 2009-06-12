//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldBinByBin.h,v 1.3 2009-06-12 00:44:42 adye Exp $
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

  // Standard methods

  RooUnfoldBinByBin(); // default constructor
  RooUnfoldBinByBin (const char*    name, const char*    title); // named constructor
  RooUnfoldBinByBin (const TString& name, const TString& title); // named constructor
  RooUnfoldBinByBin (const RooUnfoldBinByBin& rhs); // copy constructor
  RooUnfoldBinByBin& operator= (const RooUnfoldBinByBin& rhs); // assignment operator

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

// Inline method definitions

inline RooUnfoldBinByBin::RooUnfoldBinByBin()                                           : RooUnfoldBayes()           {}
inline RooUnfoldBinByBin::RooUnfoldBinByBin (const char* name, const char* title)       : RooUnfoldBayes(name,title) {}
inline RooUnfoldBinByBin::RooUnfoldBinByBin (const TString& name, const TString& title) : RooUnfoldBayes(name,title) {}
inline RooUnfoldBinByBin& RooUnfoldBinByBin::operator= (const RooUnfoldBinByBin& rhs) {Assign(rhs); return *this;}

#endif
