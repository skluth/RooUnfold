//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBinByBinOld.h,v 1.1 2010-08-18 12:58:04 fwx38934 Exp $
//
// Description:
//      Unfolding bin-by-bin. Just an interface to RooUnfoldBayesImpl.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDBINBYBINOLD_HH
#define ROOUNFOLDBINBYBINOLD_HH

#include "RooUnfoldBayes.h"

class RooUnfoldResponse;
class TH1;

class RooUnfoldBinByBinOld : public RooUnfoldBayes {

public:

  // Standard methods

  RooUnfoldBinByBinOld(); // default constructor
  RooUnfoldBinByBinOld (const char*    name, const char*    title); // named constructor
  RooUnfoldBinByBinOld (const TString& name, const TString& title); // named constructor
  RooUnfoldBinByBinOld (const RooUnfoldBinByBinOld& rhs); // copy constructor
  RooUnfoldBinByBinOld& operator= (const RooUnfoldBinByBinOld& rhs); // assignment operator

  // Special constructors
  RooUnfoldBinByBinOld (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit= false,
                     const char* name= 0, const char* title= 0);

protected:

  void Init();
  virtual Int_t unfold (vector<Double_t>& causes);
  virtual Int_t getCovariance() const;
  virtual void GetSettings();

  // instance variables

public:

  ClassDef (RooUnfoldBinByBinOld, 0) // Bin-by-bin Unfolding
};

// Inline method definitions

inline void RooUnfoldBinByBinOld::Init() {GetSettings();}
inline RooUnfoldBinByBinOld::RooUnfoldBinByBinOld()                                           : RooUnfoldBayes()           {Init();}
inline RooUnfoldBinByBinOld::RooUnfoldBinByBinOld (const char* name, const char* title)       : RooUnfoldBayes(name,title) {Init();}
inline RooUnfoldBinByBinOld::RooUnfoldBinByBinOld (const TString& name, const TString& title) : RooUnfoldBayes(name,title) {Init();}
inline RooUnfoldBinByBinOld& RooUnfoldBinByBinOld::operator= (const RooUnfoldBinByBinOld& rhs) {Assign(rhs); return *this;}

#endif
