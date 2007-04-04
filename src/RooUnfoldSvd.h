//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldSvd.h,v 1.1.1.1 2007-04-04 21:27:02 adye Exp $
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

#ifndef ROOUNFOLDSVD_HH
#define ROOUNFOLDSVD_HH

#include "RooUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TUnfHisto;

class RooUnfoldSvd : public RooUnfold {

public:

  // Standard methods

  RooUnfoldSvd(); // default constructor
  RooUnfoldSvd (const char*    name, const char*    title); // named constructor
  RooUnfoldSvd (const TString& name, const TString& title); // named constructor
  RooUnfoldSvd (const RooUnfoldSvd& rhs); // copy constructor

  // Special constructors
  RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t tau= 1, Int_t ntoys= 1000,
                const char* name= 0, const char* title= 0);

  // Set up an existing object
  virtual RooUnfoldSvd& Clear ();
  virtual RooUnfoldSvd& Setup (const RooUnfoldSvd& rhs);
  virtual RooUnfoldSvd& Setup (const RooUnfoldResponse* res, const TH1* meas, Int_t tau= 1, Int_t ntoys= 1000);

protected:

  virtual RooUnfoldSvd& Setup();

  // instance variables
  TUnfHisto* _svd;
  Int_t _tau;
  Int_t _ntoys;

  TH1D *_meas1d, *_train1d, *_truth1d;

public:

  ClassDef (RooUnfoldSvd, 0) // SVD Unfolding
};

// Inline method definitions

inline RooUnfoldSvd::RooUnfoldSvd()                                           : RooUnfold()           {Setup();}
inline RooUnfoldSvd::RooUnfoldSvd (const char* name, const char* title)       : RooUnfold(name,title) {Setup();}
inline RooUnfoldSvd::RooUnfoldSvd (const TString& name, const TString& title) : RooUnfold(name,title) {Setup();}

#endif
