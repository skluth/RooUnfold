//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldBayes.h,v 1.1.1.1 2007-04-04 21:27:02 adye Exp $
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

#ifndef ROOUNFOLDBAYES_HH
#define ROOUNFOLDBAYES_HH

#include "RooUnfold.h"
#include <vector>
using std::vector;

class RooUnfoldResponse;
class TH1;
class TH2D;
class RooUnfoldBayesImpl;
class Array2D;

class RooUnfoldBayes : public RooUnfold {

public:

  // Standard methods

  RooUnfoldBayes(); // default constructor
  RooUnfoldBayes (const char*    name, const char*    title); // named constructor
  RooUnfoldBayes (const TString& name, const TString& title); // named constructor
  RooUnfoldBayes (const RooUnfoldBayes& rhs); // copy constructor

  // Special constructors
  RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas, Int_t niter= 4, Bool_t smoothit= false,
                     const char* name= 0, const char* title= 0);

  // Set up an existing object
  virtual RooUnfoldBayes& Clear ();
  virtual RooUnfoldBayes& Setup (const RooUnfoldBayes& rhs);
  virtual RooUnfoldBayes& Setup (const RooUnfoldResponse* res, const TH1* meas, Int_t niter= 4, Bool_t smoothit= false);

  virtual void Print (Option_t* o= 0) const;

  static vector<Double_t>& H2VD (const TH1*  h, vector<Double_t>& v);
  static Array2D&          H2AD (const TH2D* h, Array2D& m, const TH1* norm= 0);
  static TVectorD&         VD2V (const vector<Double_t>& vd, TVectorD& v);
  static TMatrixD&         AD2M (const Array2D& ad, TMatrixD& m);

protected:

  virtual RooUnfoldBayes& Setup();

  // instance variables
  RooUnfoldBayesImpl* _bayes;
  Int_t _niter;
  Int_t _smoothit;

public:

  ClassDef (RooUnfoldBayes, 0) // Bayesian Unfolding
};

// Inline method definitions

inline RooUnfoldBayes::RooUnfoldBayes()                                           : RooUnfold()           {Setup();}
inline RooUnfoldBayes::RooUnfoldBayes (const char* name, const char* title)       : RooUnfold(name,title) {Setup();}
inline RooUnfoldBayes::RooUnfoldBayes (const TString& name, const TString& title) : RooUnfold(name,title) {Setup();}

#endif
