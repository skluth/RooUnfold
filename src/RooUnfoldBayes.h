//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBayes.h,v 1.6 2010-01-26 00:53:17 adye Exp $
//
// Description:
//      Bayesian unfolding. Just an interface to RooUnfoldBayesImpl.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
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
  virtual ~RooUnfoldBayes(); // destructor
  RooUnfoldBayes& operator= (const RooUnfoldBayes& rhs); // assignment operator

  // Special constructors
  RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas, Int_t niter= 4, Bool_t smoothit= false,
                  const char* name= 0, const char* title= 0);

  virtual void SetIterations (Int_t niter= 4)         { _niter=    niter;    }
  virtual void SetSmoothing  (Bool_t smoothit= false) { _smoothit= smoothit; }

  virtual void Reset();
  virtual void Print (Option_t* option= "") const;

  static vector<Double_t>& H2VD (const TH1*  h, vector<Double_t>& v);
  static Array2D&          H2AD (const TH2D* h, Array2D& m, const TH1* norm= 0);
  static TVectorD&         VD2V (const vector<Double_t>& vd, TVectorD& v);
  static TMatrixD&         AD2M (const Array2D& ad, TMatrixD& m);

protected:

  void Init();
  void Destroy();
  virtual void Unfold();
  virtual Int_t unfold (vector<Double_t>& causes);
  virtual Int_t getCovariance() const;
  virtual void GetCov();
  void Assign   (const RooUnfoldBayes& rhs); // implementation of assignment operator
  void CopyData (const RooUnfoldBayes& rhs);

  // instance variables
  RooUnfoldBayesImpl* _bayes;
  Int_t _niter;
  Int_t _smoothit;

private:

public:

  ClassDef (RooUnfoldBayes, 0) // Bayesian Unfolding
};

// Inline method definitions

inline RooUnfoldBayes::RooUnfoldBayes()                                           : RooUnfold()           {Init();}
inline RooUnfoldBayes::RooUnfoldBayes (const char* name, const char* title)       : RooUnfold(name,title) {Init();}
inline RooUnfoldBayes::RooUnfoldBayes (const TString& name, const TString& title) : RooUnfold(name,title) {Init();}
inline RooUnfoldBayes::~RooUnfoldBayes() {Destroy();}
inline RooUnfoldBayes& RooUnfoldBayes::operator= (const RooUnfoldBayes& rhs) {Assign(rhs); return *this;}

#endif
