//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
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
class TH2;
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
  virtual RooUnfoldBayes* Clone (const char* newname= 0) const;

  // Special constructors

  RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas, Int_t niter= 4, Bool_t smoothit= false,
                  const char* name= 0, const char* title= 0);

  void SetIterations (Int_t niter= 4)         { _niter=    niter;    }
  void SetSmoothing  (Bool_t smoothit= false) { _smoothit= smoothit; }
  Int_t GetIterations() const { return _niter;    }
  Int_t GetSmoothing()  const { return _smoothit; }

  virtual void  SetRegParm (Double_t parm) { SetIterations(Int_t(parm+0.5)); }
  virtual Double_t GetRegParm() const { return GetIterations(); }
  virtual void Reset();
  virtual void Print (Option_t* option= "") const;
  RooUnfoldBayesImpl* Impl();

  static vector<Double_t>& H2VD (const TH1*  h, vector<Double_t>& v, Bool_t overflow= kFALSE);
  static Array2D&          H2AD (const TH2* h, Array2D& m, const TH1* norm= 0, Bool_t overflow= kFALSE);
  static TVectorD&         VD2V (const vector<Double_t>& vd, TVectorD& v);
  static TMatrixD&         AD2M (const Array2D& ad, TMatrixD& m);
  static TVectorD&         AD2V (const Array2D& ad, TVectorD& m);

protected:
  void Assign (const RooUnfoldBayes& rhs); // implementation of assignment operator
  virtual void Unfold();
  virtual void GetCov();
  virtual void GetErrors();
  virtual void GetSettings();

private:
  void Init();
  void CopyData (const RooUnfoldBayes& rhs);

protected:
  // instance variables
  RooUnfoldBayesImpl* _bayes;
  Int_t _niter;
  Int_t _smoothit;

public:
  ClassDef (RooUnfoldBayes, 0) // Bayesian Unfolding
};

// Inline method definitions

inline RooUnfoldBayes::RooUnfoldBayes()                                           : RooUnfold()           {Init();}
inline RooUnfoldBayes::RooUnfoldBayes (const char* name, const char* title)       : RooUnfold(name,title) {Init();}
inline RooUnfoldBayes::RooUnfoldBayes (const TString& name, const TString& title) : RooUnfold(name,title) {Init();}
inline RooUnfoldBayes& RooUnfoldBayes::operator= (const RooUnfoldBayes& rhs) {Assign(rhs); return *this;}

#endif
