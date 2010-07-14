//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldSvd.h,v 1.8 2010-07-14 21:57:44 adye Exp $
//
// Description:
//      SVD unfolding. Just an interface to RooUnfHistoSvd.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDSVD_HH
#define ROOUNFOLDSVD_HH

#include "RooUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;
class TUnfHisto;

class RooUnfoldSvd : public RooUnfold {

public:

  // Standard methods

  RooUnfoldSvd(); // default constructor
  RooUnfoldSvd (const char*    name, const char*    title); // named constructor
  RooUnfoldSvd (const TString& name, const TString& title); // named constructor
  RooUnfoldSvd (const RooUnfoldSvd& rhs); // copy constructor
  RooUnfoldSvd& operator= (const RooUnfoldSvd& rhs); // assignment operator

  // Special constructors
  RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kterm= 1, Int_t ntoys= 1000,
                const char* name= 0, const char* title= 0);

  void SetKterm (Int_t kterm) { _kterm= kterm; }
  void SetNtoys (Int_t ntoys) { _ntoys= ntoys; }
  Int_t GetKterm() const { return _kterm; }
  Int_t GetNtoys() const { return _ntoys; }

  virtual void  SetRegParm (Int_t parm) { SetKterm(parm); }
  virtual Int_t GetRegParm() const { return GetKterm(); }

  virtual void Reset();
  virtual TObject* Impl();

protected:

  void Init();
  void Destroy();
  virtual void Unfold();
  virtual void GetCov();
  void Assign   (const RooUnfoldSvd& rhs); // implementation of assignment operator
  void CopyData (const RooUnfoldSvd& rhs);
  TH1D* CopyOverflow   (const TH1* h) const;
  TH2D* CopyOverflow2D (const TH1* h) const;

  // instance variables
  TUnfHisto* _svd;
  Int_t _kterm;
  Int_t _ntoys;

  TH1D *_meas1d, *_train1d, *_truth1d;
  TH2D *_reshist;

public:

  ClassDef (RooUnfoldSvd, 0) // SVD Unfolding
};

// Inline method definitions

inline RooUnfoldSvd::RooUnfoldSvd()                                           : RooUnfold()           {Init();}
inline RooUnfoldSvd::RooUnfoldSvd (const char* name, const char* title)       : RooUnfold(name,title) {Init();}
inline RooUnfoldSvd::RooUnfoldSvd (const TString& name, const TString& title) : RooUnfold(name,title) {Init();}
inline RooUnfoldSvd& RooUnfoldSvd::operator= (const RooUnfoldSvd& rhs) {Assign(rhs); return *this;}

#endif
