//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
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
class TSVDUnfold;

class RooUnfoldSvd : public RooUnfold {

public:

  // Standard methods

  RooUnfoldSvd(); // default constructor
  RooUnfoldSvd (const char*    name, const char*    title); // named constructor
  RooUnfoldSvd (const TString& name, const TString& title); // named constructor
  RooUnfoldSvd (const RooUnfoldSvd& rhs); // copy constructor
  virtual ~RooUnfoldSvd(); // destructor
  RooUnfoldSvd& operator= (const RooUnfoldSvd& rhs); // assignment operator
  virtual RooUnfoldSvd* Clone (const char* newname= 0) const;

  // Special constructors

  RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kreg= 0, Int_t ntoyssvd= 1000,
                const char* name= 0, const char* title= 0);

  void SetKterm (Int_t kreg) { _kreg= kreg; }
  void SetNtoysSVD (Int_t ntoyssvd) { _ntoyssvd= ntoyssvd; }
  Int_t GetKterm() const { return _kreg; }
  Int_t GetNtoysSVD() const { return _ntoyssvd; }

  virtual void  SetRegParm (Double_t parm) { SetKterm(Int_t(parm+0.5)); }
  virtual Double_t GetRegParm() const { return GetKterm(); }
  virtual void Reset();
  TSVDUnfold* Impl();

protected:
  void Assign (const RooUnfoldSvd& rhs); // implementation of assignment operator
  virtual void Unfold();
  virtual void GetCov();
  virtual void GetSettings();

private:
  void Init();
  void Destroy();
  void CopyData (const RooUnfoldSvd& rhs);

protected:
  // instance variables
  TSVDUnfold* _svd;
  Int_t _kreg;
  Int_t _ntoyssvd;

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
inline RooUnfoldSvd::~RooUnfoldSvd() {Destroy();}

#endif
