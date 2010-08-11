//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldSvd.h,v 1.14 2010-08-11 19:27:37 adye Exp $
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
  virtual ~RooUnfoldSvd(); // destructor
  RooUnfoldSvd& operator= (const RooUnfoldSvd& rhs); // assignment operator
  virtual RooUnfoldSvd* Clone (const char* newname= 0) const;

  // Special constructors

  RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kterm= 0, Int_t ntoyssvd= 1000,
                const char* name= 0, const char* title= 0);

  void SetKterm (Int_t kterm) { _kterm= kterm; }
  void SetNtoysSVD (Int_t ntoyssvd) { _ntoyssvd= ntoyssvd; }
  Int_t GetKterm() const { return _kterm; }
  Int_t GetNtoysSVD() const { return _ntoyssvd; }

  virtual void  SetRegParm (Double_t parm) { SetKterm(Int_t(parm+0.5)); }
  virtual Double_t GetRegParm() const { return GetKterm(); }
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
  virtual void GetSettings();

  // instance variables
  TUnfHisto* _svd;
  Int_t _kterm;
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
