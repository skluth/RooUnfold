//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.h,v 1.8 2010-05-20 22:50:59 adye Exp $
//
// Description:
//      Unfolding framework base class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLD_HH
#define ROOUNFOLD_HH

#include "TNamed.h"
#include "TVectorD.h"
#include "TMatrixD.h"

class TH1;
class RooUnfoldResponse;

class RooUnfold : public TNamed {

public:

  // Standard methods

  RooUnfold(); // default constructor
  RooUnfold (const char*    name, const char*    title); // named constructor
  RooUnfold (const TString& name, const TString& title); // named constructor
  RooUnfold (const RooUnfold& rhs); // copy constructor
  virtual ~RooUnfold(); // destructor
  RooUnfold& operator= (const RooUnfold& rhs); // assignment operator

  // Special constructors

  RooUnfold (const RooUnfoldResponse* res, const TH1* meas, const char* name= 0, const char* title= 0);

  // Set up an existing object

  virtual RooUnfold& Setup (const RooUnfoldResponse* res, const TH1* meas);
  virtual void Reset ();

  // Accessors

  virtual const RooUnfoldResponse* response() const;
  virtual const TH1*               Hmeasured() const;
  virtual TH1*                     Hreco (Bool_t withError= true);

  virtual TVectorD&                Vreco();
  virtual TMatrixD&                Ereco();

  virtual Int_t                    verbose() const;
  virtual void SetVerbose (Int_t level);

  virtual void PrintTable (std::ostream& o, const TH1* hTrue= 0, Bool_t withError= true);

  virtual TObject* Impl();


protected:

  void Init();
  virtual void Unfold();
  virtual void GetCov();
  virtual void SetNameTitleDefault();
  void Assign   (const RooUnfold& rhs); // implementation of assignment operator
  void CopyData (const RooUnfold& rhs);

  // instance variables

  Int_t _verbose;  // Debug print level
  Int_t _nm;   // Total number of measured bins
  Int_t _nt;   // Total number of truth    bins
  mutable Bool_t _unfolded, _haveCov, _fail;
  const RooUnfoldResponse* _res;   // Response matrix (not owned)
  const TH1* _meas;                // Measured distribution (not owned)
  mutable TVectorD _rec;  // Reconstructed distribution
  mutable TMatrixD _cov;  // Reconstructed distribution covariance

private:

public:

  ClassDef (RooUnfold, 0) // Unfold
};

// Inline method definitions

inline RooUnfold::RooUnfold()                                           : TNamed()           {Init();}
inline RooUnfold::RooUnfold (const char*    name, const char*    title) : TNamed(name,title) {Init();}
inline RooUnfold::RooUnfold (const TString& name, const TString& title) : TNamed(name,title) {Init();}
inline RooUnfold::~RooUnfold() {}
inline RooUnfold& RooUnfold::operator= (const RooUnfold& rhs) {Assign(rhs); return *this;}

inline Int_t                    RooUnfold::verbose()   const { return _verbose; }
inline const RooUnfoldResponse* RooUnfold::response()  const { return _res;     }
inline const TH1*               RooUnfold::Hmeasured() const { return _meas;    }
inline TVectorD&                RooUnfold::Vreco()           { if (!_unfolded) Unfold(); return _rec; }
inline TMatrixD&                RooUnfold::Ereco()           { if (!_haveCov)  GetCov(); return _cov; }
inline TObject*                 RooUnfold::Impl()            { return 0; };
inline void RooUnfold::SetVerbose (Int_t level)              { _verbose= level; }

#endif
