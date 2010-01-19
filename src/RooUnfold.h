//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.h,v 1.6 2010-01-19 15:33:47 adye Exp $
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

  virtual RooUnfold& Clear ();
  virtual RooUnfold& Setup (const RooUnfold& rhs);
  virtual RooUnfold& Setup (const RooUnfoldResponse* res, const TH1* meas);

  // Accessors

  virtual const RooUnfoldResponse* response() const;
  virtual const TH1*               Hmeasured() const;
  virtual TH1*                     Hreco (Bool_t withError= true) const;

  virtual const TVectorD&          Vreco() const;
  virtual TVectorD&                Vreco();
  virtual const TMatrixD&          Ereco() const;
  virtual TMatrixD&                Ereco();

  virtual Int_t                    verbose() const;

  virtual void PrintTable (std::ostream& o, const TH1* hTrue= 0, Bool_t withError= true) const;

protected:

  virtual RooUnfold& Setup();
  virtual void SetNameTitleDefault();
  virtual void GetCov() const;  // actually updates mutable _cov
  virtual void Assign (const RooUnfold& rhs); // implementation of assignment operator

  // instance variables

  Int_t _verbose;  // Debug print level
  Int_t _nm;   // Total number of measured bins
  Int_t _nt;   // Total number of truth    bins
  mutable Bool_t _haveCov;
  const RooUnfoldResponse* _res;   // Response matrix (not owned)
  const TH1* _meas;                // Measured distribution (not owned)
  TVectorD _rec;  // Reconstructed distribution
  mutable TMatrixD _cov;  // Reconstructed distribution covariance

public:

  ClassDef (RooUnfold, 0) // Unfold
};

// Inline method definitions

inline RooUnfold::RooUnfold()                                           : TNamed()           {Setup();}
inline RooUnfold::RooUnfold (const char*    name, const char*    title) : TNamed(name,title) {Setup();}
inline RooUnfold::RooUnfold (const TString& name, const TString& title) : TNamed(name,title) {Setup();}
inline RooUnfold::~RooUnfold()                                                               {Clear();}
inline RooUnfold& RooUnfold::operator= (const RooUnfold& rhs) {Assign(rhs); return *this;}

inline Int_t                    RooUnfold::verbose()   const { return _verbose; }
inline const RooUnfoldResponse* RooUnfold::response()  const { return _res;     }
inline const TH1*               RooUnfold::Hmeasured() const { return _meas;    }
inline const TVectorD&          RooUnfold::Vreco()     const { return _rec;     }
inline TVectorD&                RooUnfold::Vreco()           { return _rec;     }
inline const TMatrixD&          RooUnfold::Ereco()     const { if (!_haveCov) GetCov(); return _cov; }
inline TMatrixD&                RooUnfold::Ereco()           { if (!_haveCov) GetCov(); return _cov; }

#endif
