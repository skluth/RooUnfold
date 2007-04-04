//==============================================================================
// File and Version Information:
//      $Id: RooUnfold.h,v 1.1.1.1 2007-04-04 21:27:02 adye Exp $
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

#ifndef ROOUNFOLD_HH
#define ROOUNFOLD_HH

#include "TNamed.h"
#include "TVectorD.h"
#include "TMatrixD.h"

class RooUnfoldResponse;
class TH1;

class RooUnfold : public TNamed {

public:

  // Standard methods

  RooUnfold(); // default constructor
  RooUnfold (const char*    name, const char*    title); // named constructor
  RooUnfold (const TString& name, const TString& title); // named constructor
  RooUnfold (const RooUnfold& rhs); // copy constructor
  virtual ~RooUnfold(); // destructor
  virtual RooUnfold& operator= (const RooUnfold& rhs); // assignment operator

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

protected:

  virtual RooUnfold& Setup();
  virtual void SetNameTitleDefault();

  // instance variables

  Int_t _verbose;  // Debug print level
  Int_t _nm;   // Total number of measured bins
  Int_t _nt;   // Total number of truth    bins
  const RooUnfoldResponse* _res;   // Response matrix (not owned)
  const TH1* _meas;                // Measured distribution (not owned)
  TVectorD _rec;  // Reconstructed distribution
  TMatrixD _cov;  // Reconstructed distribution covariance

public:

  ClassDef (RooUnfold, 0) // Unfold
};

// Inline method definitions

inline RooUnfold::RooUnfold()                                           : TNamed()           {Setup();}
inline RooUnfold::RooUnfold (const char*    name, const char*    title) : TNamed(name,title) {Setup();}
inline RooUnfold::RooUnfold (const TString& name, const TString& title) : TNamed(name,title) {Setup();}
inline RooUnfold::~RooUnfold()                                                               {Clear();}

inline Int_t                    RooUnfold::verbose()   const { return _verbose; }
inline const RooUnfoldResponse* RooUnfold::response()  const { return _res;     }
inline const TH1*               RooUnfold::Hmeasured() const { return _meas;    }
inline const TVectorD&          RooUnfold::Vreco()     const { return _rec;     }
inline TVectorD&                RooUnfold::Vreco()           { return _rec;     }
inline const TMatrixD&          RooUnfold::Ereco()     const { return _cov;     }
inline TMatrixD&                RooUnfold::Ereco()           { return _cov;     }

#endif
