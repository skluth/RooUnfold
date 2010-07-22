//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.h,v 1.14 2010-07-22 16:15:30 fwx38934 Exp $
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

  enum Algorithm { kNone, kBayes, kSVD, kBinByBin };

  static RooUnfold* New (Algorithm alg, const RooUnfoldResponse* res, const TH1* meas, Int_t regparm= -1,
                         const char* name= 0, const char* title= 0);

  // Standard methods

  RooUnfold(); // default constructor
  RooUnfold (const char*    name, const char*    title); // named constructor
  RooUnfold (const TString& name, const TString& title); // named constructor
  RooUnfold (const RooUnfold& rhs); // copy constructor
  virtual ~RooUnfold(); // destructor
  RooUnfold& operator= (const RooUnfold& rhs); // assignment operator
  virtual RooUnfold* Clone (const char* newname= 0) const;

  // Special constructors

  RooUnfold (const RooUnfoldResponse* res, const TH1* meas, const char* name= 0, const char* title= 0);

  // Set up an existing object

  virtual RooUnfold& Setup (const RooUnfoldResponse* res, const TH1* meas);
  virtual void Reset ();

  // Accessors

  virtual const RooUnfoldResponse* response() const;
  virtual const TH1*               Hmeasured() const;
  virtual TH1*                     Hreco (Int_t withError= 1);

  virtual TVectorD&                Vreco();
  virtual TMatrixD&                Ereco();
  virtual TMatrixD&                Freco();

  virtual Int_t                    verbose() const;
  virtual void SetVerbose (Int_t level);
  virtual Int_t                    Nits() const;
  virtual void SetNits (Int_t iterations);

  virtual void PrintTable (std::ostream& o, const TH1* hTrue= 0, Int_t withError= 1);

  virtual TObject* Impl();

  virtual void  SetRegParm (Int_t parm);
  virtual Int_t GetRegParm() const;

  Double_t Chi2 (const TH1* hTrue);
protected:

  void Init();
  virtual void Unfold();
  virtual void GetCov();
  virtual void SetNameTitleDefault();
  virtual void Get_err_mat();
  void Assign   (const RooUnfold& rhs); // implementation of assignment operator
  void CopyData (const RooUnfold& rhs);

  // instance variables

  Int_t _verbose;  // Debug print level
  Int_t _nm;   // Total number of measured bins
  Int_t _nt;   // Total number of truth    bins
  Int_t _overflow;   // Use histogram under/overflows if 1 (set from RooUnfoldResponse)
  Int_t _Nits;
  mutable Bool_t _unfolded, _haveCov,_have_err_mat, _fail;
  const RooUnfoldResponse* _res;   // Response matrix (not owned)
  const TH1* _meas;                // Measured distribution (not owned)
  mutable TVectorD _rec;  // Reconstructed distribution
  mutable TMatrixD _cov;  // Reconstructed distribution covariance
  mutable TMatrixD _err_mat;
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
inline Int_t                    RooUnfold::Nits()   const { return _Nits; }
inline const RooUnfoldResponse* RooUnfold::response()  const { return _res;     }
inline const TH1*               RooUnfold::Hmeasured() const { return _meas;    }
inline TVectorD&                RooUnfold::Vreco()           { if (!_unfolded) Unfold(); return _rec; }
inline TMatrixD&                RooUnfold::Ereco()           { if (!_haveCov)  GetCov(); return _cov; }
inline TMatrixD&                RooUnfold::Freco()           { if (!_have_err_mat)  Get_err_mat(); return _err_mat; }
inline TObject*                 RooUnfold::Impl()            { return 0; };
inline void  RooUnfold::SetVerbose (Int_t level)             { _verbose= level; }
inline void  RooUnfold::SetNits (Int_t iterations)             { _Nits= iterations; }
inline void  RooUnfold::SetRegParm (Int_t)                   {}
inline Int_t RooUnfold::GetRegParm() const                   {return -1;}

#endif
