//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.h,v 1.35 2010-09-13 21:19:09 adye Exp $
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
#include "RooUnfoldResponse.h"

class TH1;
class TH1D;

class RooUnfold : public TNamed {

public:

  enum Algorithm { kNone, kBayes, kSVD, kBinByBin, kTUnfold,kInvert}; // Selection of unfolding algorithm.
  enum ErrorTreatment {kNoError, kErrors, kCovariance, kCovToy};
  static RooUnfold* New (Algorithm alg, const RooUnfoldResponse* res, const TH1* meas, Double_t regparm= -1e30,
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
  virtual void SetMeasured (const TH1* meas);
  virtual void SetResponse (const RooUnfoldResponse* res);
  virtual void Reset ();

  // Accessors

  virtual const RooUnfoldResponse* response() const;
  virtual const TH1* Hmeasured() const;
  virtual       TH1* Hreco (ErrorTreatment withError=kErrors);
  const    TVectorD& Vmeasured() const;   // Measured distribution as a TVectorD
  const    TVectorD& Emeasured() const;   // Measured distribution errors as a TVectorD

  virtual TVectorD&  Vreco();
  virtual TMatrixD   Ereco  (ErrorTreatment witherror=kCovariance);
  virtual TVectorD   ErecoV (ErrorTreatment witherror=kErrors);

  virtual Int_t      verbose() const;
  virtual void       SetVerbose (Int_t level);
  virtual Int_t      NToys() const;         // Number of toys
  virtual void       SetNToys (Int_t toys); // Set number of toys
  virtual Int_t      Overflow() const;
  virtual void       PrintTable (std::ostream& o, const TH1* hTrue= 0, ErrorTreatment=kNoError);
  virtual void       SetRegParm (Double_t parm);
  virtual Double_t   GetRegParm() const; // Get Regularisation Parameter
  Double_t Chi2 (const TH1* hTrue,ErrorTreatment DoChi2=kCovariance);
  Double_t GetMinParm() const;
  Double_t GetMaxParm() const;
  Double_t GetStepSizeParm() const;
  Double_t GetDefaultParm() const;
  TH1* Runtoy(ErrorTreatment doerror=kNoError,double* chi2=0,const TH1* hTrue=0) const;
  void Print(Option_t *opt="")const;

protected:
  void Assign (const RooUnfold& rhs); // implementation of assignment operator
  virtual void SetNameTitleDefault(); 
  virtual void Unfold(); 
  virtual void GetErrors();
  virtual void GetCov(); // Get covariance matrix using errors on measured distribution
  virtual void GetErrMat(); // Get covariance matrix using errors from residuals on reconstructed distribution
  virtual void GetSettings();
  virtual Bool_t UnfoldWithErrors (ErrorTreatment withError);
  TH1D* HmeasuredNoOverflow1D();

  static TH1*     Add_Random   (const TH1* meas);
  static TMatrixD CutZeros     (const TMatrixD& ereco);
  static TH1D*    HistNoOverflow (const TH1* h, Bool_t overflow);

private:
  void Init();
  void Destroy();
  void CopyData (const RooUnfold& rhs);

protected:
  // instance variables
  Double_t _minparm; //Minimum value to be used in RooUnfoldParms
  Double_t _maxparm; //Maximum value to be used in RooUnfoldParms
  Double_t _stepsizeparm; //StepSize value to be used in RooUnfoldParms
  Double_t _defaultparm; //Recommended value for regularisation parameter
  Int_t _verbose;  // Debug print level
  Int_t _nm;   // Total number of measured bins (including under/overflows if _overflow set)
  Int_t _nt;   // Total number of truth    bins (including under/overflows if _overflow set)
  Int_t _overflow;   // Use histogram under/overflows if 1 (set from RooUnfoldResponse)
  Int_t _NToys; // Number of toys to be used
  Bool_t _unfolded, _haveCov,_have_err_mat, _fail, _haveErrors;
  const RooUnfoldResponse* _res;   // Response matrix (not owned)
  const TH1* _meas;                // Measured distribution (not owned)
  TVectorD _rec;       // Reconstructed distribution
  TMatrixD _cov;       // Reconstructed distribution covariance
  TVectorD _variances; // Error matrix diagonals
  TMatrixD _err_mat;   // Error matrix (from toys)
  mutable TVectorD *_vMes, *_eMes; //! Cached measured vector and error
  
public:

  ClassDef (RooUnfold, 0) // Unfold
};

// Inline method definitions

inline RooUnfold::RooUnfold()                                           : TNamed()           {Init();}
inline RooUnfold::RooUnfold (const char*    name, const char*    title) : TNamed(name,title) {Init();}
inline RooUnfold::RooUnfold (const TString& name, const TString& title) : TNamed(name,title) {Init();}
inline RooUnfold::~RooUnfold() {Destroy();}
inline RooUnfold& RooUnfold::operator= (const RooUnfold& rhs) {Assign(rhs); return *this;}

inline Int_t                    RooUnfold::verbose()   const { return _verbose; } // Controls amount of information to be printed
inline Int_t                    RooUnfold::NToys()     const { return _NToys;   } // Sets Number of toys
inline Int_t                    RooUnfold::Overflow()  const { return _overflow;}
inline const RooUnfoldResponse* RooUnfold::response()  const { return _res;     } // Response object
inline const TH1*               RooUnfold::Hmeasured() const { return _meas;    } // Measured Distribution
inline TVectorD&                RooUnfold::Vreco()           { if (!_unfolded) Unfold(); return _rec; } // Vector of reconstructed points
inline const TVectorD&          RooUnfold::Vmeasured() const { if (!_vMes) _vMes= RooUnfoldResponse::H2V  (_meas, _res->GetNbinsMeasured(), _overflow); return *_vMes; }
inline const TVectorD&          RooUnfold::Emeasured() const { if (!_eMes) _eMes= RooUnfoldResponse::H2VE (_meas, _res->GetNbinsMeasured(), _overflow); return *_eMes; }
inline void  RooUnfold::SetVerbose (Int_t level)             { _verbose= level; } // Set verbose
inline void  RooUnfold::SetNToys (Int_t toys)                { _NToys= toys; } // Set toys
inline void  RooUnfold::SetRegParm (Double_t)                {} // Set Regularisation parameter
inline Double_t RooUnfold::GetRegParm() const                {return -1;}

#endif
