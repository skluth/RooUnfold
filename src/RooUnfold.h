//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.h,v 1.17 2010-08-06 15:37:25 fwx38934 Exp $
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

  enum Algorithm { kNone, kBayes, kSVD, kBinByBin, kTUnfold}; // Selection of unfolding algorithm.

  static RooUnfold* New (Algorithm alg, const RooUnfoldResponse* res, const TH1* meas, Double_t regparm= 1,
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
  virtual Int_t                    Nits() const; // Number of iterations
  virtual void SetNits (Int_t iterations); // Set number of iterations

  virtual void PrintTable (std::ostream& o, const TH1* hTrue= 0, Int_t withError= 1);

  virtual TObject* Impl();

  virtual void  SetRegParm (Double_t parm);
  virtual Double_t GetRegParm() const; // Get Regularisation Parameter
 virtual void Get_settings();
  Double_t Chi2 (const TH1* hTrue,Int_t DoChi2=1);
  Double_t Get_minparm();
  Double_t Get_maxparm();
  Double_t Get_stepsizeparm();
  Double_t Get_defaultparm();
  TH1* Runtoy(Int_t doerror=0,double* chi2=0,const TH1* hTrue=0);
  
  Double_t _minparm; //Minimum value to be used in RooUnfoldParms
  Double_t _maxparm; //Maximum value to be used in RooUnfoldParms
  Double_t _stepsizeparm; //StepSize value to be used in RooUnfoldParms
  Double_t _defaultparm; //Recommended value for regularisation parameter
  
  protected:

  void Init();
  virtual void Unfold(); 
  virtual void GetCov(); // Get covariance matrix using errors on measured distribution
  virtual void SetNameTitleDefault(); 
  virtual void Get_err_mat(); // Get covariance matrix using errors from residuals on reconstructed distribution
  void Assign   (const RooUnfold& rhs); // implementation of assignment operator
  void CopyData (const RooUnfold& rhs);
  TH1* Add_Random(const TH1* hMeas_AR);
  // instance variables

  Int_t _verbose;  // Debug print level
  Int_t _nm;   // Total number of measured bins
  Int_t _nt;   // Total number of truth    bins
  Int_t _overflow;   // Use histogram under/overflows if 1 (set from RooUnfoldResponse)
  Int_t _Nits; // Number of iterations to be used
  mutable Bool_t _unfolded, _haveCov,_have_err_mat, _fail;
  const RooUnfoldResponse* _res;   // Response matrix (not owned)
  const TH1* _meas;                // Measured distribution (not owned)
  mutable TVectorD _rec;  // Reconstructed distribution
  mutable TMatrixD _cov;  // Reconstructed distribution covariance
  mutable TMatrixD _err_mat; // Error Matrix
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

inline Int_t                    RooUnfold::verbose()   const { return _verbose; } // Controls amount of information to be printed
inline Int_t                    RooUnfold::Nits()   const { return _Nits; } // Sets Number of iterations
inline const RooUnfoldResponse* RooUnfold::response()  const { return _res;     } // Response object
inline const TH1*               RooUnfold::Hmeasured() const { return _meas;    } // Measured Distribution
inline TVectorD&                RooUnfold::Vreco()           { if (!_unfolded) Unfold(); return _rec; } // Vector or reconstructed points
inline TMatrixD&                RooUnfold::Ereco()           { if (!_haveCov)  GetCov(); return _cov; } // Covariance matrix from measured distribution
inline TMatrixD&                RooUnfold::Freco()           { if (!_have_err_mat)  Get_err_mat(); return _err_mat; } // Covariance matrix from residuals in reconstructed distribution 
inline TObject*                 RooUnfold::Impl()            { return 0; };
inline void  RooUnfold::SetVerbose (Int_t level)             { _verbose= level; } // Set verbose
inline void  RooUnfold::SetNits (Int_t iterations)             { _Nits= iterations; } // Set iterations
inline void  RooUnfold::SetRegParm (Double_t)                   {} // Set Regularisation parameter
inline Double_t RooUnfold::GetRegParm() const                   {return -1;}

#endif
