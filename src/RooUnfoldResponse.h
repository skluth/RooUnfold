//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldResponse.h,v 1.1.1.1 2007-04-04 21:27:02 adye Exp $
//
// Description:
//      Response Matrix
//
// Author List:
//      Tim Adye <T.J.Adye@rl.ac.uk>
//
// Copyright Information:
//      Copyleft (C) 2006 Rutherford Appleton Laboratory
//
//==============================================================================

#ifndef ROOUNFOLDRESPONSE_HH
#define ROOUNFOLDRESPONSE_HH

#include "TNamed.h"
#include "TMatrixD.h"

class TH1;
class TH1D;
class TH2D;
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
#include "TVectorDfwd.h"
#else
class TVectorD;
#endif

class RooUnfoldResponse : public TNamed {

public:

  // Standard methods

  RooUnfoldResponse(); // default constructor
  RooUnfoldResponse (const char* name, const char* title); // named constructor
  RooUnfoldResponse (const TString& name, const TString& title); // named constructor
  RooUnfoldResponse (const RooUnfoldResponse& rhs); // copy constructor
  virtual ~RooUnfoldResponse(); // destructor
  virtual RooUnfoldResponse& operator= (const RooUnfoldResponse& rhs); // assignment operator

  // Special constructors

  RooUnfoldResponse (Int_t nb, Double_t xlo, Double_t xhi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case with same binning, measured vs truth
  RooUnfoldResponse (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case
  RooUnfoldResponse (const TH1* measured, const TH1* truth, const char* name= 0, const char* title= 0);  // constructor - measured and truth only used for shape
  RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2D* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms

  // Set up an existing object

  virtual RooUnfoldResponse& Clear ();  // clear an existing object
  virtual RooUnfoldResponse& Setup (const RooUnfoldResponse& rhs);  // set up based on another instance
  virtual RooUnfoldResponse& Setup (Int_t nb, Double_t xlo, Double_t xhi);  // set up simple 1D case with same binning, measured vs truth
  virtual RooUnfoldResponse& Setup (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi);  // set up simple 1D case
  virtual RooUnfoldResponse& Setup (const TH1* measured, const TH1* truth);  // set up - measured and truth only used for shape
  virtual RooUnfoldResponse& Setup (const TH1* measured, const TH1* truth, const TH2D* response);  // set up from already-filled histograms

  // Fill with training data

  virtual Int_t Fill (Double_t xr, Double_t xt);  // Fill 1D Response Matrix
  virtual Int_t Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt);  // Fill 2D Response Matrix
  virtual Int_t Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt);  // Fill 3D Response Matrix

  virtual Int_t Miss (Double_t xt);  // Fill missed event into 1D Response Matrix
  virtual Int_t Miss (Double_t xt, Double_t yt);  // Fill missed event into 2D Response Matrix
  virtual Int_t Miss (Double_t xt, Double_t yt, Double_t zt);  // Fill missed event into 3D Response Matrix

  // Accessors

  Int_t        GetDimensionMeasured() const;
  Int_t        GetDimensionTruth()    const;
  Int_t        GetNbinsMeasured()     const;
  Int_t        GetNbinsTruth()        const;

  const TH1*   Hmeasured()            const;
  TH1*         Hmeasured();
  const TH1*   Htruth()               const;
  TH1*         Htruth();
  const TH2D*  Hresponse()            const;
  TH2D*        Hresponse();

  TH1D*        Hmeasured1D()          const;
  TH1D*        Htruth1D()             const;

  const TVectorD& Vmeasured()         const;
  const TVectorD& Emeasured()         const;
  const TVectorD& Vtruth()            const;
  const TVectorD& Etruth()            const;
  const TMatrixD& Mresponse()         const;
  const TMatrixD& Eresponse()         const;

  Double_t operator() (Int_t r, Int_t t) const;

  static TH1D*     H2H1D(const TH1*  h, Int_t nb);
  static TVectorD* H2V  (const TH1*  h, Int_t nb);
  static TVectorD* H2VE (const TH1*  h, Int_t nb);
  static TMatrixD* H2M  (const TH2D* h, Int_t nx, Int_t ny, const TH1* norm= 0);
  static TMatrixD* H2ME (const TH2D* h, Int_t nx, Int_t ny, const TH1* norm= 0);

private:

  virtual RooUnfoldResponse& Setup();
  virtual void ClearCache();
  virtual void SetNameTitleDefault();

  // instance variables

  Int_t _mdim; // Number of measured  dimensions
  Int_t _tdim; // Number of truth dimensions
  Int_t _nm;   // Total number of measured  bins
  Int_t _nt;   // Total number of truth bins
  TH1*  _mes;  // Measured histogram
  TH1*  _tru;  // Truth histogram
  TH2D* _res;  // Response histogram

  mutable TVectorD *_vMes, *_eMes; //! Cached measured     vector and error
  mutable TVectorD *_vTru, *_eTru; //! Cached truth    vector and error
  mutable TMatrixD *_mRes, *_eRes; //! Cached response matrix and error
  mutable Bool_t _cached;          //! We are using cached vectors/matrices

public:

  ClassDef (RooUnfoldResponse, 1) // Respose Matrix
};

// Inline method definitions

inline RooUnfoldResponse::RooUnfoldResponse()                                           : TNamed()           {Setup();}
inline RooUnfoldResponse::RooUnfoldResponse (const char*    name, const char*    title) : TNamed(name,title) {Setup();}
inline RooUnfoldResponse::RooUnfoldResponse (const TString& name, const TString& title) : TNamed(name,title) {Setup();}
inline RooUnfoldResponse::~RooUnfoldResponse()                                                               {Clear();}

inline RooUnfoldResponse& RooUnfoldResponse::Setup (Int_t nb, Double_t xlo, Double_t xhi) { return Setup (nb, xlo, xhi, nb, xlo, xhi); }

inline Int_t        RooUnfoldResponse::GetDimensionMeasured() const { return _mdim; }
inline Int_t        RooUnfoldResponse::GetDimensionTruth()    const { return _tdim; }
inline Int_t        RooUnfoldResponse::GetNbinsMeasured()     const { return _nm;   }
inline Int_t        RooUnfoldResponse::GetNbinsTruth()        const { return _nt;   }

inline const TH1*   RooUnfoldResponse::Hmeasured()            const { return _mes;  }
inline TH1*         RooUnfoldResponse::Hmeasured()                  { return _mes;  }
inline const TH1*   RooUnfoldResponse::Htruth()               const { return _tru;  }
inline TH1*         RooUnfoldResponse::Htruth()                     { return _tru;  }
inline const TH2D*  RooUnfoldResponse::Hresponse()            const { return _res;  }
inline TH2D*        RooUnfoldResponse::Hresponse()                  { return _res;  }

inline TH1D*        RooUnfoldResponse::Hmeasured1D()          const { return H2H1D (_mes, _nm); }
inline TH1D*        RooUnfoldResponse::Htruth1D()             const { return H2H1D (_tru, _nt); }

inline const TVectorD& RooUnfoldResponse::Vmeasured()         const { if (!_vMes) _cached= _vMes= H2V  (_mes, _nm); return *_vMes; }
inline const TVectorD& RooUnfoldResponse::Emeasured()         const { if (!_eMes) _cached= _eMes= H2VE (_mes, _nm); return *_eMes; }
inline const TVectorD& RooUnfoldResponse::Vtruth()            const { if (!_vTru) _cached= _vTru= H2V  (_tru, _nt); return *_vTru; }
inline const TVectorD& RooUnfoldResponse::Etruth()            const { if (!_eTru) _cached= _eTru= H2VE (_tru, _nt); return *_eTru; }
inline const TMatrixD& RooUnfoldResponse::Mresponse()         const { if (!_mRes) _cached= _mRes= H2M  (_res, _nm, _nt, _tru); return *_mRes; }
inline const TMatrixD& RooUnfoldResponse::Eresponse()         const { if (!_eRes) _cached= _eRes= H2ME (_res, _nm, _nt, _tru); return *_eRes; }

inline Double_t RooUnfoldResponse::operator() (Int_t r, Int_t t) const { return Mresponse()(r,t); }

#endif
