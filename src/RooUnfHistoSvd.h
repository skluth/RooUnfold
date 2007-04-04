///////////////////////////////////////////////////////////////////////
//Kerstin Tackmann, Heiko Lacker (TU Dresden)
//based on 
//Andreas Hoecker, Vakhtang Kartvelishvili, hep-ph/9509307
//$Id: RooUnfHistoSvd.h,v 1.1.1.1 2007-04-04 21:27:01 adye Exp $
///////////////////////////////////////////////////////////////////////

#ifndef TUNFHISTO_HH
#define TUNFHISTO_HH

#include "TObject.h"
#include "TMatrixD.h"
#include "TVectorD.h"

class TH1D;
class TH2D;
class TFile;

class TUnfHisto : public TObject 
{

  public :
    TUnfHisto(Int_t gDim, Int_t rDim);
  ~TUnfHisto(); 

  void init(const TH1D *b, const TH1D *bini, const TH1D *xini, const TH2D *A, Bool_t nofile= false); //Unfolding of data
  void init(const TH1D *b, const TH1D *bini, const TH1D *xini, const TH2D *A, const TH1D *xref, Bool_t nofile= false); //Unfolding of MC 
  TVectorD Unfold(Int_t tau=1, Int_t toy=0, Int_t mattoy=0);
  TMatrixD GetCov(const TMatrixD& cov, const TH1D *bref, Int_t ntoys, Int_t tau, TH1D *toydist[16]=NULL, Int_t fg=0);
  void GetS2(const TMatrixD& cov, const TMatrixD& mcov, const TH1D *xref, Int_t ntoys, Int_t tau, Int_t fg, TH1D *s2dist);
  TMatrixD GetMatStatCov(Int_t ntoys, Int_t tau, Int_t fg=0);


  private : 

  static Double_t VecMax(const TVectorD& vec);
  static Double_t VecMin(const TVectorD& vec);
  void H2M(const TH2D *histo, TMatrixD& mat) const;
  static void H2V(const TH1D* histo, TVectorD& vec);
  static void H2Verr(const TH1D* histo, TVectorD& vec);
  static void V2H(const TVectorD& vec, TString hname);
  static Double_t SNorm(const TVectorD& vec, Int_t dim);
  static TVectorD VecDiv(const TVectorD& vec1, const TVectorD& vec2, Int_t zero=0);
  TMatrixD MatDivVec(const TMatrixD& mat, const TVectorD& vec, Int_t zero=0) const;
  TVectorD MatMultVec(const TMatrixD& mat, const TVectorD& vec, Int_t resdim) const;
  static TVectorD CompProd(const TVectorD& vec1, const TVectorD& vec2);
  void fillC(TMatrixD& tCurv, TMatrixD& tC, Int_t matdim, Int_t ndim) const;
  Double_t GetCurvature(const TVectorD& vec, const TMatrixD& curv) const;
  Double_t GetChi2(const TVectorD& vec, const TVectorD& refvec, const TVectorD& err, Int_t n) const;
  void WriteResults(Int_t htau);
  void InitHistos(Int_t htau);

  //Members
  const Int_t _gDim, _rDim;
  Int_t _ddim, _MC; 
  TFile *_resfile;
  //The input histos
  const TH1D *_b;
  const TH1D *_bini;
  const TH1D *_xini;
  const TH1D *_xref;
  const TH2D *_A;
  TH1D *_toyhisto;
  TH2D *_toymat;

  ClassDef(TUnfHisto,0) // Unfolding Histograms Using Singular Value Decomposition
};

#endif
