///////////////////////////////////////////////////////////////////////
//Kerstin Tackmann, Heiko Lacker (TU Dresden)
//based on 
//Andreas Hoecker, Vakhtang Kartvelishvili, hep-ph/9509307
//$Id: RooUnfHistoSvd.cxx,v 1.1.1.1 2007-04-04 21:27:25 adye Exp $
///////////////////////////////////////////////////////////////////////

#include "RooUnfHistoSvd.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TRandom.h"
#include <iostream>

ClassImp(TUnfHisto);

using namespace std;

TUnfHisto::TUnfHisto(Int_t gDim, Int_t rDim): 
  _gDim(gDim),
  _rDim(gDim),
  _ddim(2),
  _MC(0),
  _resfile(NULL),
  _b(NULL), 
  _bini(NULL),
  _xini(NULL),
  _xref(NULL),  
  _A(NULL), 
  _toyhisto(NULL),
  _toymat(NULL) {}

TUnfHisto::~TUnfHisto()
{
  delete _resfile; _resfile = 0;
}

//For the unfolding of data
void TUnfHisto::init(const TH1D *b, const TH1D *bini, const TH1D *xini, const TH2D *A, Bool_t nofile)
{
  if(b->GetNbinsX() != xini->GetNbinsX()){
    cout << "Only same numbers of bins for reconstructed and true spectrum are supported, results will be useless!" << endl;
  }
  // Get the input histos
  _b = b;
  _bini = bini;
  _xini = xini;
  _A = A;
  _ddim = 2;
  if (!nofile)
    _resfile = new TFile("UnfHisto.root","RECREATE");
}

//For the unfolding of MC, i.e. expected result is xref
void TUnfHisto::init(const TH1D *b, const TH1D *bini, const TH1D *xini, const TH2D *A, const TH1D *xref, Bool_t nofile)
{
  //  cout << "Initialization for MC unfolding" << endl;
  if(b->GetNbinsX() != xini->GetNbinsX()){
    cout << "Only same numbers of bins for reconstructed and true spectrum are supported, results will be useless!" << endl;
  }
  // Get the input histos
  _b = b;
  _bini = bini;
  _xini = xini;
  _A = A;
  _xref = xref;
  _MC = 1;
  _ddim = 2;
  if (!nofile)
    _resfile = new TFile("UnfHisto.root","RECREATE");
}

TVectorD TUnfHisto::Unfold(Int_t tau, Int_t toy, Int_t mattoy)
{
  //  cout << "Unfolding with regularization parameter " << tau << endl;

  //Make the histos
  if(!toy && !mattoy){
    InitHistos(tau);
  }

  //The vectors and matrices
  TVectorD vb(_rDim), vbini(_rDim), vxini(_gDim), vberr(_rDim), vxref(_gDim),
           vdini(_rDim), vd(_rDim), vdz(_gDim), vzini(_gDim), vz(_gDim), 
           vwini(_gDim), vw(_gDim), vx(_gDim); 
  TMatrixD mA(_rDim, _gDim), mCurv(_gDim, _gDim), mC(_gDim, _gDim), 
           Uort(_gDim, _gDim), Vort(_rDim, _rDim), Vreg(_rDim, _rDim);

  Double_t eps=1e-12;
  Double_t sreg;
  char hname[200];

  if(toy){
    H2V(_toyhisto, vb);
    H2Verr(_toyhisto, vberr);
  }
  else{
    H2V(_b, vb);
    H2Verr(_b, vberr);
  }
  H2V(_bini, vbini);
  H2V(_xini, vxini);
  if(mattoy){
    H2M(_toymat, mA);
  } else {
    H2M(_A, mA);
  }
  if(_MC){H2V(_xref, vxref);}

  //Scale the MC vectors to data norm
  Double_t Scale = SNorm(vb, _rDim)/SNorm(vbini, _rDim);
  for(Int_t i=0; i<_rDim; i++){
    vbini(i) = vbini(i) * Scale;
  }
  for(Int_t i=0; i<_gDim; i++){
    vxini(i) = vxini(i) * Scale;
    for(Int_t j=0; j<_rDim; j++){
      mA(j,i) = mA(j,i) * Scale;
    }
  }

  //Normalize xref to one if necessary
  if(_MC){
    Scale = SNorm(vxref, _gDim);
    if(Scale!=0){
      for(Int_t i=0; i<_gDim; i++){
        vxref(i) = vxref(i) / Scale;
      }
    }
    else{
      cout << "Norm of reference distribution in zero!" << endl;
    }
  }

  //Fill and invert the second derivative matrix
  fillC(mCurv, mC, _gDim, _ddim);
  TMatrixD mCinv(_gDim, _gDim);

  //Inversion of mC by help of SVD
  TDecompSVD CSVD(mC);
  TMatrixD CUort = CSVD.GetU();
  TMatrixD CVort = CSVD.GetV();
  TVectorD CSV = CSVD.GetSig();

  TMatrixD CSVM(_gDim, _gDim);
  for(Int_t i=0; i<_gDim; i++){
    CSVM(i,i) = 1/CSV(i);
  }

  CUort.Transpose(CUort);
  TMatrixD test(_gDim, _gDim);
  test.Mult(CVort,CSVM);
  mCinv.Mult(test,CUort);

  //Rescale matrix and vectors by error of data vector
  vbini = VecDiv(vbini, vberr);
  vb = VecDiv(vb, vberr, 1);
  mA = MatDivVec(mA, vberr, 1);
  vberr = VecDiv(vberr, vberr, 1);

  //Singular value decomposition and matrix operations
  Uort.Mult(mA,mCinv);

  TDecompSVD ASVD(Uort);
  Uort = ASVD.GetU();
  Vort = ASVD.GetV();
  TVectorD ASV = ASVD.GetSig();

  if(!toy && !mattoy){
    sprintf(hname, "%s", "sv");
    V2H(ASV, hname);
  }

  Vreg.Mult(mCinv, Vort);
  Uort.Transpose(Uort);
  vdini = MatMultVec(Uort, vbini, _rDim);
  vd = MatMultVec(Uort, vb, _rDim);

  if(!toy && !mattoy){
    sprintf(hname, "%s", "dd");
    V2H(vd, hname);
  }

  //The damping coefficient(s)
  Int_t mintau, maxtau;
  mintau = tau-1; 
  maxtau = tau;
  if(tau==0){ //If tau==0 then loop over all singular values
    mintau = 0;
    maxtau = _gDim;
  }

  for(Int_t k=mintau; k<maxtau; k++){
    //The damping factors
    for(Int_t i=0; i<_gDim; i++){
      if(ASV(i)<ASV(0)*eps){
        sreg = ASV(0)*eps;
      }
      else{
        sreg = ASV(i);
      }
	vdz(i) = sreg/(sreg*sreg + ASV(k)*ASV(k));
    }
    vzini = CompProd(vdini, vdz);
    vz = CompProd(vd, vdz);

    if(!toy && !mattoy){
      sprintf(hname, "%s%d", "z", k+1);
      V2H(vz, hname);
    }

    //Compute the weights
    vwini = MatMultVec(Vreg, vzini, _gDim);
    vw = MatMultVec(Vreg, vz, _gDim);

    if(!toy && !mattoy){
      sprintf(hname, "%s%d", "w", k+1);
      V2H(vw, hname);
    }

    //Rescale by xini
    vx = CompProd(vw, vxini);

    //Scale result to unit area
    Scale = SNorm(vx, _gDim);
    if(Scale!=0){
      for(Int_t i=0; i<_gDim; i++){
        vx(i) = vx(i)/Scale;
      }
    }
    if(!toy && !mattoy){
      sprintf(hname, "%s%d", "x", k+1);
      V2H(vx, hname);
    }

    //Get Curvature and also chi2 in case of MC unfolding
    if(!toy && !mattoy){
      if(_MC){
        cout << "Unfolding param " << k+1 << endl;
        cout << "******************" << endl;
        cout << "Curvature of weight distribution " << GetCurvature(vw, mCurv) << endl;
        cout << "Reduced chi2 result " << GetChi2(vx, vxref, vberr, _gDim) << endl;
        cout << endl;
      }
      else{
        cout << "Unfolding param " << k+1 << endl;
        cout << "******************" << endl;
        cout << "Curvature of weight distribution " << GetCurvature(vw, mCurv) << endl;
        cout << endl;
      }
    }
  }
  if(!toy && !mattoy){
    WriteResults(tau);
  }
  return vx;
}

TMatrixD TUnfHisto::GetCov(const TMatrixD& cov, const TH1D *bref,  Int_t ntoys, Int_t tau, TH1D* toydist[16], Int_t fg)
{
//Code for generation of toys taken from RooResult and modified

  TVectorD unfres(_gDim);
  TMatrixD unfcov(_gDim, _gDim);

  // calculate the elements of the upper-triangular matrix L that
  // gives Lt*L = C
  // where Lt is the transpose of L (the "square-root method")  
  Int_t nPar = bref->GetNbinsX();
  TMatrixD L(nPar,nPar);
  for(int i=0; i<nPar; i++)
    for(int j=0; j<nPar; j++)
      L(i,j) = 0.;

  for(Int_t iPar= 0; iPar < nPar; iPar++) {
    // calculate the diagonal term first
    L(iPar,iPar)= cov(iPar,iPar);
    for(Int_t k= 0; k < iPar; k++) {
      Double_t tmp= L(k,iPar);
      L(iPar,iPar)-= tmp*tmp;
    }
    if(L(iPar,iPar)>=0.0){
      L(iPar,iPar)= sqrt(L(iPar,iPar));
    }
    else{
      L(iPar,iPar)=0.0;
    }
    // then the off-diagonal terms
    for(Int_t jPar= iPar+1; jPar < nPar; jPar++) {
      L(iPar,jPar)= cov(iPar,jPar);
      for(Int_t k= 0; k < iPar; k++) {
        L(iPar,jPar)-= L(k,iPar)*L(k,jPar);
      }
      if(L(iPar,iPar)!=0.)
        L(iPar,jPar)/= L(iPar,iPar);
      else
        L(iPar,jPar) = 0.;
    }
  }

  // remember Lt
  TMatrixD *Lt= new TMatrixD(TMatrixD::kTransposed,L);
  TRandom random(123);
  for(Int_t i=1; i<16*fg; i++){
    random.Gaus(0.,1.);
  }

  _toyhisto = new TH1D(*bref);
  TH1D *toymean = new TH1D(*bref);
  for(Int_t j=0; j<_gDim; j++){
    toymean->SetBinContent(j+1,0.);
  }

  //Get the mean of the toys first
  for(int i=1; i<=ntoys; i++){

    // create a vector of unit Gaussian variables
    TVectorD g(nPar);
    for(Int_t k= 0; k < nPar; k++) {
      g(k) = random.Gaus(0.,1.);
    }
    // multiply this vector by Lt to introduce the appropriate
    // correlations
    g*= (*Lt);
    // add the mean value offsets and store the results
    for(int j=1; j<=nPar; j++){
      _toyhisto->SetBinContent(j,bref->GetBinContent(j)+g(j-1));
      _toyhisto->SetBinError(j,bref->GetBinError(j));
    }
    unfres = Unfold(tau, 1, 0);
    if(VecMax(unfres)>1.){
      cout << "Max of toy " << i << "  " << VecMax(unfres) << endl;
      if(toydist != NULL){
        for(int j=0; j<_gDim; j++)
          toydist[j]->Fill(unfres(j));
	i--;
	continue;
      }
    }
    for(Int_t j=0; j<_gDim; j++){
      toymean->SetBinContent(j+1, toymean->GetBinContent(j+1) + unfres(j)/ntoys);
      if(toydist != NULL){
        toydist[j]->Fill(unfres(j));
      }
    }
  }

  //Reset the random seed
  random.SetSeed(123);
  for(Int_t i=1; i<16*fg; i++){
    random.Gaus(0.,1.);
  }


  //Now the toys for the covariance matrix
  for(int i=1; i<=ntoys; i++){

    // create a vector of unit Gaussian variables
    TVectorD g(nPar);
    for(Int_t k= 0; k < nPar; k++) g(k) = random.Gaus(0.,1.);
    // multiply this vector by Lt to introduce the appropriate
    // correlations
    g*= (*Lt);
    // add the mean value offsets and store the results
    for(int j=1; j<=nPar; j++){
      _toyhisto->SetBinContent(j,bref->GetBinContent(j)+g(j-1));
      _toyhisto->SetBinError(j,bref->GetBinError(j));
    }
    unfres = Unfold(tau, 1, 0);
    if(VecMax(unfres)>1.){
      cout << "Max of toy " << i << "  " << VecMax(unfres) << endl;
      i--;
      continue;
    }
    for(Int_t j=0; j<_gDim; j++){
      for(Int_t k=0; k<_gDim; k++){
        unfcov(j,k) = unfcov(j,k) + (unfres(j) - toymean->GetBinContent(j+1))*(unfres(k) - toymean->GetBinContent(k+1))/(ntoys-1);
      }
    }
  }
  delete Lt;
  return unfcov;
}


void TUnfHisto::GetS2(const TMatrixD& cov, const TMatrixD& mcov, const TH1D *xref, Int_t ntoys, Int_t tau, Int_t fg, TH1D *s2dist)
{
//Code for generation of toys taken from RooResult and modified

  TVectorD unfres(_gDim); 
  TVectorD truthtoy(_gDim);
  TVectorD bias(_gDim);
  Double_t S2=0.;
  Double_t truthtoynorm=0.;

  // calculate the elements of the upper-triangular matrix L that
  // gives Lt*L = C
  // where Lt is the transpose of L (the "square-root method")  
  Int_t nPar = xref->GetNbinsX();
  TMatrixD L(nPar,nPar);
  for(int i=0; i<nPar; i++)
    for(int j=0; j<nPar; j++)
      L(i,j) = 0.;

  for(Int_t iPar= 0; iPar < nPar; iPar++) {
    // calculate the diagonal term first
    L(iPar,iPar)= cov(iPar,iPar);
    for(Int_t k= 0; k < iPar; k++) {
      Double_t tmp= L(k,iPar);
      L(iPar,iPar)-= tmp*tmp;
    }
    if(L(iPar,iPar)<0){
      L(iPar,iPar) = 0.;
    }
    else{
      L(iPar,iPar)= sqrt(L(iPar,iPar));
    }
    // then the off-diagonal terms
    for(Int_t jPar= iPar+1; jPar < nPar; jPar++) {
      L(iPar,jPar)= cov(iPar,jPar);
      for(Int_t k= 0; k < iPar; k++) {
        L(iPar,jPar)-= L(k,iPar)*L(k,jPar);
      }
      if(L(iPar,iPar)!=0.)
        L(iPar,jPar)/= L(iPar,iPar);
      else
        L(iPar,jPar) = 0.;
    }
  }

  // remember Lt
  TMatrixD *Lt= new TMatrixD(TMatrixD::kTransposed,L);

  TRandom random(123);
  for(Int_t i=1; i<16*fg; i++){
    random.Gaus(0.,1.);
  }

  _toyhisto = new TH1D(*xref);
  TH1D *biashisto = new TH1D(*xref);
  TH1D *truthhisto = new TH1D(*xref);

  //Get the mean of the toys first
  for(int i=1; i<=ntoys; i++){
    // create a vector of unit Gaussian variables
    TVectorD g(nPar);
    for(Int_t k= 0; k < nPar; k++) g(k) = random.Gaus(0.,1.);
    // Was: g(k)= RooRandom::gaussian();
    // multiply this vector by Lt to introduce the appropriate
    // correlations
    g*= (*Lt);
    // add the mean value offsets and store the results
    for(int j=0; j<nPar; j++){
      truthtoy(j) = xref->GetBinContent(j+1)+g(j);
    }
    for(int k=0; k<nPar; k++){
      _toyhisto->SetBinContent(i+1, 0.);
      for(int j=0; j<nPar; j++){
        _toyhisto->SetBinContent(k+1, _toyhisto->GetBinContent(k+1) + _A->GetBinContent(k+1,j+1)*truthtoy(j));
      }
      _toyhisto->SetBinError(k+1, sqrt(mcov(k,k)));
    }

    unfres = Unfold(tau, 1, 0);
    if(VecMax(unfres) > 1.){
      cout << "Max of toy " << i << "  " << VecMax(unfres) << endl;
      for(int j=0; j<nPar; j++)
        cout << j << "  " << _toyhisto->GetBinContent(j+1) << "  " << unfres(j) << endl;
      i--;
      continue;
    }
    //Evaluate the bias...
    truthtoynorm=0.;
    for(int k=0; k<nPar; k++){
      truthtoynorm = truthtoynorm + truthtoy(k);
      biashisto->SetBinContent(k+1, 0.);
    }
    for(int k=0; k<nPar; k++){
      for(int j=0; j<nPar; j++){
        if(truthtoy(j)!=0.){
          biashisto->SetBinContent(k+1, biashisto->GetBinContent(k+1) + _A->GetBinContent(k+1,j+1) * unfres(j) * truthtoynorm / truthtoy(j));
          biashisto->SetBinError(k+1, sqrt(mcov(k,k)));
        }
      }
      truthhisto->SetBinContent(k+1, truthtoy(k)/truthtoynorm);
    }
    if(biashisto->Integral()<0.){
      i--;
      continue;
    }
    S2 = 0.;
    for(Int_t j=0; j<_gDim; j++){
      S2 = S2 + (xref->GetBinContent(j+1)/xref->Integral() - unfres(j))*(xref->GetBinContent(j+1)/xref->Integral() - unfres(j))*xref->Integral()*xref->Integral()/cov(j,j);
    }
    s2dist->Fill(S2);
  }
  delete Lt;
}

TMatrixD TUnfHisto::GetMatStatCov(Int_t ntoys, Int_t tau, Int_t fg)
{
  TVectorD unfres(_gDim);
  TMatrixD unfcov(_gDim, _gDim);

  //Now the toys for the detector response matrix
  TRandom random(123);
  for(Int_t i=1; i<256*fg; i++){
    random.Gaus(0.,1.);
  }

  _toymat = new TH2D(*_A);
  TH1D *toymean = new TH1D(*_xini);
  for(Int_t j=0; j<_gDim; j++){
    toymean->SetBinContent(j+1,0.);
  }

  Int_t nPar = _b->GetNbinsX();

  for(int i=1; i<=ntoys; i++){
    
    for(Int_t k=0; k<nPar; k++){
      for(Int_t m=0; m<nPar; m++){
	if(_A->GetBinContent(k+1,m+1)){
	  _toymat->SetBinContent(k+1, m+1, random.Poisson(_A->GetBinContent(k+1,m+1)));
	}
      }
    }

    unfres = Unfold(tau, 0, 1);
    if(VecMax(unfres)>1.){
      cout << "Max of toy " << i << "  " << VecMax(unfres) << endl;
      i--;
      continue;
    }
    for(Int_t j=0; j<_gDim; j++){
      toymean->SetBinContent(j+1, toymean->GetBinContent(j+1) + unfres(j)/ntoys);
    }
  }

  //Reset the random seed
  random.SetSeed(123);
  for(Int_t i=1; i<256*fg; i++){
    random.Gaus(0.,1.);
  }
  for(int i=1; i<=ntoys; i++){
    for(Int_t k=0; k<nPar; k++){
      for(Int_t m=0; m<nPar; m++){
	if(_A->GetBinContent(k+1,m+1))
	  _toymat->SetBinContent(k+1, m+1, random.Poisson(_A->GetBinContent(k+1,m+1)));
      }
    }

    unfres = Unfold(tau, 0, 1);
    if(VecMax(unfres)>1.){
      cout << "Max of toy " << i << "  " << VecMax(unfres) << endl;
      i--;
      continue;
    }
    for(Int_t j=0; j<_gDim; j++){
      for(Int_t k=0; k<_gDim; k++){
        unfcov(j,k) = unfcov(j,k) + (unfres(j) - toymean->GetBinContent(j+1))*(unfres(k) - toymean->GetBinContent(k+1))/(ntoys-1);
      }
    }
  }
  return unfcov;
}

Double_t TUnfHisto::VecMax(const TVectorD& vec)
{
  Double_t max=0;
  for(Int_t i=0; i<vec.GetNrows(); i++){
    if(max<vec(i)){
      max = vec(i);
    }
  }
  return max;
}

Double_t TUnfHisto::VecMin(const TVectorD& vec)
{
  Double_t min=0;
  for(Int_t i=0; i<vec.GetNrows(); i++){
    if(min>vec(i)){
      min = vec(i);
    }
  }
  return min;
}

void TUnfHisto::H2M(const TH2D* histo, TMatrixD& mat) const
{
  for(Int_t j=0; j<_gDim; j++){
    for(Int_t i=0; i<_rDim; i++){
      mat(i,j) = histo->GetBinContent(i+1,j+1);
    }
  }
}

void TUnfHisto::H2V(const TH1D* histo, TVectorD& vec)
{
  for(Int_t i=0; i<histo->GetNbinsX(); i++){
    vec(i) = histo->GetBinContent(i+1);
  }
}

void TUnfHisto::H2Verr(const TH1D* histo, TVectorD& vec)
{
  for(Int_t i=0; i<histo->GetNbinsX(); i++){
    vec(i) = histo->GetBinError(i+1);
  }
}

Double_t TUnfHisto::SNorm(const TVectorD& vec, Int_t dim)
{
  Double_t snorm = 0;
  for(Int_t i=0; i<dim; i++){
    snorm = snorm + vec(i);
  }
  return snorm;
}

TVectorD TUnfHisto::VecDiv(const TVectorD& vec1, const TVectorD& vec2, Int_t zero)
{
  TVectorD quot(vec1.GetNrows());
  for(Int_t i=0; i<vec1.GetNrows(); i++){
    if(vec2(i) != 0.){
      quot(i) = vec1(i) / vec2(i);
    }
    else{
      if(zero){
        quot(i) = 0.;
      }
      else{
        quot(i) = vec1(i);
      }
    }
  }
  return quot;
}

TMatrixD TUnfHisto::MatDivVec(const TMatrixD& mat, const TVectorD& vec, Int_t zero) const
{
  TMatrixD quotmat(_rDim, _gDim);
  for(Int_t i=0; i<_rDim; i++){
    for(Int_t j=0; j<_gDim; j++){
      if(vec(i) != 0.){
        quotmat(i,j) = mat(i,j) / vec(i);
      }
      else{
        if(zero){
          quotmat(i,j) = 0.;
        }
        else{
          quotmat(i,j) = mat(i,j);
        }
      }
    }
  }
  return quotmat;
}

TVectorD TUnfHisto::MatMultVec(const TMatrixD& mat, const TVectorD& vec, Int_t resdim) const
{
  TVectorD res(resdim);
  for(Int_t i=0; i<resdim; i++){
    for(Int_t j=0; j<vec.GetNrows(); j++){
      res(i) = res(i) + mat(i,j) * vec(j);
    }
  }
  return res;
}

TVectorD TUnfHisto::CompProd(const TVectorD& vec1, const TVectorD& vec2)
{
  TVectorD res(vec1.GetNrows());
  if(vec1.GetNrows() != vec2.GetNrows()){
    cout << "CompProd does not make sense for unequal dimensions of vectors.." << endl;
  }
  for(Int_t i=0; i<vec1.GetNrows(); i++){
    res(i) = vec1(i) * vec2(i);
  }
  return res;
}

void TUnfHisto::fillC(TMatrixD& tCurv, TMatrixD& tC, Int_t matdim, Int_t ndim) const
{
  Double_t eps = 0.00001;

  for(Int_t i=0; i<matdim; i++){
    for(Int_t j=0; j<matdim; j++){
      tC(i,j) = 0;
      tCurv(i,j) = 0;
    }
  }

  if(ndim == 0){
    for(Int_t i=0; i<matdim; i++){
      tC(i,i) = 1.;
    }
  }
  else if(ndim == 1){
    for(Int_t i=0; i<matdim; i++){
      if(i < matdim-1) {tC(i,i+1) = 1.0;}
      tC(i,i) = 1.0;
    }
  }
  else if(ndim == 2){
    for(Int_t i=0; i<matdim; i++){
      if(i > 0) {tC(i,i-1) = 1.0;}
      if(i < matdim-1) {tC(i,i+1) = 1.0;}
      tC(i,i) = -2.0;
    }
    tC(0,0) = -1.0;
    tC(matdim-1,matdim-1) = -1.0;
  }
  else if(ndim == 3){
    for(Int_t i=1; i<matdim-2; i++){
      tC(i,i-1) =  1.0;
      tC(i,i)   = -3.0;
      tC(i,i+1) =  3.0;
      tC(i,i+2) = -1.0;
    }
  }
  else if(ndim==4){
    for(Int_t i=0; i<matdim; i++){
      if(i > 0) {tC(i,i-1) = -4.0;}
      if(i < matdim-1) {tC(i,i+1) = -4.0;}
      if(i > 1) {tC(i,i-2) =  1.0;}
      if(i < matdim-2) {tC(i,i+2) =  1.0;}
      tC(i,i) = 6.0;
    }
    tC(0,0) = 2.0;
    tC(matdim-1,matdim-1) = 2.0;
    tC(0,1) = -3.0;
    tC(matdim-2,matdim-1) = -3.0;
    tC(1,0) = -3.0;
    tC(matdim-1,matdim-2) = -3.0;
    tC(1,1) =  6.0;
    tC(matdim-2,matdim-2) =  6.0;
  }
  else if(ndim == 5){
    for(Int_t i=2; i < matdim-3; i++){
      tC(i,i-2) = 1.0;
      tC(i,i-1) = -5.0;
      tC(i,i)   = 10.0;
      tC(i,i+1) = -10.0;
      tC(i,i+2) = 5.0;
      tC(i,i+3) = -1.0;
    }
  }
  else if(ndim == 6){
    for(Int_t i = 3; i < matdim - 3; i++){
      tC(i,i-3) = 1.0;
      tC(i,i-2) = -6.0;
      tC(i,i-1) = 15.0;
      tC(i,i)   = -20.0;
      tC(i,i+1) = 15.0;
      tC(i,i+2) = -6.0;
      tC(i,i+3) = 1.0;
    }
  }

  //Add epsilon to avoid singularities
  for(Int_t i=0; i<matdim; i++){
    tC(i,i) = tC(i,i) + eps;
  }

  //Get curvature matrix
  for(Int_t i=0; i<matdim; i++){
    for(Int_t j=0; j<matdim; j++){
      tCurv(i,j) = 0.0;
      for(Int_t k=0; k<matdim; k++){
        tCurv(i,j) = tCurv(i,j) + tC(k,i)*tC(k,j);
      }
    }
  }
}

Double_t TUnfHisto::GetCurvature(const TVectorD& vec, const TMatrixD& curv) const
{      
  Double_t cv(0);

  for(Int_t i=0; i<_gDim; i++){
    for(Int_t j=0; j<_gDim; j++){
      cv = cv + vec(i) * curv(i,j) * vec(j);
    }
  }
  return cv;
}

Double_t TUnfHisto::GetChi2(const TVectorD& vec, const TVectorD& refvec, const TVectorD& err, Int_t n) const
{
  Double_t chi2=0;
  Int_t ndof=0;

  for(Int_t i=0; i<n; i++){
    if(err(i)>0.){
      chi2 = chi2 + (vec(i) - refvec(i))*(vec(i) - refvec(i))/err(i)/err(i);
      ndof++;
    }
  }

  if(ndof>0){
    chi2 = chi2/ndof;
  }
  else{
    cout << "Chi2 calculation: No degrees of freedom" << endl;
  }

  return chi2;
}


void TUnfHisto::WriteResults(Int_t htau)
// Write the histos into a rootfile
{
  if (!_resfile) return;

  char name[200];

  //The histos that were read in 
  _b->Write();
  _xini->Write();
  _bini->Write();
  _A->Write();

  sprintf(name,"dd");
  if((TH1D*)gDirectory->Get(name)==NULL) cout << name << " is nil" << endl;
  ((TH1D*)gDirectory->Get(name))->Write();

  Int_t mintau, maxtau;
  mintau = htau-1; 
  maxtau = htau;
  if(htau==0){ //If htau==0 then loop over all singular values
    mintau = 0;
    maxtau = _gDim;
  }

  for(int j=mintau; j<maxtau; j++) {
    sprintf(name,"%s%d", "z", j+1);  
    ((TH1D*)gDirectory->Get(name))->Write();
    sprintf(name,"%s%d", "w", j+1);  
    ((TH1D*)gDirectory->Get(name))->Write();
    sprintf(name,"%s%d", "x", j+1);  
    ((TH1D*)gDirectory->Get(name))->Write();
  }

  sprintf(name,"sv");
  ((TH1D*)gDirectory->Get(name))->Write();

  _resfile->Close();
  delete _resfile; _resfile = 0;

}

void TUnfHisto::InitHistos(Int_t htau)
{
  TH1 *h;  

  char name[200], title[200];	

  //Delete the existing histos
  sprintf(name, "dd");  
  sprintf(title, "d vector after orthogonal transformation");  
  if((TH1D*)gDirectory->Get(name)!=NULL)
    ((TH1D*)gDirectory->Get(name))->~TH1D();

  Int_t mintau, maxtau;
  mintau = 0;
  maxtau = _gDim;

  for(int j=mintau; j<maxtau; j++) {
    sprintf(name,"%s%d", "z", j+1);  
    sprintf(title, "%s%d", "Regularized vector z with tau=", j+1);  
    if((TH1D*)gDirectory->Get(name)!=NULL)
      ((TH1D*)gDirectory->Get(name))->~TH1D();
    sprintf(name,"%s%d", "w", j+1);  
    sprintf(title, "%s%d", "Regularized weight w with tau=", j+1);  
    if((TH1D*)gDirectory->Get(name)!=NULL)
      ((TH1D*)gDirectory->Get(name))->~TH1D();
    sprintf(name,"%s%d", "x", j+1);  
    sprintf(title, "%s%d", "Regularized result x with tau=", j+1);  
    if((TH1D*)gDirectory->Get(name)!=NULL)
      ((TH1D*)gDirectory->Get(name))->~TH1D();
  }

  sprintf(name, "sv");  
  sprintf(title, "Singular values of AC^-1");  
  if((TH1D*)gDirectory->Get(name)!=NULL)
    ((TH1D*)gDirectory->Get(name))->~TH1D();

  //Now make the needed new histos
  sprintf(name, "dd");  
  sprintf(title, "d vector after orthogonal transformation");  
  h = new TH1D(name, title, _b->GetNbinsX(), 0, _b->GetNbinsX());  
  h->Sumw2();

  mintau = htau-1; 
  maxtau = htau;
  if(htau==0){ //If htau==0 then loop over all singular values
    mintau = 0;
    maxtau = _gDim;
  }

  for(int j=mintau; j<maxtau; j++) {
    sprintf(name,"%s%d", "z", j+1);  
    sprintf(title, "%s%d", "Regularized vector z with tau=", j+1);  
    h = new TH1D(name, title, _xini->GetNbinsX(), 0, _xini->GetNbinsX());  
    h->Sumw2();
    sprintf(name,"%s%d", "w", j+1);  
    sprintf(title, "%s%d", "Regularized weight w with tau=", j+1);  
    h = new TH1D(name, title, _xini->GetNbinsX(), _xini->GetXaxis()->GetXmin(), _xini->GetXaxis()->GetXmax());  
    h->Sumw2();
    sprintf(name,"%s%d", "x", j+1);  
    sprintf(title, "%s%d", "Regularized result x with tau=", j+1);  
    h = new TH1D(name, title, _xini->GetNbinsX(), _xini->GetXaxis()->GetXmin(), _xini->GetXaxis()->GetXmax());  
    h->Sumw2();
  }

  sprintf(name, "sv");  
  sprintf(title, "Singular values of AC^-1");  
  h = new TH1D(name, title, _gDim, 0, _gDim);  
  h->Sumw2();

  return;
}

void TUnfHisto::V2H(const TVectorD& vec, TString hname)
{
  TH1D *histo;
  histo = ((TH1D*)gDirectory->Get(hname));
  for(Int_t i=0; i<vec.GetNrows(); i++){
    histo->SetBinContent(i+1, vec(i));
  }
  return;
}


