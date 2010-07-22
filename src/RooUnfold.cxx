//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.cxx,v 1.18 2010-07-22 10:30:25 fwx38934 Exp $
//
// Description:
//      Unfolding framework base class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfold.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "TMatrixD.h"
#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "RooUnfoldResponse.h"

// Need subclasses just for RooUnfold::New()
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;
using std::sqrt;
using std::fabs;

ClassImp (RooUnfold);

RooUnfold*
RooUnfold::New (Algorithm alg, const RooUnfoldResponse* res, const TH1* meas, Int_t regparm,
                const char* name, const char* title)
{
  RooUnfold* unfold;
  switch (alg) {
    case kNone:
      unfold= new RooUnfold         (res, meas);   // dummy unfolding
      break;
    case kBayes:
      unfold= new RooUnfoldBayes    (res, meas, regparm);
      break;
    case kSVD:
      unfold= new RooUnfoldSvd      (res, meas, regparm);
      break;
    case kBinByBin:
      unfold= new RooUnfoldBinByBin (res, meas);
      break;
    default:
      cerr << "Unknown RooUnfold method " << Int_t(alg) << endl;
      return 0;
  }
  if (name || title) unfold->SetNameTitle (name, title);
  return unfold;
}

RooUnfold*
RooUnfold::Clone (const char* newname) const
{
  RooUnfold* unfold= new RooUnfold(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

RooUnfold::RooUnfold (const RooUnfold& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  Init();
  CopyData (rhs);
}

RooUnfold::RooUnfold (const RooUnfoldResponse* res, const TH1* meas, const char* name, const char* title)
  : TNamed (name, title)
{
  Init();
  Setup (res, meas);
}

void
RooUnfold::Assign (const RooUnfold& rhs)
{
  if (this == &rhs) return;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  CopyData (rhs);
}

void
RooUnfold::CopyData (const RooUnfold& rhs)
{
  Setup (rhs.response(), rhs.Hmeasured());
}

void
RooUnfold::Reset()
{
  Init();
}

void
RooUnfold::Init()
{
  _res= 0;
  _meas= 0;
  _nm= _nt= 0;
  _verbose= 1;
  _overflow= 0;
  _unfolded= _haveCov= _fail= false;
}

RooUnfold&
RooUnfold::Setup (const RooUnfoldResponse* res, const TH1* meas)
{
  Reset();
  _res= res;
  _meas= meas;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  SetNameTitleDefault();
  return *this;
}

// Dummy unfolding - just copies input
void
RooUnfold::Unfold()
{
  cout << "********************** " << ClassName() << ": dummy unfolding - just copy input **********************" << endl;
  Int_t first= _overflow ? 0 : 1;
  _rec.ResizeTo (_nt + (_overflow ? 2 : 0));
  Int_t nb= _nm < _nt ? _nm : _nt;
  if (_overflow) nb += 2;
  for (Int_t i= 0; i < nb; i++) {
    _rec[i]= _meas->GetBinContent (i+first);
  }
  _unfolded= true;
}

// Dummy covariance matrix - just uses input
void
RooUnfold::GetCov()
{
  Int_t first= _overflow ? 0 : 1;
  Int_t nt= _nt + (_overflow ? 2 : 0);
  _cov.ResizeTo (nt, nt);
  Int_t nb= _nm < _nt ? _nm : _nt;
  if (_overflow) nb += 2;
  for (Int_t i= 0; i < nb; i++) {
    Double_t err= _meas->GetBinError (i+first);
    _cov(i,i)= err*err;
  }
  _haveCov= true;
}


Double_t RooUnfold::Chi2(const TH1* hTrue)
{
	const TH1* hReco=Hreco ();
	Int_t nt= _nt+(_overflow ? 2 : 0);
	TMatrixD reco_matrix(nt,1);
  	for (Int_t i = 0 ; i < nt; i++) {
    	Int_t it= RooUnfoldResponse::GetBin (hReco, i, _overflow);
    	if ((hReco->GetBinContent(it)!=0.0 || (hReco->GetBinError(it)>0.0)) &&
            (hTrue->GetBinContent(it)!=0.0 || (hTrue->GetBinError(it)>0.0))) {
           	reco_matrix(i,0)    = hReco->GetBinContent(it) - hTrue->GetBinContent(it);
        }
  	}
  	TMatrixD Ereco_copy=Ereco();
	TMatrixD Ereco_copy2=Ereco();
  	Double_t Ereco_det = Ereco_copy.Determinant();
  	TMatrixD reco_matrix_t=reco_matrix;
  	TMatrixD reco_matrix_copy=reco_matrix;
	reco_matrix_t.T();
	if (fabs(Ereco_det)<1e-5){
		cerr << "Warning: Small Determinant of Covariance Matrix =" << Ereco_det << endl;
		cerr << "Chi^2 may be invalid due to small determinant" << endl;
	}
	TDecompSVD svd(Ereco_copy);
	svd.MultiSolve(reco_matrix);
	Double_t cond=svd.Condition();
	TMatrixD chisq_nw = reco_matrix_t*reco_matrix;
	
	double cond_max=1e17;
	if (cond>=cond_max){
		cerr << "Warning, very large matrix condition= "<< cond<<" chi^2 may be inaccurate"<<endl;}
	return chisq_nw(0,0);	
}


void
RooUnfold::PrintTable (std::ostream& o, const TH1* hTrue, Bool_t withError)
{
  const TH1* hReco=      Hreco (withError);
  if (_fail) return;
  const TH1* hMeas=      Hmeasured();
  const TH1* hTrainTrue= response()->Htruth();
  const TH1* hTrain=     response()->Hmeasured();

  std::ostringstream fmt;
  fmt.copyfmt (o);   // save original ostream format
  Int_t dim= hReco->GetDimension(), ntxb= hReco->GetNbinsX()+2, ntyb= hReco->GetNbinsY()+2;
  if (hMeas->GetDimension() != dim || hMeas->GetNbinsX()+2 != ntxb || hMeas->GetNbinsY()+2 != ntyb) dim= 1;
  Int_t iwid= (dim==3) ? 8 : (dim==2) ? 7 : 5;
  const char* xwid= (dim==3) ? "===" : (dim==2) ? "==" : "";
  o << "====================================================================" << xwid << endl
    << setw(iwid) << ""      << setw(9) << "Train" << setw(9) << "Train"    << setw(9) << "Test"  << setw(9) << "Test"  << setw(9) << "Unfolded" << setw(9) << "Diff" << setw(9) << "Pull" << endl
    << setw(iwid) << "Bin"   << setw(9) << "Truth" << setw(9) << "Measured" << setw(9) << "Truth" << setw(9) << "Input" << setw(9) << "Output"   << endl
    << "====================================================================" << xwid << endl;

  Double_t chi2= 0.0;
  Int_t ndf= 0, first= (_overflow ? 0 : 1);
  Int_t overflow= (_overflow ? 2 : 0), nt= _nt+overflow, nm= _nm+overflow;
  Int_t maxbin= nt < nm ? nm : nt;
  for (Int_t i = 0 ; i < maxbin; i++) {
    Int_t it= RooUnfoldResponse::GetBin (hReco, i, _overflow);
    Int_t im= RooUnfoldResponse::GetBin (hMeas, i, _overflow);
	
    if (dim==2 || dim==3) {
      // ROOT 5.26 has GetBinXYZ to do this.
      Int_t iw= (dim==2) ? 3 : 2;
      Int_t ix= it%ntxb;
      Int_t iy= ((it-ix)/ntxb)%ntyb;
      o << setw(iw) << ix << ',' << setw(iw) << iy;
      if (dim==3) o << ',' << setw(iw) << ((it-ix)/ntxb - iy)/ntyb;
    } else
      o << setw(iwid) << i+first;
    o << std::fixed << setprecision(0);
    if (i<nt)
      o << ' ' << setw(8) << hTrainTrue->GetBinContent(it);
    else
      o << setw(9) << " ";
    if (i<nm)
      o << ' ' << setw(8) << hTrain->GetBinContent(im);
    else
      o << setw(9) << " ";
    if (hTrue && i<nt)
      o << ' ' << setw(8) << hTrue->GetBinContent(it);
    else
      o << setw(9) << " ";
    if (i<nm)
      o << ' ' << setw(8) << hMeas->GetBinContent(im);
    else
      o << setw(9) << " ";
    o << setprecision(1);
    if (i<nt) {
      o << ' ' << setw(8) << hReco->GetBinContent(it);
      if (hTrue &&
          ((hReco->GetBinContent(it)!=0.0 || (withError && hReco->GetBinError(it)>0.0)) &&
           (hTrue->GetBinContent(it)!=0.0 || (withError && hTrue->GetBinError(it)>0.0)))) {
        Double_t ydiff    = hReco->GetBinContent(it) - hTrue->GetBinContent(it);
        Double_t ydiffErr = hReco->GetBinError(it);
        o << ' ' << setw(8) << ydiff;
        if (ydiffErr>0.0) {
          Double_t ypull = ydiff/ydiffErr;
          chi2 += ypull*ypull;//calculates chi^2 value without matrix method. Only valid for diagonal covariance matrices
          ndf++;
          o << setprecision(3) << ' ' << setw(8) << ypull;
        }
      }
    }
    o << endl;
  }
  
  o << "====================================================================" << xwid << endl
    << setw(iwid) << "" << std::fixed << setprecision(0)
    << ' ' << setw(8) << hTrainTrue->Integral()
    << ' ' << setw(8) << hTrain->Integral();
  if (hTrue)
    o << ' ' << setw(8) << hTrue->Integral();
  else
    o << setw(9) << " ";
  o << ' ' << setw(8) << hMeas->Integral() << setprecision(1)
    << ' ' << setw(8) << hReco->Integral()
    << endl
    << "====================================================================" << xwid << endl;
  o.copyfmt (fmt);  // restore original ostream format
  Double_t chi_squ;
  if (withError){
    chi_squ = Chi2(hTrue);
  }
  else{
    chi_squ = chi2;
  }
  o << "Chi^2/NDF=" << chi_squ << "/" << ndf;
  if (withError) o << " (bin-by-bin Chi^2=" << chi2 << ")";
  o << endl;
  if (chi_squ<=0){
  	cerr << "Warning: Invalid Chi^2 Value" << endl;
  }
}

void
RooUnfold::SetNameTitleDefault()
{
  if (!_res) return;
  const char* s= GetName();
  if (s[0] == '\0') SetName  (_res->GetName());
  s= GetTitle();
  if (s[0] == '\0') {
    TString title= "Unfold ";
    title += _res->GetTitle();
    SetTitle (title);
  }
}

TH1*
RooUnfold::Hreco (Bool_t withError)
{
  TH1* reco= (TH1*) _res->Htruth()->Clone(GetName());
  reco->Reset();
  reco->SetTitle (GetTitle());
  if (!_unfolded)             Unfold();
  if (_fail)                  return 0;
  if (withError && !_haveCov) GetCov();
  if (!_haveCov) withError= false;
  Int_t nt= _nt + (_overflow ? 2 : 0);
  for (Int_t i= 0; i < nt; i++) {
    Int_t j= RooUnfoldResponse::GetBin (reco, i, _overflow);
    reco->SetBinContent (j,             _rec(i));
    if (withError)
      reco->SetBinError (j, sqrt (fabs (_cov(i,i))));
  }
  return reco;
}
