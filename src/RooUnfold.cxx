//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.cxx,v 1.22 2010-08-04 22:05:31 adye Exp $
//
// Description:
//      Unfolding framework base class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p> Main class for unfolding of distributions.</p>
<p> Permits the use of several different unfolding techniques which can be initiated using the Algorithm parameter, 
returns a table of values for each bin giving the true, measured and reconstructed values of a distribution and a chi squared value. 
Also sets the errors on a reconstructed distribution (Hreco())</p>
<p> The errors in both chi squared and in the reconstructed distribution can both be calculated in one of 3 ways. 
These can be set using their respective parameters.</p>
<p> The calculation of chi squared can be done with a simple sum of the sum of the squares of the residuals/error,
 or using a covariance matrix which can be set either using the error matrix from the unfolding or using the spread of the reconstructed values. 
 The errors on the reconstructed distribution are calculated in the same way except that the simplest option is simply to have no errors. </p>  
END_HTML */

/////////////////////////////////////////////////////////////

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
#include "RooUnfoldAll.h"

// Need subclasses just for RooUnfold::New()
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBinByBin.h"
#ifndef NOTUNFOLD
#include "RooUnfoldTUnfold.h"
#endif

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
	/*Unfolds according to the value of alg:
	0: a dummy unfold
	1: Unfold via a Bayes method
	2: Unfold using singlar value decomposition
	3: Unfold bin by bin.
	*/
  RooUnfold* unfold;
  switch (alg) {
    case kNone:
      unfold= new RooUnfold         (res, meas); 
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
#ifndef NOTUNFOLD
    case kTUnfold:
      unfold= new RooUnfoldTUnfold (res,meas,regparm);
      break;
#endif
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
	//Creates a copy of the unfold object
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
	//Calls Init()
  Init();
}

void
RooUnfold::Init()
{
	//Sets al variables to 0
  _res= 0;
  _meas= 0;
  _nm= _nt= 0;
  _verbose= 1;
  _overflow= 0;
  _unfolded= _haveCov= _fail=_have_err_mat= false;
  _Nits=50;
}

RooUnfold&
RooUnfold::Setup (const RooUnfoldResponse* res, const TH1* meas)
{
	//Sets parameters
  Reset();
  _res= res;
  _meas= meas;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  SetNameTitleDefault();
  Get_settings();
  return *this;
}


void
RooUnfold::Unfold()
{
	// Dummy unfolding - just copies input
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

void
RooUnfold::GetCov()
{
	//Creates covariance matrix using bin error on measured distribution.
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

void
RooUnfold::Get_err_mat()
{
	//Get error matrix based on residuals (uses RooUnfoldAll)
  Int_t nt= _nt;
  _err_mat.ResizeTo(nt,nt);
  RooUnfoldAll All(_Nits,this);
  All.Plotting();
  _err_mat=All.True_err();
  _have_err_mat= true;
}
Double_t RooUnfold::Chi2(const TH1* hTrue,Int_t DoChi2)
{
	/*Calculates Chi squared. Method depends on value of DoChi2
	0: sum of (residuals/error)squared
	1: use covariance matrix returned from unfolding
	2: use covariance matrix returned from Get_err_mat() 
	Returns warnings for small determinants of covariance matrices and if the condition is very large.*/
	const TH1* hReco=Hreco (DoChi2);
	Int_t nt= _nt+(_overflow ? 2 : 0);
	if (DoChi2==1||DoChi2==2){
		TMatrixD reco_matrix(nt,1);
	  	for (Int_t i = 0 ; i < nt; i++) {
	    	Int_t it= RooUnfoldResponse::GetBin (hReco, i, _overflow);
	    	if ((hReco->GetBinContent(it)!=0.0 || (hReco->GetBinError(it)>0.0)) &&
	            (hTrue->GetBinContent(it)!=0.0 || (hTrue->GetBinError(it)>0.0))) {
	           	reco_matrix(i,0)    = hReco->GetBinContent(it) - hTrue->GetBinContent(it);
	        }
	  	}
	  	
	  	TMatrixD Ereco_copy(nt,nt);
	  	if (DoChi2==1){
	  	Ereco_copy=Ereco();
	  	}
	  	if (DoChi2==2){
	  	Ereco_copy=Freco();
	  	}
		TMatrixD Ereco_copy2=Ereco_copy;
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
			cerr << "Warning, very large matrix condition= "<< cond<<" chi^2 may be inaccurate"<<endl;
		}
		return chisq_nw(0,0);	
	}
	else{
		Double_t chi2=0;
		for (Int_t i = 0 ; i < nt; i++) {
	    	Int_t it= RooUnfoldResponse::GetBin (hReco, i, _overflow);
			if((hReco->GetBinContent(it)!=0.0 || (hReco->GetBinError(it)>0.0)) &&
           (hTrue->GetBinContent(it)!=0.0 || (hTrue->GetBinError(it)>0.0))) {
       			Double_t ydiff    = hReco->GetBinContent(it) - hTrue->GetBinContent(it);
        		Double_t ydiffErr = hReco->GetBinError(it);
        		if (ydiffErr>0.0) {
	          		Double_t ypull = ydiff/ydiffErr;
    	      		 chi2 += ypull*ypull;
        		}
        	}
		}
        return chi2;
	}
}


void
RooUnfold::PrintTable (std::ostream& o, const TH1* hTrue, Int_t withError)
{
	//Prints data from truth, measured and reconstructed data for each bin. 
  const TH1* hReco=      Hreco (withError);
  if (_fail) return;
  const TH1* hMeas=      Hmeasured();
  const TH1* hTrainTrue= response()->Htruth();
  const TH1* hTrain=     response()->Hmeasured();

  std::ostringstream fmt;
  fmt.copyfmt (o);  
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
          ndf++;
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
  o.copyfmt (fmt);
  Double_t chi_squ;
  chi_squ = Chi2(hTrue,withError);
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
RooUnfold::Hreco (Int_t withError)
{
	/*Creates reconstructed distribution. Error calculation varies by withError:
	0: No errors
	1: Errors from the square root of the covariance matrix given by the unfolding
	2: Errors from the square root of the covariance matrix given by Get_err_mat()
	*/

  TH1* reco= (TH1*) _res->Htruth()->Clone(GetName());
  reco->Reset();
  reco->SetTitle (GetTitle());
  if (!_unfolded)             Unfold();
  if (_fail)                  return 0;
  if (withError==1 && !_haveCov) GetCov();
  if (!_haveCov) withError= 0;
  if (!_have_err_mat && withError==2) Get_err_mat();
  if (!_err_mat.GetNoElements() && withError==2){cout<<"it breaks here..."<<endl;}
  Int_t nt= _nt + (_overflow ? 2 : 0);
  for (Int_t i= 0; i < nt; i++) {
    Int_t j= RooUnfoldResponse::GetBin (reco, i, _overflow);
    reco->SetBinContent (j,             _rec(i));
    if (withError==1){
      reco->SetBinError (j, sqrt (fabs (_cov(i,i))));
  	}
  	if (withError==2){
  		reco->SetBinError(j,sqrt(fabs(_err_mat(i,i))));
  	}
  }
  return reco;
}

void
RooUnfold::Get_settings(){
	//Gets maximum and minimum parameters and step size
	_minparm=0;
	_maxparm=0;
	_stepsizeparm=0;
	_defaultparm=0;
}

Double_t
RooUnfold::Get_minparm(){
	return _minparm;
}

Double_t
RooUnfold::Get_maxparm(){
	return _maxparm;
}

Double_t
RooUnfold::Get_stepsizeparm(){
	return _stepsizeparm;
}

Double_t
RooUnfold::Get_defaultparm(){
	return _defaultparm;
}

