//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfold.cxx,v 1.41 2010-08-25 22:30:53 adye Exp $
//
// Description:
//      Unfolding framework base class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>A base class for several unfolding methods.
<p>The unfolding method can either use the constructors for individual unfolding algorithms or the New() method, specifiying the algorithm to be used.
<p>The resultant distribution can be displayed as a plot (Hreco) or as a bin by bin breakdown of the true, measured and reconstructed values (PrintTable)
<p>A covariance matrix can be returned using the Ereco() method. A vector of its diagonals can be returned with the ErecoV() method.
<p>A summary of the unfolding algorithms which inherit from this class is below: 
<ul>
<li>RooUnfoldBayes: Uses the Bayes method of unfolding based on the method written by D'Agostini (<a href="http://www.slac.stanford.edu/spires/find/hep/www?j=NUIMA,A362,487">NIM A 362 (1995) 487</a>).
<ul>
<li>Works for 1 & 2 dimensional distributions
<li>Returned errors can be either as a diagonal matrix or as a full matrix of covariances
<li>Regularisation parameter sets the number of iterations used in the unfolding (default=4)
<li>Is able to account for bin migration and smearing
<li>Can unfold if test and measured distributions have different binning. 
<li>Returns covariance matrices with conditions approximately that of the machine precision. This occasionally leads to very large chi squared values
</ul>
<li> RooUnfoldSVD: Uses the singular value decomposition method of Hocker and Kartvelishvili (<a href="http://arxiv.org/abs/hep-ph/9509307">NIM A 372 (1996) 469</a>)
<ul>
<li>Regularisation parameter defines the level at which values are deemed to be due to statistical fluctuations and are cut out. (Default= number of bins/2)
<li>Returns errors as a full matrix of covariances
<li>Error processing is much the same as with the kCovToy setting with 1000 toys. This is quite slow but can be switched off.
<li>Can only handle 1 dimensional distributions
<li>True and measured distributions must have the same binning
<li>Can account for both smearing and biasing
<li>Returns near singular covariance matrices, again leading to very large chi squared values
</ul>
<li> RooUnfoldBinByBin: Unfolds using the method of correction factors.
<ul>
<li>Returns errors as a diagonal matrix. 
<li>Is not able to handle bin migration caused by bias/smearing of the distribution
<li>Can only handle 1 dimensional distributions
<li>True and measured distributions must have the same binning
</ul>
<li> RooUnfoldTUnfold: Uses the unfolding method implemented in ROOT's <a href="http://root.cern.ch/root/html/TUnfold.html">TUnfold</a> class
<ul>
<li>Only included in ROOT versions 5.22 and higher 
<li>Only able to reconstruct 1 dimensional distributions
<li>Can account for bin migration and smearing
<li>Errors come as a full covariance matrix. 
<li>Will sometimes warn of "unlinked" bins. These are bins with 0 entries and do not effect the results of the unfolding
<li>Regularisation parameter can be either optimised internally by plotting log10(chi2 squared) against log10(tau). The 'kink' in this curve is deemed the optimum tau value. This value can also be set manually (FixTau)
<li>The latest version (TUnfold v15) requires that RooUnfoldResponse::SetOverflow=0. ROOT versions 5.26 or below use v13 and so should be safe to use overflows
</ul>
<li> RooUnfoldInvert: The simplest method of unfolding works by simply inverting the response matrix. 
<ul>
<li>For small statistics, this method does not produce useful results. 
<li>The inversion method is included largely to illustrate the necessity of a more effective method of unfolding</ul>
</ul>     
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfold.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>

#include "TMatrixD.h"
#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldErrors.h"
#include "TRandom.h"
// Need subclasses just for RooUnfold::New()
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldInvert.h"
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
RooUnfold::New (Algorithm alg, const RooUnfoldResponse* res, const TH1* meas,Double_t regparm,
                const char* name, const char* title)
{
    /*Unfolds according to the value of the alg enum:
    0: a dummy unfold
    1: Unfold via a Bayes method
    2: Unfold using singlar value decomposition
    3: Unfold bin by bin.
    4: Unfold with TUnfold
    5: Unfold using inversion of response matrix
    */
  RooUnfold* unfold;
  switch (alg) {
    case kNone:
      unfold= new RooUnfold         (res, meas); 
      break;
    case kBayes:
      unfold= new RooUnfoldBayes    (res, meas);
      break;
    case kSVD:
      unfold= new RooUnfoldSvd      (res, meas);
      break;
    case kBinByBin:
      unfold= new RooUnfoldBinByBin (res, meas);
      break;
    case kTUnfold:
#ifndef NOTUNFOLD
      unfold= new RooUnfoldTUnfold  (res,meas);
      break;
#else
      cerr << "TUnfold library is not available" << endl;
      return 0;
#endif
    case kInvert:
      unfold = new RooUnfoldInvert (res,meas);
      break;
    default:
      cerr << "Unknown RooUnfold method " << Int_t(alg) << endl;
      return 0;
  }
  if (name || title) unfold->SetNameTitle (name, title);
  if (regparm != -1e30){
    unfold->SetRegParm(regparm);
  }
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
  SetVerbose (rhs.verbose());
  SetNToys   (rhs.NToys());
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
  _unfolded= _haveCov= _fail=_have_err_mat=_haveErrors= false;
  _NToys=50;
  GetSettings();
}

RooUnfold&
RooUnfold::Setup (const RooUnfoldResponse* res, const TH1* meas)
{
  Reset();
  SetResponse (res);
  SetMeasured (meas);
  return *this;
}

void
RooUnfold::SetResponse (const RooUnfoldResponse* res)
{
  _res= res;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  SetNameTitleDefault();
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
RooUnfold::GetErrors()
{
    //Creates vector of diagonals of covariance matrices
    if (!_unfolded) Unfold();
    if (_fail)      return;
    Int_t nt= _nt + (_overflow ? 2 : 0);
    _errors.ResizeTo(nt);
    if (!_haveCov) GetCov();
    for (Int_t i= 0; i < nt; i++) {
      _errors(i)= _cov(i,i);
    }
    _haveErrors= true;
}

void
RooUnfold::GetCov()
{
    //Creates covariance matrix using bin error on measured distribution.
  if (!_unfolded) Unfold();
  if (_fail)      return;
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
RooUnfold::GetErrMat()
{
    //Get error matrix based on residuals found using the Runtoy method. 
    if (!_unfolded) Unfold();
    if (_fail)      return;
    Int_t total= _nt + (_overflow ? 2 : 0);
    vector<double> bc_vec(total);
    TMatrixD bc_mat(total,total);
    _err_mat.ResizeTo(total,total);
    for (int k=0; k<_NToys; k++){
        TH1* res=this->Runtoy();
        for (int i=0; i<total;i++){
            Double_t res_bci=res->GetBinContent(i);
            bc_vec[i]+=res_bci;
            for (int j=0; j<total; j++){
            bc_mat(i,j)+=(res_bci*res->GetBinContent(j));
            }
        }
    }
    for (unsigned int i=0; i<bc_vec.size();i++){
        for (unsigned int j=0; j<bc_vec.size(); j++){
            _err_mat(i,j)=(bc_mat(i,j)-(bc_vec[i]*bc_vec[j])/_NToys)/_NToys;
        }
    }
    _have_err_mat=true;
}

Double_t RooUnfold::Chi2(const TH1* hTrue,ErrorTreatment DoChi2)
{
    /*Calculates Chi squared. Method depends on value of DoChi2
    0: sum of (residuals/error)squared
    1: use errors propagated through the unfolding
    2: use covariance matrix returned from unfolding
    3: use covariance matrix returned from GetErrMat() 
    Returns warnings for small determinants of covariance matrices and if the condition is very large.
    If a matrix has to be inverted also removes rows/cols with all their elements equal to 0*/

    const TH1* hReco=Hreco (DoChi2);
    Int_t nt= _nt+(_overflow ? 2 : 0);
    TMatrixD Ereco_copy;
    if (DoChi2==kCovariance||DoChi2==kCovToy){
        Ereco_copy.ResizeTo(Ereco(DoChi2).GetNrows(),Ereco(DoChi2).GetNcols());
        Ereco_copy=Ereco(DoChi2);
        nt=Ereco(DoChi2).GetNrows();
        
        TMatrixD reco_matrix(nt,1);
        for (Int_t i = 0 ; i < nt; i++) {
            Int_t it= RooUnfoldResponse::GetBin (hReco, i, _overflow);
            if ((hReco->GetBinContent(it)!=0.0 || (hReco->GetBinError(it)>0.0)) &&
                (hTrue->GetBinContent(it)!=0.0 || (hTrue->GetBinError(it)>0.0))) {
                reco_matrix(i,0)    = hReco->GetBinContent(it) - hTrue->GetBinContent(it);
            }
        }
        TMatrixD Ereco_copy2=Ereco_copy;
        //cutting out 0 elements    
        TMatrixD Ereco_copy_cut=CutZeros(Ereco_copy);
        TMatrixD reco_matrix_cut(Ereco_copy_cut.GetNrows(),1);
        int v=0;
        for (int i=0;i<Ereco_copy.GetNrows();i++){
            if(Ereco_copy(i,i)==0){
                v++;
            }
            else{
                reco_matrix_cut(i-v,0)=reco_matrix(i,0);
            }
        }
        int zeros;
        for (int i=0; i<Ereco_copy_cut.GetNrows(); i++){
            if (Ereco_copy (i,i) ==0){
                zeros++;
            }
        }
        Double_t Ereco_det=Ereco_copy_cut.Determinant();
        TMatrixD reco_matrix_t=reco_matrix_cut;
        reco_matrix_t.T();
        if (fabs(Ereco_det)<1e-5 && _verbose>=1){
            cerr << "Warning: Small Determinant of Covariance Matrix =" << Ereco_det << endl;
            cerr << "Chi^2 may be invalid due to small determinant" << endl;
        }
        TDecompSVD svd(Ereco_copy_cut);
        Double_t cond=svd.Condition();
        if (_verbose>=1){
            cout<<"For Covariance matrix condition= "<<cond<<" determinant= "<<Ereco_det<<endl;
        }
        svd.MultiSolve(reco_matrix_cut);
        TMatrixD chisq_nw = reco_matrix_t*reco_matrix_cut;
        double cond_max=1e17;
        if (cond>=cond_max && _verbose>=1){
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
RooUnfold::PrintTable (std::ostream& o, const TH1* hTrue, ErrorTreatment withError)
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
  Double_t true_train_tot=0;
  Double_t meas_train_tot=0;
  Double_t true_test_tot=0;
  Double_t meas_test_tot=0;
  Double_t unf_tot=0;
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
    true_train_tot+=hTrainTrue->GetBinContent(it);
    meas_train_tot+=hTrain->GetBinContent(im);
    true_test_tot+=hTrue->GetBinContent(it);
    meas_test_tot+=hMeas->GetBinContent(im);
    unf_tot+=hReco->GetBinContent(it);
    if (i<nt){
      o << ' ' << setw(8) << hTrainTrue->GetBinContent(it);
    }
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
          Double_t ypull = ydiff/ydiffErr;
          chi2 += ypull*ypull;
           o << ' ' << setw(8) << ypull;
        }
      }
    }
    o << endl;
    
  }
  
  o << "====================================================================" << xwid << endl
    << setw(iwid) << "" << std::fixed << setprecision(0)
    << ' ' << setw(8) << true_train_tot
    << ' ' << setw(8) << meas_train_tot;
  if (hTrue)
    o << ' ' << setw(8) << true_test_tot;
  else
    o << setw(9) << " ";
  o << ' ' << setw(8) << meas_test_tot << setprecision(1)
    << ' ' << setw(8) << unf_tot
    << ' ' << setw(8) << unf_tot-true_test_tot;
    if(hMeas->Integral()>0){
    o<< ' ' << setw(8) <<(unf_tot-true_test_tot)/sqrt(meas_test_tot);
    }
    o<< endl
    << "====================================================================" << xwid << endl;
  o.copyfmt (fmt);
  Double_t chi_squ = chi2;
  if (withError==kCovariance || withError==kCovToy) {
    chi_squ = Chi2(hTrue,withError);
    o << "Chi^2/NDF=" << chi_squ << "/" << ndf << " (bin-by-bin Chi^2=" << chi2 << ")";
  } else {
    o << "Bin-by-bin Chi^2/NDF=" << chi_squ << "/" << ndf;
  }
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
RooUnfold::Hreco (ErrorTreatment withError)
{
    /*Creates reconstructed distribution. Error calculation varies by withError:
    0: No errors
    1: Errors from the square root of the diagonals of the covariance matrix given by the unfolding
    2: Errors from the square root of of the covariance matrix given by the unfolding
    3: Errors from the square root of the covariance matrix given by GetErrMat()
    */
  TH1* reco= (TH1*) _res->Htruth()->Clone(GetName());
  reco->Reset();
  reco->SetTitle (GetTitle());
  if (!_unfolded) Unfold();
  if (_fail)      return 0;
  if (withError==kErrors && !_haveErrors){
    GetErrors();
    if (!_haveErrors) withError= kNoError;
  }
  if (withError==kCovariance && !_haveCov) {
    GetCov();
    if (!_haveCov) withError= kNoError;
  }
  
  if (!_have_err_mat && withError==kCovToy) GetErrMat();
  Int_t nt= _nt + (_overflow ? 2 : 0);
  for (Int_t i= 0; i < nt; i++) {
    Int_t j= RooUnfoldResponse::GetBin (reco, i, _overflow);
    reco->SetBinContent (j,             _rec(i));
    if (withError==kErrors){
      reco->SetBinError (j, sqrt (fabs (_errors(i))));  
    }
    if (withError==kCovariance){
        reco->SetBinError(j,sqrt(fabs(_cov(i,i))));
    }
    if (withError==kCovToy){
        reco->SetBinError(j,sqrt(fabs(_err_mat(i,i))));
    }
  }
  
  return reco;
}

void
RooUnfold::GetSettings(){
    //Gets maximum and minimum parameters and step size
    _minparm=0;
    _maxparm=0;
    _stepsizeparm=0;
    _defaultparm=0;
}

Double_t
RooUnfold::GetMinParm() const{
    //Get minimum regularisation parameter for unfolding method
    return _minparm;
}

Double_t
RooUnfold::GetMaxParm() const{
    //Get maximum regularisation parameter for unfolding method
    return _maxparm;
}

Double_t
RooUnfold::GetStepSizeParm() const{
    //Get suggested step size for unfolding distribution
    return _stepsizeparm;
}

Double_t
RooUnfold::GetDefaultParm() const{
    //Get suggested regularisation parameter.
    return _defaultparm;
}

TH1*
RooUnfold::Runtoy(ErrorTreatment witherror,double* chi2, const TH1* hTrue) const{
    /*
    Returns unfolded distribution for one iteration of unfolding. Use multiple toys to find residuals
    Can also return chi squared if a truth distribution is available. 
     */
    RooUnfold* unfold_copy = Clone("unfold_toy");
    const TH1* hMeas_AR = Hmeasured();
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    TH1* hMeas=Add_Random(hMeas_AR);
    TH1::AddDirectory (oldstat);
    unfold_copy->SetMeasured(hMeas);
    TH1* hReco;
    if (witherror==kCovToy){
        hReco= unfold_copy->Hreco(kNoError);
    }
    else{
        hReco= unfold_copy->Hreco(witherror);
    }
    if (chi2 && !hTrue){
        cerr<<"Error: can't calculate chi^2 without a truth distribution"<<endl;
    }
    if (chi2 && hTrue) {*chi2=unfold_copy->Chi2(hTrue,witherror);}
    delete hMeas;
    return hReco;
}


TH1*
RooUnfold::Add_Random(const TH1* hMeas_AR)
{
    //Adds a random number to the measured distribution before the unfolding//
    TH1* hMeas2=   dynamic_cast<TH1*>(hMeas_AR->Clone());
    for (Int_t i=0; i<hMeas_AR->GetNbinsX()+2 ; i++){
            Double_t err=hMeas_AR->GetBinError(i);
            Double_t new_x = hMeas_AR->GetBinContent(i) + gRandom->Gaus(0,err);
            hMeas2->SetBinContent(i,new_x);
            hMeas2->SetBinError(i,err);
    }
    return hMeas2;
}

void 
RooUnfold::Print(Option_t *opt)const
{
    TNamed::Print(opt);
    cout <<"regularisation parameter = "<<this->GetRegParm()<<endl;
    cout <<"ntoys = "<<this->NToys()<<endl;
}

TMatrixD
RooUnfold::CutZeros(const TMatrixD& Ereco_copy)
{
    //Removes row & column if all their elements are 0. 
    vector<int> diags;
        int missed=0;
        for (int i=0; i<Ereco_copy.GetNrows(); i++){
            double coltot=0;
            for (int j=0;j<Ereco_copy.GetNrows();j++){
                coltot+=Ereco_copy(i,j);
            }
            if (coltot==0){
                diags.push_back(i);
                missed++;
            }
        }
        int x=Ereco_copy.GetNrows()-missed; 
        int y=Ereco_copy.GetNcols()-missed;
        TMatrixD Ereco_copy_cut(x,y);
        unsigned int v=0;
        for (int i=0;i<Ereco_copy.GetNrows();i++){
            if(v<diags.size() && diags[v]==i){
                v++;
            }
            else{
                for (int j=0; j<Ereco_copy_cut.GetNcols();j++){                 
                    Ereco_copy_cut(i-v,j)=Ereco_copy(i,j+v);
                    }
                }
        }
    return Ereco_copy_cut;
}

TMatrixD
RooUnfold::Ereco(ErrorTreatment witherror)
{
    /*Returns covariance matrices for error calculation of type witherror
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix given by GetErrMat()
    */
    TMatrixD Ereco_m;
    switch(witherror){
      case kNoError:
        {
          TH1* HR=Hreco(kNoError);
          Ereco_m.ResizeTo(HR->GetNbinsX(),HR->GetNbinsX());
          for (int i=0; i<HR->GetNbinsX(); i++){
            Ereco_m(i,i)=HR->GetBinError(i);
          }
        }
        break;
      case kErrors:
        if (!_haveErrors) GetErrors();
        Ereco_m.ResizeTo(_errors.GetNrows(),_errors.GetNrows());
        for (int i=0; i<_errors.GetNrows();i++){
            Ereco_m(i,i)=_errors(i);
        }
        break;
      case kCovariance:
        if (!_haveCov) GetCov();
        Ereco_m.ResizeTo(_cov.GetNrows(),_cov.GetNcols());
        Ereco_m=_cov;
        break;
      case kCovToy:
        if (!_have_err_mat) GetErrMat();
        Ereco_m.ResizeTo(_err_mat.GetNrows(),_err_mat.GetNcols());
        Ereco_m=_err_mat;
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<witherror<<endl;
    }
    return Ereco_m;
}

TVectorD
RooUnfold::ErecoV(ErrorTreatment witherror)
{
    /*Returns vectors of the diagonals of the covariance matrices for error calculation of type witherror
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix given by GetErrMat()
    */
    TVectorD Ereco_m;
    switch(witherror){
      case kNoError:
        {
          TH1* HR=Hreco(kNoError);
          Ereco_m.ResizeTo(HR->GetNbinsX());
          for (int i=0; i<HR->GetNbinsX(); i++){
            Ereco_m(i)=HR->GetBinError(i);
          }
        }
        break;
      case kErrors:
        if (!_haveErrors) GetErrors();
        Ereco_m.ResizeTo(_errors.GetNrows());
        Ereco_m=_errors;
        break;
      case kCovariance:
        if (!_haveCov) GetCov();
        Ereco_m.ResizeTo(_cov.GetNrows());
        for (int i=0; i<_cov.GetNrows(); i++){
            Ereco_m(i)=_cov(i,i);
        }
        break;
      case kCovToy:
        if (!_have_err_mat) GetErrMat();
        Ereco_m.ResizeTo(_err_mat.GetNrows(),_err_mat.GetNcols());
        for (int i=0; i<_err_mat.GetNrows(); i++){
            Ereco_m(i)=_err_mat(i,i);
        }
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<witherror<<endl;
    }
    return Ereco_m;
}

Bool_t
RooUnfold::HaveErrors(ErrorTreatment witherror)
{
    if (!_unfolded) return false;
    switch(witherror){
      case kNoError:
        return true;
      case kErrors:
        return _haveErrors;
      case kCovariance:
        return _haveCov;
      case kCovToy:
        return _have_err_mat;
      default:
        cerr<<"Error, unrecognised error method= "<<witherror<<endl;
    }
    return false;
}
