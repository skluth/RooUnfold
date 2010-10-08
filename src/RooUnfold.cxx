//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
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
#ifdef HAVE_DAGOSTINI
#include "RooUnfoldDagostini.h"
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;
using std::sqrt;
using std::fabs;

ClassImp (RooUnfold);

RooUnfold::RooUnfold (const RooUnfoldResponse* res, const TH1* meas, const char* name, const char* title)
  : TNamed (name, title)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  // Should not normally be used directly - instead, create an instance of one of RooUnfold's subclasses,
  // or use the New() static constructor.
  Init();
  Setup (res, meas);
}

RooUnfold* RooUnfold::New (Algorithm alg, const RooUnfoldResponse* res, const TH1* meas,Double_t regparm,
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
      unfold = new RooUnfoldInvert  (res,meas);
      break;
    case kDagostini:
#ifdef HAVE_DAGOSTINI
      unfold = new RooUnfoldDagostini (res,meas);
      break;
#else
      cerr << "RooUnfoldDagostini is not available" << endl;
      return 0;
#endif
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

RooUnfold* RooUnfold::Clone (const char* newname) const
{
    //Creates a copy of the unfold object
  RooUnfold* unfold= new RooUnfold(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void RooUnfold::Destroy()
{
  delete _vMes;
  delete _eMes;
}

RooUnfold::RooUnfold (const RooUnfold& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

void RooUnfold::Assign (const RooUnfold& rhs)
{
  if (this == &rhs) return;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  CopyData (rhs);
}

void RooUnfold::CopyData (const RooUnfold& rhs)
{
  Setup (rhs.response(), rhs.Hmeasured());
  SetVerbose (rhs.verbose());
  SetNToys   (rhs.NToys());
}

void RooUnfold::Reset()
{
  Destroy();
  Init();
}

void RooUnfold::Init()
{
  _res= 0;
  _vMes= _eMes= 0;
  _meas= 0;
  _nm= _nt= 0;
  _verbose= 1;
  _overflow= 0;
  _unfolded= _haveCov= _fail= _have_err_mat= _haveErrors= false;
  _NToys=50;
  GetSettings();
}

RooUnfold& RooUnfold::Setup (const RooUnfoldResponse* res, const TH1* meas)
{
  Reset();
  SetResponse (res);
  SetMeasured (meas);
  return *this;
}

void RooUnfold::SetMeasured (const TH1* meas)
{
  _meas= meas;
  delete _vMes; _vMes= 0;
  delete _eMes; _eMes= 0;
}

void RooUnfold::SetResponse (const RooUnfoldResponse* res)
{
  _res= res;
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  if (_overflow) {
    _nm += 2;
    _nt += 2;
  }
  SetNameTitleDefault();
}

void RooUnfold::Unfold()
{
    // Dummy unfolding - just copies input
  cout << "********************** " << ClassName() << ": dummy unfolding - just copy input **********************" << endl;
  _rec.ResizeTo (_nt);
  Int_t nb= _nm < _nt ? _nm : _nt;
  for (Int_t i= 0; i < nb; i++) {
    _rec(i)= Vmeasured()(i);
  }
  _unfolded= true;
}

void RooUnfold::GetErrors()
{
    //Creates vector of diagonals of covariance matrices.
    //This may be overridden if it can be computed more quickly without the covariance matrix.
    if (!_haveCov) GetCov();
    if (!_haveCov) return;
    _variances.ResizeTo(_nt);
    for (Int_t i= 0; i < _nt; i++) {
      _variances(i)= _cov(i,i);
    }
    _haveErrors= true;
}

void RooUnfold::GetCov()
{
    //Dummy routine to get covariance matrix. It should be overridden by derived classes.
  _cov.ResizeTo (_nt, _nt);
  Int_t nb= _nm < _nt ? _nm : _nt;
  for (Int_t i= 0; i < nb; i++) {
    Double_t err= Emeasured()(i);
    _cov(i,i)= err*err;
  }
  _haveCov= true;
}

void RooUnfold::GetErrMat()
{
    //Get covariance matrix from the variation of the results in toy MC tests, using the Runtoy method.
    TVectorD bc_vec(_nt);
    TMatrixD bc_mat(_nt,_nt);
    _err_mat.ResizeTo(_nt,_nt);
    for (Int_t k=0; k<_NToys; k++){
        TH1* res=Runtoy();
        for (Int_t i=0; i<_nt;i++){
            Double_t res_bci=RooUnfoldResponse::GetBinContent (res, i, _overflow);
            bc_vec[i]+=res_bci;
            for (Int_t j=0; j<_nt; j++){
                bc_mat(i,j) += res_bci * RooUnfoldResponse::GetBinContent (res, j, _overflow);
            }
        }
    }
    for (Int_t i=0; i<_nt;i++){
        for (Int_t j=0; j<_nt; j++){
            _err_mat(i,j)=(bc_mat(i,j)-(bc_vec[i]*bc_vec[j])/_NToys)/_NToys;

        }
    }
    _have_err_mat=true;
}

Bool_t RooUnfold::UnfoldWithErrors (ErrorTreatment withError)
{
  if (!_unfolded) {
    if (_fail) return false;
    Unfold();
    if (!_unfolded) {
      _fail= true;
      return false;
    }
  }
  Bool_t ok;
  switch (withError) {
    case kErrors:
      if   (!_haveErrors)   GetErrors();
      ok= _haveErrors;
      break;
    case kCovariance:
      if   (!_haveCov)      GetCov();
      ok= _haveCov;
      break;
    case kCovToy:
      if   (!_have_err_mat) GetErrMat();
      ok= _have_err_mat;
      break;
    default:
      ok= true;
  }
  if (!ok) _fail= true;
  return ok;
}

Double_t RooUnfold::Chi2(const TH1* hTrue,ErrorTreatment DoChi2)
{
    /*Calculates Chi squared. Method depends on value of DoChi2
    0: sum of (residuals/error)squared
    1: use errors propagated through the unfolding
    2: use covariance matrix returned from unfolding
    3: use covariance matrix from the variation of the results in toy MC tests
    Returns warnings for small determinants of covariance matrices and if the condition is very large.
    If a matrix has to be inverted also removes rows/cols with all their elements equal to 0*/

    if (DoChi2==kCovariance||DoChi2==kCovToy){
        TMatrixD ereco(_nt,_nt);
        ereco=Ereco(DoChi2);
        if (!_unfolded) return -1;

        TVectorD res(_nt);
        for (Int_t i = 0 ; i < _nt; i++) {
            Int_t it= RooUnfoldResponse::GetBin (hTrue, i, _overflow);
            if (hTrue->GetBinContent(it)!=0.0 || hTrue->GetBinError(it)>0.0) {
              res(i) = _rec(i) - hTrue->GetBinContent(it);
            }
        }
        //cutting out 0 elements
        TMatrixD ereco_cut=CutZeros(ereco);
        if (ereco_cut.GetNrows()<=0) return 0.0;
        TMatrixD res_matrix_cut(ereco_cut.GetNrows(),1);
        for (Int_t i=0, v=0; i<_nt; i++){
            if (ereco(i,i) != 0.0) res_matrix_cut(v++,0)= res(i);
        }
        Double_t Ereco_det=ereco_cut.Determinant();
        TMatrixD res_matrix_t=res_matrix_cut;
        res_matrix_t.T();
        if (fabs(Ereco_det)<1e-5 && _verbose>=1){
            cerr << "Warning: Small Determinant of Covariance Matrix =" << Ereco_det << endl;
            cerr << "Chi^2 may be invalid due to small determinant" << endl;
        }
        TDecompSVD svd(ereco_cut);
        Double_t cond=svd.Condition();
        if (_verbose>=1){
            cout<<"For Covariance matrix condition= "<<cond<<" determinant= "<<Ereco_det<<endl;
        }
        svd.MultiSolve(res_matrix_cut);
        TMatrixD chisq_nw = res_matrix_t*res_matrix_cut;
        double cond_max=1e17;
        if (cond>=cond_max && _verbose>=1){
            cerr << "Warning, very large matrix condition= "<< cond<<" chi^2 may be inaccurate"<<endl;
        }
        return chisq_nw(0,0);
    }
    else{
        Double_t chi2=0;
        TVectorD ereco(_nt);
        ereco= ErecoV(DoChi2);
        if (!_unfolded) return -1;
        for (Int_t i = 0 ; i < _nt; i++) {
            Int_t it= RooUnfoldResponse::GetBin (hTrue, i, _overflow);
            if (ereco(i)>0.0 &&
                hTrue->GetBinContent(it)!=0.0 || hTrue->GetBinError(it)>0.0) {
                Double_t ypull = (_rec(i) - hTrue->GetBinContent(it)) / ereco(i);
                chi2 += ypull*ypull;
            }
        }

        return chi2;
    }
}


void RooUnfold::PrintTable (std::ostream& o, const TH1* hTrue, ErrorTreatment withError)
{
    //Prints data from truth, measured and reconstructed data for each bin.
  const TH1* hReco=      Hreco (withError);
  if (!_unfolded) return;
  const TH1* hMeas=      Hmeasured();
  const TH1* hTrainTrue= response()->Htruth();
  const TH1* hTrain=     response()->Hmeasured();
  std::ostringstream fmt;
  fmt.copyfmt (o);
  Int_t dim= hReco->GetDimension(), ntxb= hReco->GetNbinsX()+2, ntyb= hReco->GetNbinsY()+2;
  if (hMeas->GetDimension() != dim || hMeas->GetNbinsX()+2 != ntxb || hMeas->GetNbinsY()+2 != ntyb) dim= 1;
  Int_t iwid= (dim==3) ? 8 : (dim==2) ? 7 : 5;
  const char* xwid= (dim==3) ? "===" : (dim==2) ? "==" : "";
  o << "===============================================================================" << xwid << endl
    << setw(iwid) << ""      << setw(9) << "Train" << setw(9) << "Train"    << setw(9) << "Test"  << setw(9) << "Test"  << setw(9) << "Unfolded" << setw(10)<<"Error on"<<setw(9) << "Diff" << setw(9) << "Pull" << endl
    << setw(iwid) << "Bin"   << setw(9) << "Truth" << setw(9) << "Measured" << setw(9) << "Truth" << setw(9) << "Input" << setw(9) << "Output"   << setw(10)<<"Unfolding"<<endl
    << "===============================================================================" << xwid << endl;
  Double_t true_train_tot=0;
  Double_t meas_train_tot=0;
  Double_t true_test_tot=0;
  Double_t meas_test_tot=0;
  Double_t unf_tot=0;
  Double_t chi2= 0.0;
  Int_t ndf= 0, first= (_overflow ? 0 : 1);
  Int_t maxbin= _nt < _nm ? _nm : _nt;
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
    if (i<_nt){
      o << ' ' << setw(8) << hTrainTrue->GetBinContent(it);
    }
    else
      o << setw(9) << " ";
    if (i<_nm)
      o << ' ' << setw(8) << hTrain->GetBinContent(im);
    else
      o << setw(9) << " ";
    if (hTrue && i<_nt)
      o << ' ' << setw(8) << hTrue->GetBinContent(it);
    else
      o << setw(9) << " ";
    if (i<_nm)
      o << ' ' << setw(8) << hMeas->GetBinContent(im);
    else
      o << setw(9) << " ";
    o << setprecision(1);
    if (i<_nt) {
      o << ' ' << setw(8) << hReco->GetBinContent(it);
      if (hTrue &&
          ((hReco->GetBinContent(it)!=0.0 || (withError && hReco->GetBinError(it)>0.0)) &&
           (hTrue->GetBinContent(it)!=0.0 || (withError && hTrue->GetBinError(it)>0.0)))) {
        Double_t ydiff    = hReco->GetBinContent(it) - hTrue->GetBinContent(it);
        Double_t ydiffErr = hReco->GetBinError(it);
        o << " " << setw(9) << ydiffErr;
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

  o << "===============================================================================" << xwid << endl
    << setw(iwid) << "" << std::fixed << setprecision(0)
    << ' ' << setw(8) << true_train_tot
    << ' ' << setw(8) << meas_train_tot;
  if (hTrue)
    o << ' ' << setw(8) << true_test_tot;
  else
    o << setw(9) << " ";
  o << ' ' << setw(8) << meas_test_tot << setprecision(1)
    << ' ' << setw(8) << unf_tot
    << ' ' << setw(9) << sqrt(meas_test_tot)*(true_train_tot/meas_train_tot)
    << ' ' << setw(8) << unf_tot-true_test_tot;
    if(hMeas->Integral()>0){
    o<< ' ' << setw(8) <<(unf_tot-true_test_tot)/sqrt(meas_test_tot);
    }
    o<< endl
    << "===============================================================================" << xwid << endl;
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

void RooUnfold::SetNameTitleDefault()
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

TH1* RooUnfold::Hreco (ErrorTreatment withError)
{
    /*Creates reconstructed distribution. Error calculation varies by withError:
    0: No errors
    1: Errors from the square root of the diagonals of the covariance matrix given by the unfolding
    2: Errors from the square root of of the covariance matrix given by the unfolding
    3: Errors from the square root of the covariance matrix from the variation of the results in toy MC tests
    */
  TH1* reco= (TH1*) _res->Htruth()->Clone(GetName());
  reco->Reset();
  reco->SetTitle (GetTitle());
  if (!UnfoldWithErrors (withError)) withError= kNoError;
  if (!_unfolded) return reco;

  for (Int_t i= 0; i < _nt; i++) {
    Int_t j= RooUnfoldResponse::GetBin (reco, i, _overflow);
    reco->SetBinContent (j,             _rec(i));
    if        (withError==kErrors){
      reco->SetBinError (j, sqrt (fabs (_variances(i))));
    } else if (withError==kCovariance){
      reco->SetBinError (j, sqrt (fabs (_cov(i,i))));
    } else if (withError==kCovToy){
      reco->SetBinError (j, sqrt (fabs (_err_mat(i,i))));
    }
  }

  return reco;
}

void RooUnfold::GetSettings()
{
    //Gets maximum and minimum parameters and step size
    _minparm=0;
    _maxparm=0;
    _stepsizeparm=0;
    _defaultparm=0;
}

Double_t RooUnfold::GetMinParm() const
{
    //Get minimum regularisation parameter for unfolding method
    return _minparm;
}

Double_t RooUnfold::GetMaxParm() const
{
    //Get maximum regularisation parameter for unfolding method
    return _maxparm;
}

Double_t RooUnfold::GetStepSizeParm() const
{
    //Get suggested step size for unfolding distribution
    return _stepsizeparm;
}

Double_t RooUnfold::GetDefaultParm() const
{
    //Get suggested regularisation parameter.
    return _defaultparm;
}

TH1* RooUnfold::Runtoy(ErrorTreatment withError,double* chi2, const TH1* hTrue) const
{
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
    TH1* hReco= unfold_copy->Hreco ((withError==kCovToy) ? kNoError : withError);
    if (chi2 && !hTrue){
        cerr<<"Error: can't calculate chi^2 without a truth distribution"<<endl;
    }
    if (chi2 && hTrue) *chi2=unfold_copy->Chi2(hTrue,withError);
    delete hMeas;
    return hReco;
}


TH1* RooUnfold::Add_Random(const TH1* hMeas_AR)
{
    //Adds a random number to the measured distribution before the unfolding//
    TH1* hMeas2= dynamic_cast<TH1*>(hMeas_AR->Clone());
    hMeas2->Reset();
    Int_t nb= (hMeas_AR->GetNbinsX()+2) * (hMeas_AR->GetNbinsY()+2) * (hMeas_AR->GetNbinsZ()+2);
    for (Int_t i=0; i<nb ; i++){
            Double_t err=hMeas_AR->GetBinError(i);
            Double_t new_x = hMeas_AR->GetBinContent(i);
            if (err>0.0) new_x += gRandom->Gaus(0,err);
            hMeas2->SetBinContent(i,new_x);
            hMeas2->SetBinError(i,err);
    }
    return hMeas2;
}

void RooUnfold::Print(Option_t *opt) const
{
    TNamed::Print(opt);
    cout <<"regularisation parameter = "<<GetRegParm()<<", ntoys = "<<NToys()<<endl;
}

TMatrixD RooUnfold::CutZeros(const TMatrixD& ereco)
{
    //Removes row & column if all their elements are 0.
    vector<int> diags;
        int missed=0;
        for (int i=0; i<ereco.GetNrows(); i++){
            double coltot=0;
            for (int j=0;j<ereco.GetNrows();j++){
                coltot+=ereco(i,j);
            }
            if (coltot==0){
                diags.push_back(i);
                missed++;
            }
        }
        int x=ereco.GetNrows()-missed;
        int y=ereco.GetNcols()-missed;
        TMatrixD ereco_cut(x,y);
        unsigned int v=0;
        for (int i=0;i<ereco.GetNrows();i++){
            if(v<diags.size() && diags[v]==i){
                v++;
            }
            else{
                for (int j=0; j<ereco_cut.GetNcols();j++){
                    ereco_cut(i-v,j)=ereco(i,j+v);
                    }
                }
        }
    return ereco_cut;
}

TMatrixD RooUnfold::Ereco(ErrorTreatment withError)
{
    /*Returns covariance matrices for error calculation of type withError
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */
    TMatrixD Ereco_m(_nt,_nt);
    if (!UnfoldWithErrors (withError)) return Ereco_m;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          Ereco_m(i,i)=_rec(i);
        }
        break;
      case kErrors:
        for (int i=0; i<_nt;i++){
          Ereco_m(i,i)=_variances(i);
        }
        break;
      case kCovariance:
        Ereco_m=_cov;
        break;
      case kCovToy:
        Ereco_m=_err_mat;
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<withError<<endl;
    }
    return Ereco_m;
}

TVectorD RooUnfold::ErecoV(ErrorTreatment withError)
{
    /*Returns vector of unfolding errors computed according to the withError flag:
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */
    TVectorD Ereco_v(_nt);
    if (!UnfoldWithErrors (withError)) return Ereco_v;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_rec(i)));
        }
        break;
      case kErrors:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_variances(i)));
        }
        break;
      case kCovariance:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_cov(i,i)));
        }
        break;
      case kCovToy:
        for (int i=0; i<_nt; i++){
          Ereco_v(i)=sqrt (fabs (_err_mat(i,i)));
        }
        break;
      default:
        cerr<<"Error, unrecognised error method= "<<withError<<endl;
    }
    return Ereco_v;
}

TH1D* RooUnfold::HistNoOverflow (const TH1* h, Bool_t overflow)
{
  if (!overflow) {   // also for 2D+
    TH1D* hx= RooUnfoldResponse::H2H1D (h, h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ());
    if (!hx) return hx;
    // clear under/overflow bins for cloned TH1D
    hx->SetBinContent (0,                 0.0);
    hx->SetBinContent (hx->GetNbinsX()+1, 0.0);
    return hx;
  }
  Int_t nb= h->GetNbinsX();
  Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax(), xb= (xhi-xlo)/nb;
  nb += 2;
  TH1D* hx= new TH1D (h->GetName(), h->GetTitle(), nb, xlo-xb, xhi+xb);
  for (Int_t i= 0; i < nb; i++) {
    hx->SetBinContent (i+1, h->GetBinContent (i));
    hx->SetBinError   (i+1, h->GetBinError   (i));
  }
  return hx;
}
