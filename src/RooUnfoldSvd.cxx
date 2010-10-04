//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldSvd.cxx,v 1.24 2010-09-10 23:58:02 adye Exp $
//
// Description:
//      SVD unfolding. Just an interface to RooUnfHistoSvd.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>Links to RooUnfHistoSvd class which unfolds using Singular Value Decomposition (SVD).</p>
<p>Regularisation parameter defines the level at which values are deemed to be due to statistical fluctuations and are cut out. (Default= number of bins/2)
<p>Returns errors as a full matrix of covariances
<p>Error processing is much the same as with the kCovToy setting with 1000 toys. This is quite slow but can be switched off.
<p>Can only handle 1 dimensional distributions
<p>True and measured distributions must have the same binning
<p>Can account for both smearing and biasing
<p>Returns near singular covariance matrices, again leading to very large chi squared values
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldSvd.h"

#include <iostream>
#include <iomanip>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "RooUnfoldResponse.h"
#include "RooUnfHistoSvd.h"

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldSvd);

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldSvd& rhs)
  : RooUnfold (rhs)
{
  Init();
  CopyData (rhs);
}

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kterm, Int_t ntoyssvd,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kterm(kterm ? kterm : meas->GetNbinsX()/2), _ntoyssvd(ntoyssvd)
{
  Init();
}

RooUnfoldSvd*
RooUnfoldSvd::Clone (const char* newname) const
{
  RooUnfoldSvd* unfold= new RooUnfoldSvd(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void
RooUnfoldSvd::Reset()
{
  Destroy();
  Init();
  RooUnfold::Reset();
}

void
RooUnfoldSvd::Destroy()
{
  delete _svd;
  delete _meas1d;
  delete _train1d;
  delete _truth1d;
  delete _reshist;
}

void
RooUnfoldSvd::Init()
{
  _svd= 0;
  _meas1d= _train1d= _truth1d= 0;
  _reshist= 0;
  _prop_errors=false;
  GetSettings();
}

void
RooUnfoldSvd::Assign (const RooUnfoldSvd& rhs)
{
  RooUnfold::Assign (rhs);
  CopyData (rhs);
}

void
RooUnfoldSvd::CopyData (const RooUnfoldSvd& rhs)
{
  _kterm= rhs._kterm;
  _ntoyssvd= rhs._ntoyssvd;
  _prop_errors= rhs._prop_errors;
}

TUnfHisto*
RooUnfoldSvd::Impl()
{
  return _svd;
}

void
RooUnfoldSvd::Unfold()
{
  if (_res->GetDimensionTruth() != 1 || _res->GetDimensionMeasured() != 1) {
    cerr << "RooUnfoldSvd may not work very well for multi-dimensional distributions" << endl;
  }

  _svd= new TUnfHisto (_nt, _nm);
  if (_prop_errors) _svd->UsePropErrors();

  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _meas1d=  HmeasuredNoOverflow1D();
  _train1d= HistNoOverflow (_res->Hmeasured1D(), _overflow);
  _truth1d= HistNoOverflow (_res->Htruth1D(),    _overflow);
  _reshist= _res->HresponseNoOverflow();
  TH1::AddDirectory (oldstat);
  if (_nt != _nm) {
    cerr << "RooUnfoldSvd requires the same number of bins in the truth and measured distributions" << endl;
    return;
  }
  if (_kterm < 0) {
    cerr << "RooUnfoldSvd invalid kterm: " << _kterm << endl;
    return;
  }
  Int_t nb= _nm < _nt ? _nm : _nt;
  if (_kterm > nb) {
    cerr << "RooUnfoldSvd invalid kterm=" << _kterm << " with " << nb << " bins" << endl;
    return;
  }

  if (_verbose>=1) cout << "SVD init " << _reshist->GetNbinsX() << " x " << _reshist->GetNbinsY() << " bins" << endl;
  _svd->init (_meas1d, _train1d, _truth1d, _reshist, false);

  _rec.ResizeTo (_nt);
  if (_verbose>=1) cout << "SVD unfold kterm=" << _kterm << endl;
  _rec= _svd->Unfold (_kterm);
  Double_t sf= (_truth1d->Integral() / _train1d->Integral()) * _meas1d->Integral();
  for (Int_t i= 0; i<_nt; i++) {
    _rec[i] *= sf;
  }
  _unfolded= true;
  _haveCov=  false;
}

void
RooUnfoldSvd::GetCov()
{
  TMatrixD covMeas(_nm,_nm);
  Int_t first= _overflow ? 0 : 1;
  for (Int_t i= 0; i<_nm; i++) {
    Double_t err= _meas->GetBinError(i+first);
    covMeas(i,i)= err*err;
  }

  //Get the covariance matrix for statistical uncertainties on measured spectrum
  _cov.ResizeTo (_nt, _nt);
  TMatrixD ucovTrain(_nt,_nt);
  if (!_prop_errors){
      _cov= _svd->GetCov(covMeas, _meas1d, _ntoyssvd, _kterm);
      //Get the covariance matrix for statistical uncertainties on signal MC
      ucovTrain=_svd->GetMatStatCov (_ntoyssvd, _kterm);
  }
  else{
    _cov= _svd->GetCovProp();
  }
  Double_t sf= (_truth1d->Integral() / _train1d->Integral()) * _meas1d->Integral();
  Double_t sf2= sf*sf;
  for (Int_t i= 0; i<_nt; i++) {
    for (Int_t j= 0; j<_nt; j++) {
        if (!_prop_errors){
            _cov(i,j) += ucovTrain(i,j);
        }
      _cov(i,j) *= sf2;
    }
  }
  _haveCov= true;
}

void
RooUnfoldSvd::GetSettings(){
    _minparm=0;
    _maxparm=_meas->GetNbinsX();
    _stepsizeparm=1;
    _defaultparm=_meas->GetNbinsX()/2;
}

void 
RooUnfoldSvd::UsePropErrors(Bool_t PE)
{
    //Allows covariance matrix calculation by propagation of errors
    //At present this produces very large error bars
    _prop_errors=PE;
}
