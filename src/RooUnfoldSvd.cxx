//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldSvd.cxx,v 1.8 2010-05-20 22:50:59 adye Exp $
//
// Description:
//      SVD unfolding. Just an interface to RooUnfHistoSvd.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldSvd.h"

#include <iostream>

#include "TNamed.h"
#include "TH1.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "RooUnfoldResponse.h"
#include "RooUnfHistoSvd.h"

using std::cerr;
using std::endl;

ClassImp (RooUnfoldSvd);

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldSvd& rhs)
  : RooUnfold (rhs)
{
  Init();
  CopyData (rhs);
}

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t kterm, Int_t ntoys,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kterm(kterm), _ntoys(ntoys)
{
  Init();
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
}

void
RooUnfoldSvd::Init()
{
  _svd= 0;
  _meas1d= _train1d= _truth1d= 0;
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
  _ntoys= rhs._ntoys;
}

TObject*
RooUnfoldSvd::Impl()
{
  return _svd;
}

void
RooUnfoldSvd::Unfold()
{
  if (_fail) return;
  if (_res->GetDimensionTruth() != 1 || _res->GetDimensionMeasured() != 1) {
    cerr << "RooUnfoldSvd may not work very well for multi-dimensional distributions" << endl;
  }

  _svd= new TUnfHisto (_nt, _nm);

  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _meas1d=  RooUnfoldResponse::H2H1D (_meas, _nm);
  _train1d= _res->Hmeasured1D();
  _truth1d= _res->Htruth1D();
  TH1::AddDirectory (oldstat);

  if (_nt != _nm) {
    cerr << "RooUnfoldSvd requires the same number of bins in the truth and measured distributions" << endl;
    _fail= true;
    return;
  }
  if (_kterm < 0) {
    cerr << "RooUnfoldSvd invalid kterm: " << _kterm << endl;
    _fail= true;
    return;
  }
  Int_t nb= _nm < _nt ? _nm : _nt;
  if (_kterm > nb) {
    cerr << "RooUnfoldSvd invalid kterm=" << _kterm << " with " << nb << " bins" << endl;
    _fail= true;
    return;
  }

  _svd->init (_meas1d, _train1d, _truth1d, _res->Hresponse(), false);

  _rec.ResizeTo (_nt);
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
  if (!_unfolded) Unfold();
  if (_fail) return;

  TMatrixD covMeas(_nm,_nm);
  for (Int_t i= 0; i<_nm; i++) {
    for (Int_t j= 0; j<_nm; j++) {
      covMeas(i,j)= 0.0;
    }
    Double_t err= _meas->GetBinError(i+1);
    covMeas(i,i)= err*err;
  }

  //Get the covariance matrix for statistical uncertainties on measured spectrum
  _cov.ResizeTo (_nt, _nt);
  _cov= _svd->GetCov (covMeas, _meas1d, _ntoys, _kterm);
  //Get the covariance matrix for statistical uncertainties on signal MC
  TMatrixD ucovTrain= _svd->GetMatStatCov (_ntoys, _kterm);

  Double_t sf= (_truth1d->Integral() / _train1d->Integral()) * _meas1d->Integral();
  Double_t sf2= sf*sf;
  for (Int_t i= 0; i<_nt; i++) {
    for (Int_t j= 0; j<_nt; j++) {
      _cov(i,j) += ucovTrain(i,j);
      _cov(i,j) *= sf2;
    }
  }
  _haveCov= true;
}
