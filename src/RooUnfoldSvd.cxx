//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldSvd.cxx,v 1.9 2010-05-25 17:34:03 adye Exp $
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
  delete _reshist;
}

void
RooUnfoldSvd::Init()
{
  _svd= 0;
  _meas1d= _train1d= _truth1d= 0;
  _reshist= 0;
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
  Int_t overflow= (_overflow ? 2 : 0), nt= _nt+overflow, nm= _nm+overflow;

  _svd= new TUnfHisto (nt, nm);

  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _meas1d=  (_meas->GetDimension() == 1) ? CopyOverflow (_meas) : RooUnfoldResponse::H2H1D (_meas, _nm);
  // FIXME: use of CopyOverflow leaks copy and can do unneccessary copy
  _train1d= CopyOverflow (_res->Hmeasured1D());
  _truth1d= CopyOverflow (_res->Htruth1D());
  _reshist= CopyOverflow2D (_res->Hresponse());
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
  if (_kterm > nb+overflow) {
    cerr << "RooUnfoldSvd invalid kterm=" << _kterm << " with " << nb << " bins" << endl;
    _fail= true;
    return;
  }

  if (_verbose>=1) cout << "SVD init " << _reshist->GetNbinsX() << " x " << _reshist->GetNbinsY() << " bins" << endl;
  _svd->init (_meas1d, _train1d, _truth1d, _reshist, false);

  _rec.ResizeTo (nt);
  if (_verbose>=1) cout << "SVD unfold kterm=" << _kterm << endl;
  _rec= _svd->Unfold (_kterm);
  Double_t sf= (_truth1d->Integral() / _train1d->Integral()) * _meas1d->Integral();
  for (Int_t i= 0; i<nt; i++) {
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

  Int_t overflow= (_overflow ? 2 : 0), nt= _nt+overflow, nm= _nm+overflow;

  TMatrixD covMeas(nm,nm);
  for (Int_t i= 0; i<nm; i++) {
    for (Int_t j= 0; j<nm; j++) {
      covMeas(i,j)= 0.0;
    }
    Double_t err= _meas->GetBinError(i+1);
    covMeas(i,i)= err*err;
  }

  //Get the covariance matrix for statistical uncertainties on measured spectrum
  _cov.ResizeTo (nt, nt);
  _cov= _svd->GetCov (covMeas, _meas1d, _ntoys, _kterm);
  //Get the covariance matrix for statistical uncertainties on signal MC
  TMatrixD ucovTrain= _svd->GetMatStatCov (_ntoys, _kterm);

  Double_t sf= (_truth1d->Integral() / _train1d->Integral()) * _meas1d->Integral();
  Double_t sf2= sf*sf;
  for (Int_t i= 0; i<nt; i++) {
    for (Int_t j= 0; j<nt; j++) {
      _cov(i,j) += ucovTrain(i,j);
      _cov(i,j) *= sf2;
    }
  }
  _haveCov= true;
}

TH1D*
RooUnfoldSvd::CopyOverflow (const TH1* h) const
{
  if (!_overflow) return (TH1D*) h->Clone();
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

TH2D*
RooUnfoldSvd::CopyOverflow2D (const TH1* h) const
{
  if (!_overflow) return (TH2D*)h->Clone();
  Int_t nx= h->GetNbinsX(), ny= h->GetNbinsX();
  Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax(), xb= (xhi-xlo)/nx;
  Double_t ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax(), yb= (yhi-ylo)/ny;
  nx += 2; ny += 2;
  TH2D* hx= new TH2D (h->GetName(), h->GetTitle(), nx, xlo-xb, xhi+xb, ny, ylo-yb, yhi+yb);
  for (Int_t i= 0; i < nx; i++) {
    for (Int_t j= 0; j < ny; j++) {
      hx->SetBinContent (i+1, j+1, h->GetBinContent (i, j));
      hx->SetBinError   (i+1, j+1, h->GetBinError   (i, j));
    }
  }
  return hx;
}
