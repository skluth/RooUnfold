//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldSvd.cxx,v 1.18 2010-08-23 15:33:55 fwx38934 Exp $
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
  _cov= _svd->GetCov (covMeas, _meas1d, _ntoyssvd, _kterm);
  //Get the covariance matrix for statistical uncertainties on signal MC
  TMatrixD ucovTrain= _svd->GetMatStatCov (_ntoyssvd, _kterm);
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

void
RooUnfoldSvd::GetSettings(){
	_minparm=0;
	_maxparm=_meas->GetNbinsX();
	_stepsizeparm=1;
	_defaultparm=_meas->GetNbinsX()/2;
}
