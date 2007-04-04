//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldSvd.cxx,v 1.1.1.1 2007-04-04 21:27:25 adye Exp $
//
// Description:
//      Unfold
//
// Author List:
//      Tim Adye <T.J.Adye@rl.ac.uk>
//
// Copyright Information:
//      Copyleft (C) 2006 Rutherford Appleton Laboratory
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

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldSvd);

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldSvd& rhs)
  : RooUnfold (rhs.GetName(), rhs.GetTitle())
{
  Setup ();
  Setup (rhs);
}

RooUnfoldSvd::RooUnfoldSvd (const RooUnfoldResponse* res, const TH1* meas, Int_t tau, Int_t ntoys,
                            const char* name, const char* title)
  : RooUnfold (name, title)
{
  Setup ();
  Setup (res, meas, tau, ntoys);
}

RooUnfoldSvd&
RooUnfoldSvd::Clear()
{
  delete _svd;
  delete _meas1d;
  delete _train1d;
  delete _truth1d;
  RooUnfold::Clear();
  return *this;
}

RooUnfoldSvd&
RooUnfoldSvd::Setup()
{
  _svd= 0;
  _tau= _ntoys= 0;
  _meas1d= _train1d= _truth1d= 0;
  RooUnfold::Setup();
  return *this;
}

RooUnfoldSvd&
RooUnfoldSvd::Setup (const RooUnfoldSvd& rhs)
{
  return Setup (rhs.response(), rhs.Hmeasured(), rhs._tau, rhs._ntoys);
}

RooUnfoldSvd&
RooUnfoldSvd::Setup (const RooUnfoldResponse* res, const TH1* meas, Int_t tau, Int_t ntoys)
{
  RooUnfold::Setup (res, meas);
  _tau= tau;
  _ntoys= ntoys;

  _svd= new TUnfHisto (_nt, _nm);

  TMatrixD covMeas(_nm,_nm);
  for (Int_t i= 0; i<_nm; i++) {
    for (Int_t j= 0; j<_nm; j++) {
      covMeas(i,j)= 0.0;
    }
    Double_t err= _meas->GetBinError(i+1);
    covMeas(i,i)= err*err;
  }

  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _meas1d=  RooUnfoldResponse::H2H1D (_meas, _nm);
  _train1d= _res->Hmeasured1D();
  _truth1d= _res->Htruth1D();
  TH1::AddDirectory (oldstat);

  _svd->init (_meas1d, _train1d, _truth1d, _res->Hresponse(), false);

  _rec.ResizeTo (_nm);
  _rec= _svd->Unfold (_tau);
  //Get the covariance matrix for statistical uncertainties on measured spectrum
  _cov.ResizeTo (_nm, _nm);
  _cov= _svd->GetCov (covMeas, _meas1d, _ntoys, _tau);
  //Get the covariance matrix for statistical uncertainties on signal MC
  TMatrixD ucovTrain= _svd->GetMatStatCov (_ntoys, _tau);

  Double_t sf= (_truth1d->Integral() / _train1d->Integral()) * _meas1d->Integral();
  Double_t sf2= sf*sf;
  for (Int_t i= 0; i<_nt; i++) {
    _rec[i] *= sf;
    for (Int_t j= 0; j<_nt; j++) {
      _cov(i,j) += ucovTrain(i,j);
      _cov(i,j) *= sf2;
    }
  }
  return *this;
}
