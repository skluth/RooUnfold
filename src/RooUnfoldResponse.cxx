//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldResponse.cxx,v 1.2 2007-04-06 00:43:41 adye Exp $
//
// Description:
//      Response Matrix
//
// Author List:
//      Tim Adye <T.J.Adye@rl.ac.uk>
//
// Copyright Information:
//      Copyleft (C) 2006 Rutherford Appleton Laboratory
//
//==============================================================================

#include "RooUnfoldResponse.h"

#include <iostream>
#include <assert.h>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorD.h"
#include "TMatrixD.h"

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldResponse);

RooUnfoldResponse::RooUnfoldResponse (const RooUnfoldResponse& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  Setup ();
  Setup (rhs);
}

RooUnfoldResponse::RooUnfoldResponse (Int_t nb, Double_t xlo, Double_t xhi,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Setup ();
  Setup (nb, xlo, xhi);
}

RooUnfoldResponse::RooUnfoldResponse (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Setup ();
  Setup (nm, mlo, mhi, nt, tlo, thi);
}

RooUnfoldResponse::RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2D* response,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Setup ();
  Setup (measured, truth, response);
}

RooUnfoldResponse::RooUnfoldResponse (const TH1* measured, const TH1* truth,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Setup ();
  Setup (measured, truth);
}

RooUnfoldResponse&
RooUnfoldResponse::operator= (const RooUnfoldResponse& rhs)
{
  if (this == &rhs) return *this;
  Clear();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  return Setup (rhs);
}

RooUnfoldResponse&
RooUnfoldResponse::Clear()
{
  ClearCache();
  delete _mes;
  delete _tru;
  delete _res;
  return Setup();
}

RooUnfoldResponse&
RooUnfoldResponse::Setup()
{
  _tru= _mes= 0;
  _res= 0;
  _vMes= _eMes= _vTru= _eTru= 0;
  _mRes= _eRes= 0;
  _nm= _nt= _mdim= _tdim= 0;
  _cached= false;
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const RooUnfoldResponse& rhs)
{
  return Setup (rhs.Hmeasured(), rhs.Htruth(), rhs.Hresponse());
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi)
{
  Clear();
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= new TH1D (GetName(), GetTitle(), nm, mlo, mhi);
  _tru= new TH1D (GetName(), GetTitle(), nt, tlo, thi);
  _mdim= _tdim= 1;
  _nm= nm;
  _nt= nt;
  _res= new TH2D (GetName(), GetTitle(), nm, mlo, mhi, nt, tlo, thi);
  TH1::AddDirectory (oldstat);
  SetNameTitleDefault();
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const TH1* measured, const TH1* truth)
{
  Clear();
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= (TH1*) measured ->Clone();
  _tru= (TH1*) truth    ->Clone();
  _mes->Reset();
  _tru->Reset();
  _mdim= _mes->GetDimension();
  _tdim= _tru->GetDimension();
  _nm= _mes->GetNbinsX() * _mes->GetNbinsY() * _mes->GetNbinsZ();
  _nt= _tru->GetNbinsX() * _tru->GetNbinsY() * _tru->GetNbinsZ();
  _res= new TH2D (GetName(), GetTitle(),
                  _nm, _mdim==1 ? _mes->GetXaxis()->GetXmin() : 0.0, _mdim==1 ? _mes->GetXaxis()->GetXmax() : Double_t(_nm),
                  _nt, _tdim==1 ? _tru->GetXaxis()->GetXmin() : 0.0, _tdim==1 ? _tru->GetXaxis()->GetXmax() : Double_t(_nt));
  TH1::AddDirectory (oldstat);
  SetNameTitleDefault();
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const TH1* measured, const TH1* truth, const TH2D* response)
{
  Clear();
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= (TH1*)  measured->Clone();
  _tru= (TH1*)  truth   ->Clone();
  _res= (TH2D*) response->Clone();
  TH1::AddDirectory (oldstat);
  _mdim= _mes->GetDimension();
  _tdim= _tru->GetDimension();
  _nm= _mes->GetNbinsX() * _mes->GetNbinsY() * _mes->GetNbinsZ();
  _nt= _tru->GetNbinsX() * _tru->GetNbinsY() * _tru->GetNbinsZ();
  if (_nm != _res->GetNbinsX() || _nt != _res->GetNbinsY()) {
    cerr << "Warning: RooUnfoldResponse measured X truth is " << _nm << " X " << _nt
         << ", but matrix is " << _res->GetNbinsX()<< " X " << _res->GetNbinsY() << endl;
  }
  SetNameTitleDefault();
  return *this;
}

void
RooUnfoldResponse::ClearCache()
{
  delete _vMes; _vMes= 0;
  delete _eMes; _eMes= 0;
  delete _vTru; _vTru= 0;
  delete _eTru; _eTru= 0;
  delete _mRes; _mRes= 0;
  delete _eRes; _eRes= 0;
  _cached= false;
}

Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t xt)
{
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (_cached) ClearCache();
  ((TH1D*)_mes)->Fill (xr);
  ((TH1D*)_tru)->Fill (xt);
  return _res->Fill (xr, xt);
}

Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt)
{
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (_cached) ClearCache();
  ((TH2D*)_mes)->Fill (xr, yr);
  ((TH2D*)_tru)->Fill (xt, yt);
  return _res->Fill (Double_t(_mes->FindBin (xr, yr))-.5, Double_t(_tru->FindBin (xt, yt))-.5);
}

Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt)
{
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (_cached) ClearCache();
  ((TH3D*)_mes)->Fill (xr, yr, zr);
  ((TH3D*)_tru)->Fill (xt, yt, zt);
  return _res->Fill (Double_t(_mes->FindBin (xr, yr, zt))-.5, Double_t(_tru->FindBin (xt, yt, zt))-.5);
}

Int_t
RooUnfoldResponse::Miss (Double_t xt)
{
  assert (_tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (_cached) ClearCache();
  return ((TH1D*)_tru)->Fill (xt);
}

Int_t
RooUnfoldResponse::Miss (Double_t xt, Double_t yt)
{
  assert (_tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (_cached) ClearCache();
  return ((TH2D*)_tru)->Fill (xt, yt);
}

Int_t
RooUnfoldResponse::Miss (Double_t xt, Double_t yt, Double_t zt)
{
  assert (_tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (_cached) ClearCache();
  return ((TH3D*)_tru)->Fill (xt, yt, zt);
}

TH1D*
RooUnfoldResponse::H2H1D(const TH1* h, Int_t nb)
{
  if (h->GetDimension() == 1) return (TH1D*) h->Clone();
  TH1D* h1d= new TH1D(h->GetName(), h->GetTitle(), nb, 0.0, 1.0);
  for (size_t i= 1; i <= nb; i++) {
    h1d->SetBinContent (i, h->GetBinContent (i));
    h1d->SetBinError   (i, h->GetBinError   (i));
  }
  return h1d;
}

TVectorD*
RooUnfoldResponse::H2V  (const TH1* h, Int_t nb)
{
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (size_t i= 0; i < nb; i++) {
    (*v)(i)= h->GetBinContent (i+1);
  }
  return v;
}

TVectorD*
RooUnfoldResponse::H2VE (const TH1* h, Int_t nb)
{
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (size_t i= 0; i < nb; i++) {
    (*v)(i)= h->GetBinError (i+1);
  }
  return v;
}

TMatrixD*
RooUnfoldResponse::H2M  (const TH2D* h, Int_t nx, Int_t ny, const TH1* norm)
{
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (size_t j= 0; j < ny; j++) {
    Double_t nTrue= norm ? norm->GetBinContent (j+1) : 1.0;
    if (nTrue == 0.0) {
      for (size_t i= 0; i < nx; i++) {
        (*m)(i,j)= 0.0;
      }
    } else {
      for (size_t i= 0; i < nx; i++) {
        (*m)(i,j)= h->GetBinContent(i+1,j+1) / nTrue;
      }
    }
  }
  return m;
}

TMatrixD*
RooUnfoldResponse::H2ME (const TH2D* h, Int_t nx, Int_t ny, const TH1* norm)
{
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (size_t j= 0; j < ny; j++) {
    Double_t nTrue= norm ? norm->GetBinContent (j+1) : 1.0;
    if (nTrue == 0.0) {
      for (size_t i= 0; i < nx; i++) {
        (*m)(i,j)= 0.0;
      }
    } else {
      for (size_t i= 0; i < nx; i++) {
        // Assume Poisson nTrue, Multinomial P(mes|tru)
        (*m)(i,j)= h->GetBinError(i+1,j+1) / nTrue;
      }
    }
  }
  return m;
}

void
RooUnfoldResponse::SetNameTitleDefault()
{
  const char* s= GetName();
  if (s[0] == '\0') {
    if (_res) s= _res->GetName();
    if (s[0] == '\0') {
      if (_mes && _tru) {
        TString n= _mes->GetName();
        n.Append ("_");
        n.Append (_tru->GetName());
        SetName (n);
      }
    } else
      SetName (s);
  }
  s= GetTitle();
  if (s[0] == '\0') {
    if (_res) s= _res->GetTitle();
    if (s[0] == '\0') {
      if (_mes && _tru) {
        TString n= "Response ";
        n.Append (_tru->GetTitle());
        n.Append (" -> ");
        n.Append (_mes->GetTitle());
        SetTitle (n);
      }
    } else
      SetTitle (s);
  }
}

void
RooUnfoldResponse::Streamer (TBuffer &R__b)
{
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldResponse::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldResponse::Class()->WriteBuffer (R__b, this);
  }
}
