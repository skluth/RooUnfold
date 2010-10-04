//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldResponse.cxx,v 1.19 2010-09-14 22:47:56 adye Exp $
//
// Description:
//      Response Matrix
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
 <p> Class to create response object as used in RooUnfold </p>
 <p> Contains measured and truth distributions as TH1s and the response matrix as a TH2D. Also contains methods for handling these data</p>
<p> Can handle 1,2 or 3 dimensional histograms and return vectors and matrices of their bin content and error (1 and 2D distributions respectively).
 Conversely can also convert these vectors and matrices into TH1Ds and TH2Ds. </p>
<p> Can also take a variety of parameters as inputs. This includes maximum and minimum values, distributions and vectors/matrices of values. </p>
<p> This class does the numerical modifications needed to allow unfolding techniques to work in the unfolding routines used in RooUnfold. </p>
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldResponse.h"

#include <iostream>
#include <assert.h>

#include "TClass.h"
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
  Init();
  Setup (rhs);
}

RooUnfoldResponse::RooUnfoldResponse (Int_t nb, Double_t xlo, Double_t xhi,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Init();
  Setup (nb, xlo, xhi);
}

RooUnfoldResponse::RooUnfoldResponse (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Init();
  Setup (nm, mlo, mhi, nt, tlo, thi);
}

RooUnfoldResponse::RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2D* response,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Init();
  Setup (measured, truth, response);
}

RooUnfoldResponse::RooUnfoldResponse (const TH1* measured, const TH1* truth,
                                      const char* name, const char* title)
  : TNamed (name, title)
{
  Init();
  Setup (measured, truth);
}

RooUnfoldResponse&
RooUnfoldResponse::operator= (const RooUnfoldResponse& rhs)
{
  if (this == &rhs) return *this;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  return Setup (rhs);
}

void
RooUnfoldResponse::Add (const RooUnfoldResponse& rhs)
{
  // Add another RooUnfoldResponse, accumulating contents
  if (_res == 0) {
    Setup (rhs);
    return;
  }
  assert (_mdim==rhs._mdim);
  assert (_tdim==rhs._tdim);
  assert (_mes != 0 && rhs._mes != 0);
  assert (_tru != 0 && rhs._tru != 0);
  assert (_res != 0 && rhs._res != 0);
  if (_cached) ClearCache();
  _mes->Add (rhs._mes);
  _tru->Add (rhs._tru);
  _res->Add (rhs._res);
}

RooUnfoldResponse&
RooUnfoldResponse::Reset()
{
    //Deletes all previous information. Sets parameters to 0
  ClearCache();
  delete _mes;
  delete _tru;
  delete _res;
  return Setup();
}

RooUnfoldResponse&
RooUnfoldResponse::Init()
{
    //Initialises all variables
  _overflow= 0;
  return Setup();
}

RooUnfoldResponse&
RooUnfoldResponse::Setup()
{
    //Initialise
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
    //Copy data from another RooUnfoldResponse
  _overflow= rhs._overflow;
  return Setup (rhs.Hmeasured(), rhs.Htruth(), rhs.Hresponse());
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi)
{
  Reset();
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= new TH1D ("measured", "Measured", nm, mlo, mhi);
  _tru= new TH1D ("truth",    "Truth",    nt, tlo, thi);
  _mdim= _tdim= 1;
  _nm= nm;
  _nt= nt;
  _res= new TH2D (GetName(), GetTitle(), nm, mlo, mhi, nt, tlo, thi);
  TH1::AddDirectory (oldstat);
  SetNameTitleDefault ("response", "Response");
  return *this;
}

void
RooUnfoldResponse::ReplaceAxis (TObject* hist, TAxis* dest, const TAxis* source)
{
    //Replaces an axis with that of a different histogram
  source->Copy (*dest);
  dest->SetParent (hist);
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const TH1* measured, const TH1* truth)
{
  Reset();
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= (TH1*) measured ->Clone();
  _tru= (TH1*) truth    ->Clone();
  _mes->Reset();
  _tru->Reset();
  _mdim= _mes->GetDimension();
  _tdim= _tru->GetDimension();
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  _nm= _mes->GetNbinsX() * _mes->GetNbinsY() * _mes->GetNbinsZ();
  _nt= _tru->GetNbinsX() * _tru->GetNbinsY() * _tru->GetNbinsZ();
  _res= new TH2D (GetName(), GetTitle(), _nm, 0.0, Double_t(_nm), _nt, 0.0, Double_t(_nt));
  if (_mdim==1) ReplaceAxis (_res, _res->GetXaxis(), _mes->GetXaxis());
  if (_tdim==1) ReplaceAxis (_res, _res->GetYaxis(), _tru->GetXaxis());
  TH1::AddDirectory (oldstat);
  SetNameTitleDefault();
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup (const TH1* measured, const TH1* truth, const TH2D* response)
{
    //Reads in measured and truth distribution as well as error matrix.
  Reset();
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  _mes= (TH1*)  measured->Clone();
  _tru= (TH1*)  truth   ->Clone();
  _res= (TH2D*) response->Clone();
  TH1::AddDirectory (oldstat);
  _mdim= _mes->GetDimension();
  _tdim= _tru->GetDimension();
  if (_overflow && (_mdim > 1 || _tdim > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
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
RooUnfoldResponse::Fill (Double_t xr, Double_t xt, Double_t w)
{
  // Fill 1D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (_cached) ClearCache();
  ((TH1D*)_mes)->Fill (xr, w);
  ((TH1D*)_tru)->Fill (xt, w);
  return _res->Fill (xr, xt, w);
}

Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt, Double_t w)
{
  // Fill 2D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (_cached) ClearCache();
  ((TH2D*)_mes)->Fill (xr, yr, w);
  ((TH2D*)_tru)->Fill (xt, yt, w);
  return _res->Fill (Double_t(FindBin (_mes, xr, yr))+.5, Double_t(FindBin (_tru, xt, yt))+.5, w);
}

Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill 3D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (_cached) ClearCache();
  ((TH3D*)_mes)->Fill (xr, yr, zr, w);
  ((TH3D*)_tru)->Fill (xt, yt, zt, w);
  return _res->Fill (Double_t(FindBin (_mes, xr, yr, zt))+.5, Double_t(FindBin (_tru, xt, yt, zt))+.5, w);
}

Int_t
RooUnfoldResponse::FindBin(const TH1* h, Double_t x, Double_t y)
{
  // Get vector index (0..nx*ny-1) for (x,y) coordinates
  Int_t nx=   h->GetNbinsX();
  Int_t ny=   h->GetNbinsY();
  Int_t binx= h->GetXaxis()->FindBin(x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny;
  Int_t biny= h->GetYaxis()->FindBin(y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny;
  return binx + nx*biny;
}

Int_t
RooUnfoldResponse::FindBin(const TH1* h, Double_t x, Double_t y, Double_t z)
{
  // Get vector index (0..nx*ny*nz-1) for (x,y,z) coordinates
  Int_t nx=   h->GetNbinsX();
  Int_t ny=   h->GetNbinsY();
  Int_t nz=   h->GetNbinsZ();
  Int_t binx= h->GetXaxis()->FindBin(x) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny*nz;
  Int_t biny= h->GetYaxis()->FindBin(y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny*nz;
  Int_t binz= h->GetZaxis()->FindBin(z) - 1;
  if (binz <  0)  return -1;
  if (binz >= nz) return nx*ny*nz;
  return binx + nx*(biny + ny*binz);
}

Int_t
RooUnfoldResponse::GetBinDim (const TH1* h, Int_t i)
{
  // Converts from vector index (0..nx*ny-1) or (0..nx*ny*nz-1) to multi-dimensional histogram
  // global bin number (0..(nx+2)*(ny+2)-1) or (0..(nx+2)*(ny+2)*(nz+2)-1), skipping under/overflow bins.
  Int_t ndim= h->GetDimension(), nx= h->GetNbinsX();
  if        (ndim == 2) {
//  cout << i << " -> " << "(" << i%nx+1 << "," << i/nx+1 << ")" << endl;
    return (i%nx+1) + (nx+2)*(i/nx+1);
  } else if (ndim == 3) {
    Int_t ny= h->GetNbinsY();
//  cout << i << " -> " << "(" << i%nx+1 << "," << (i/nx)%ny+1 << "," << i/(nx*ny)+1 << ")" << endl;
    return (i%nx+1) + (nx+2)*((i/nx)%ny+1 + (ny+2)*(i/(nx*ny)+1));
  }
  return i+1;   // not used: 1D handled by inline GetBin() (and handling UseOverflow), don't support >3D.
}

Int_t
RooUnfoldResponse::Miss1D (Double_t xt, Double_t w)
{
  // Fill missed event into 1D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_mdim==1 && _tdim==1);
  if (_cached) ClearCache();
  return ((TH1D*)_tru)->Fill (xt, w);
}

Int_t
RooUnfoldResponse::Miss2D (Double_t xt, Double_t yt, Double_t w)
{
  // Fill missed event into 2D Response Matrix (with weight)
  assert (_tru != 0);
  assert (_mdim==2 && _tdim==2);
  if (_cached) ClearCache();
  return ((TH2D*)_tru)->Fill (xt, yt, w);
}

Int_t
RooUnfoldResponse::Miss (Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  // Fill missed event into 3D Response Matrix
  assert (_tru != 0);
  assert (_mdim==3 && _tdim==3);
  if (_cached) ClearCache();
  return ((TH3D*)_tru)->Fill (xt, yt, zt, w);
}

TH1D*
RooUnfoldResponse::H2H1D(const TH1* h, Int_t nb)
{
  if (h->GetDimension() == 1) return (TH1D*) h->Clone();
  TH1D* h1d= new TH1D(h->GetName(), h->GetTitle(), nb, 0.0, 1.0);
  for (Int_t i= 0; i < nb; i++) {
    Int_t j= GetBin (h, i);  // don't bother with under/overflow bins (not supported for >1D)
    h1d->SetBinContent (i+1, h->GetBinContent (j));
    h1d->SetBinError   (i+1, h->GetBinError   (j));
  }
  return h1d;
}

TH2D*
RooUnfoldResponse::HresponseNoOverflow() const
{
  const TH2D* h= Hresponse();
  Int_t nx= h->GetNbinsX(), ny= h->GetNbinsY();
  if (!_overflow) {
    TH2D* hx= dynamic_cast<TH2D*>(h->Clone());
    if (!hx) return hx;
    // clear under/overflows
    for (Int_t i= 0; i <= nx+1; i++) {
      hx->SetBinContent (i, 0,    0.0);
      hx->SetBinContent (i, ny+1, 0.0);
    }
    for (Int_t i= 1; i <= ny;   i++) {
      hx->SetBinContent (0,    i, 0.0);
      hx->SetBinContent (nx+1, i, 0.0);
    }
    return hx;
  }
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

TVectorD*
RooUnfoldResponse::H2V  (const TH1* h, Int_t nb, Bool_t overflow)
{
    //Returns TVectorD of the bin contents of the input histogram
  if (overflow) nb += 2;
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (Int_t i= 0; i < nb; i++) {
    (*v)(i)= GetBinContent (h, i, overflow);
  }
  return v;
}

void
RooUnfoldResponse::V2H (const TVectorD& v, TH1* h, Int_t nb, Bool_t overflow)
{
    //Sets the bin content of the histogram as that element of the input vector
  h->Reset();  // in particular, ensure under/overflows are reset
  if (overflow) nb += 2;
  for (Int_t i= 0; i < nb; i++) {
    Int_t j= GetBin (h, i, overflow);
    h->SetBinContent (j, v(i));
  }
}

TVectorD*
RooUnfoldResponse::H2VE (const TH1* h, Int_t nb, Bool_t overflow)
{
    //Returns TVectorD of bin errors for input histogram
  if (overflow) nb += 2;
  TVectorD* v= new TVectorD (nb);
  if (!h) return v;
  for (Int_t i= 0; i < nb; i++) {
    (*v)(i)= GetBinError (h, i, overflow);
  }
  return v;
}

TMatrixD*
RooUnfoldResponse::H2M  (const TH2D* h, Int_t nx, Int_t ny, const TH1* norm, Bool_t overflow)
{
    //Returns Matrix of values of bins in a 2D input histogram
  Int_t first= overflow ? 0 : 1;
  if (overflow) {
    nx += 2;
    ny += 2;
  }
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (Int_t j= 0; j < ny; j++) {
    Double_t nTrue= norm ? GetBinContent (norm, j, overflow) : 1.0;
    if (nTrue == 0.0) {
      for (Int_t i= 0; i < nx; i++) {
        (*m)(i,j)= 0.0;
      }
    } else {
      for (Int_t i= 0; i < nx; i++) {
        (*m)(i,j)= h->GetBinContent(i+first,j+first) / nTrue;
      }
    }
  }
  return m;
}

TMatrixD*
RooUnfoldResponse::H2ME (const TH2D* h, Int_t nx, Int_t ny, const TH1* norm, Bool_t overflow)
{
    //Returns matrix of bin errors for a 2D histogram.
  Int_t first= overflow ? 0 : 1;
  if (overflow) {
    nx += 2;
    ny += 2;
  }
  TMatrixD* m= new TMatrixD (nx, ny);
  if (!h) return m;
  for (Int_t j= 0; j < ny; j++) {
    Double_t nTrue= norm ? GetBinContent (norm, j, overflow) : 1.0;
    if (nTrue == 0.0) {
      for (Int_t i= 0; i < nx; i++) {
        (*m)(i,j)= 0.0;
      }
    } else {
      for (Int_t i= 0; i < nx; i++) {
        // Assume Poisson nTrue, Multinomial P(mes|tru)
        (*m)(i,j)= h->GetBinError(i+first,j+first) / nTrue;
      }
    }
  }
  return m;
}

TH1*
RooUnfoldResponse::ApplyToTruth (const TH1* truth, const char* name) const
{
  // Apply the response matrix to the truth
  // Errors not set, since we assume original truth has no errors
  if (!Htruth()) return 0;  // Needed for checking binning if nothing else

  // If no truth histogram input, use training truth
  // If truth histogram input, make sure its binning is correct
  TVectorD* resultvect;
  if (truth) {
    if (truth->GetNbinsX() != _tru->GetNbinsX() ||
        truth->GetNbinsY() != _tru->GetNbinsY() ||
        truth->GetNbinsZ() != _tru->GetNbinsZ())
      cerr << "Warning: RooUnfoldResponse::ApplyToTruth truth histogram is a different size ("
           << (truth->GetNbinsX() * truth->GetNbinsY() * truth->GetNbinsZ()) << " bins) or shape from response matrix truth ("
           << ( _tru->GetNbinsX() *  _tru->GetNbinsY() *  _tru->GetNbinsZ()) << " bins)" << endl;
    resultvect= H2V (truth, GetNbinsTruth(), _overflow);
    if (!resultvect) return 0;
  } else {
    resultvect= new TVectorD (Vtruth());
  }

  (*resultvect) *= Mresponse();   // v= A*v

  // Turn results vector into properly binned histogram
  TH1* result= (TH1*) Hmeasured()->Clone (name);
  result->SetTitle (name);
  V2H (*resultvect, result, GetNbinsMeasured(), _overflow);
  delete resultvect;
  return result;
}

void
RooUnfoldResponse::SetNameTitleDefault (const char* defname, const char* deftitle)
{
  // Set object name and title
  const char* s= GetName();
  if (s[0] == '\0') {
    if (_res) s= _res->GetName();
    if (s[0] == '\0') {
      if (defname) SetName (defname);
      else if (_mes && _tru) {
        TString n= _mes->GetName();
        if (n.Length()) n.Append ("_");
        n.Append (_tru->GetName());
        if (!n.Length()) n= "response";
        SetName (n);
      }
    } else
      SetName (s);
  }
  s= GetTitle();
  if (s[0] == '\0') {
    if (_res) s= _res->GetTitle();
    if (s[0] == '\0') {
      if (deftitle) SetTitle (deftitle);
      else if (_mes && _tru) {
        TString n= _tru->GetTitle();
        if (n.Length()) n.Append (" #rightarrow ");
        n.Append (_mes->GetTitle());
        if (n.Length())
          n.Prepend ("Response ");
        else
          n= "Response";
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
