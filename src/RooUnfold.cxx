//==============================================================================
// File and Version Information:
//      $Id: RooUnfold.cxx,v 1.3 2009-05-22 17:10:20 adye Exp $
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

#include "RooUnfold.h"

#include <iostream>
#include <math.h>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorD.h"
#include "TMatrixD.h"

#include "RooUnfoldResponse.h"

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfold);

RooUnfold::RooUnfold (const RooUnfold& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  Setup ();
  Setup (rhs);
}

RooUnfold::RooUnfold (const RooUnfoldResponse* res, const TH1* meas, const char* name, const char* title)
  : TNamed (name, title)
{
  Setup ();
  Setup (res, meas);
}

RooUnfold&
RooUnfold::operator= (const RooUnfold& rhs)
{
  if (this == &rhs) return *this;
  Clear();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  return Setup (rhs);
}

RooUnfold&
RooUnfold::Setup (const RooUnfold& rhs)
{
  return Setup (rhs.response(), rhs.Hmeasured());
}

RooUnfold&
RooUnfold::Clear()
{
  return Setup();
}

RooUnfold&
RooUnfold::Setup()
{
  _res= 0;
  _meas= 0;
  _nm= _nt= 0;
  _verbose= 1;
  _haveCov= false;
  return *this;
}

RooUnfold&
RooUnfold::Setup (const RooUnfoldResponse* res, const TH1* meas)
{
  Clear();
  _res= res;
  _meas= meas;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  SetNameTitleDefault();
  return *this;
}

void
RooUnfold::SetNameTitleDefault()
{
  if (!_res) return;
  const char* s= GetName();
  if (s[0] == '\0') SetName (_res->GetName());
  s= GetTitle();
  if (s[0] == '\0') SetTitle (_res->GetTitle());
}

TH1*
RooUnfold::Hreco (Bool_t withError) const
{
  TH1* reco= (TH1*) _res->Htruth()->Clone();
  reco->Reset();
  reco->SetNameTitle (GetName(), GetTitle());
  if (withError && !_haveCov) GetCov();
  for (size_t i= 0; i < _nt; i++) {
    Int_t j= RooUnfoldResponse::GetBin(reco, i);
    reco->SetBinContent (j,             _rec(i));
    if (withError)
      reco->SetBinError (j, sqrt (fabs (_cov(i,i))));
  }
  return reco;
}

void RooUnfold::GetCov() const {
  _haveCov= true;
}
