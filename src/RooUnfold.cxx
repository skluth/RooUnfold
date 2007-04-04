//==============================================================================
// File and Version Information:
//      $Id: RooUnfold.cxx,v 1.1.1.1 2007-04-04 21:27:25 adye Exp $
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
  for (size_t i= 0; i < _nt; i++) {
    reco->SetBinContent (i+1,             _rec(i));
    if (withError)
      reco->SetBinError (i+1, sqrt (fabs (_cov(i,i))));
  }
  return reco;
}
