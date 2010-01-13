//==============================================================================
// File and Version Information:
//      $Id: RooUnfold.cxx,v 1.6 2010-01-13 00:18:21 adye Exp $
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
#include <iomanip>
#include <sstream>
#include <cmath>

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
using std::setw;
using std::setprecision;
using std::sqrt;
using std::fabs;

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

void
RooUnfold::Assign (const RooUnfold& rhs)
{
  if (this == &rhs) return;
  Clear();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  Setup (rhs);
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
RooUnfold::PrintTable (std::ostream& o, const TH1* hTrue, Bool_t withError) const
{
  const TH1* hMeas=      Hmeasured();
  const TH1* hReco=      Hreco (withError);
  const TH1* hTrainTrue= response()->Htruth();
  const TH1* hTrain=     response()->Hmeasured();

  std::ostringstream fmt;
  fmt.copyfmt (o);   // save original ostream format
  o << "==================================================================" << endl
    << setw(3) << ""      << setw(9) << "Train" << setw(9) << "Train"    << setw(9) << "Test"  << setw(9) << "Test"  << setw(9) << "Unfolded" << setw(9) << "Diff" << setw(9) << "Pull" << endl
    << setw(3) << "Bin"   << setw(9) << "Truth" << setw(9) << "Measured" << setw(9) << "Truth" << setw(9) << "Input" << setw(9) << "Output"   << endl
    << "==================================================================" << endl;

  Double_t chi2= 0.0;
  Int_t ndf= 0;
  Int_t maxbin= _nt < _nm ? _nm : _nt;
  for (Int_t i = 0 ; i <= maxbin+1; i++) {

    if (i <= 0) o << " <1";
    else if (i>maxbin) o << ">" << setw(2) << maxbin;
    else               o        << setw(3) << i;
    o << std::fixed << setprecision(0);
    if (i<=_nt+1)
      o << ' ' << setw(8) << hTrainTrue->GetBinContent(i);
    else
      o << setw(9) << " ";
    if (i<=_nm+1)
      o << ' ' << setw(8) << hTrain->GetBinContent(i);
    else
      o << setw(9) << " ";
    if (hTrue && i<=_nt+1)
      o << ' ' << setw(8) << hTrue->GetBinContent(i);
    else
      o << setw(9) << " ";
    if (i<=_nm+1)
      o << ' ' << setw(8) << hMeas->GetBinContent(i);
    else
      o << setw(9) << " ";
    o << setprecision(1);
    if (i<=_nt+1) {
      o << ' ' << setw(8) << hReco->GetBinContent(i);
      if (hTrue) {
        Double_t ydiff    = hReco->GetBinContent(i) - hTrue->GetBinContent(i);
        Double_t ydiffErr = hReco->GetBinError(i);
        o << ' ' << setw(8) << ydiff;
        if (ydiffErr!=0) {
          Double_t ypull = ydiff/ydiffErr;
          chi2 += ypull*ypull;
          ndf++;
          o << setprecision(3) << ' ' << setw(8) << ypull;
        }
      }
    }
    o << endl;
  }

  o << "==================================================================" << endl
    << setw(3) << "" << std::fixed << setprecision(0)
    << ' ' << setw(8) << hTrainTrue->Integral()
    << ' ' << setw(8) << hTrain->Integral();
  if (hTrue)
    o << ' ' << setw(8) << hTrue->Integral();
  else
    o << setw(9) << " ";
  o << ' ' << setw(8) << hMeas->Integral() << setprecision(1)
    << ' ' << setw(8) << hReco->Integral()
    << endl
    << "==================================================================" << endl;
  o.copyfmt (fmt);  // restore original ostream format
  o << "Chi^2/NDF=" << chi2 << "/" << ndf << endl;
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
  TH1* reco= (TH1*) _res->Htruth()->Clone(GetName());
  reco->Reset();
  reco->SetTitle (GetTitle());
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
