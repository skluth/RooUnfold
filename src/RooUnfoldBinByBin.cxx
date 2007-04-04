//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldBinByBin.cxx,v 1.1.1.1 2007-04-04 21:27:25 adye Exp $
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

#include "RooUnfoldBinByBin.h"

#include <iostream>
#include <vector>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayesImpl.h"
#include "Array2D.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldBinByBin);

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit,
                                      const char* name, const char* title)
  : RooUnfoldBayes (name, title)
{
  Setup (res, meas, smoothit);
}

RooUnfoldBinByBin&
RooUnfoldBinByBin::Setup (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit)
{
  RooUnfold::Setup (res, meas);
  _smoothit= smoothit;

  _bayes= new RooUnfoldBayesImpl (GetName(), GetTitle());

  _bayes->build (1, vector<Int_t>(1,_nt), vector<Double_t>(1,0.0), vector<Double_t>(1,1.0),
                 1, vector<Int_t>(1,_nm), vector<Double_t>(1,0.0), vector<Double_t>(1,1.0));

  if (verbose() >= 1) Print();

  vector<Double_t> vtruth(_nt), vtrain(_nm), vmeasured(_nm);
  Array2D aresp(_nt,_nm);
  _bayes->setupTrain  (H2VD (_res->Htruth(),    vtruth),
                       H2VD (_res->Hmeasured(), vtrain),
                       H2AD (_res->Hresponse(), aresp));

  _bayes->trainBinByBin (_smoothit);

  if (verbose() >= 1) Print();

  _bayes->setupUnfold (H2VD (_meas,             vmeasured));

  if (verbose() >= 1) Print();

  vector<Double_t> causes;
  _bayes->unfoldBinByBin (causes);

  if (verbose() >= 1) Print();

  _rec.ResizeTo (_nm);
  VD2V (causes, _rec);
  _cov.ResizeTo (_nm, _nm);
  AD2M (_bayes->covariance(), _cov);

  return *this;
}
