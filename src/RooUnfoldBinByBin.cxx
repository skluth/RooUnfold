//==============================================================================
// File and Version Information:
//      $Id: RooUnfoldBinByBin.cxx,v 1.2 2009-05-22 17:10:20 adye Exp $
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

#include <vector>

#include "TNamed.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBayesImpl.h"

class TH1;
class RooUnfoldResponse;
using std::vector;

ClassImp (RooUnfoldBinByBin);

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit,
                                      const char* name, const char* title)
  : RooUnfoldBayes (res, meas, 0, smoothit, name, title) {}

RooUnfoldBinByBin&
RooUnfoldBinByBin::Setup (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit)
{
  RooUnfoldBayes::Setup (res, meas, 0, smoothit);
  return *this;
}

// Override RooUnfoldBayes routines.
Int_t RooUnfoldBinByBin::unfold (vector<Double_t>& causes) { return _bayes->unfoldBinByBin (causes);   }
Int_t RooUnfoldBinByBin::train()                           { return _bayes->trainBinByBin (_smoothit); }
Int_t RooUnfoldBinByBin::getCovariance() const             { return _bayes->getCovarianceBinByBin();   }
