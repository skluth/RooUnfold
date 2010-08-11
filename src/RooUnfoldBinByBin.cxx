//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBinByBin.cxx,v 1.7 2010-08-11 19:27:37 adye Exp $
//
// Description:
//      Unfolding bin-by-bin. Just an interface to RooUnfoldBayesImpl.
//      Note that we use 1D distributions in RooUnfoldBayesImpl, even if we are
//      unfolding multi-dimensional distributions: RooUnfold already converted
//      to 1D. Except for smoothing (which RooUnfoldBayesImpl doesn't implement)
//      this is just a matter of bookkeeping.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldBinByBin.h"

#include <vector>

#include "RooUnfoldBayes.h"
#include "RooUnfoldBayesImpl.h"

class TH1;
class RooUnfoldResponse;
using std::vector;

ClassImp (RooUnfoldBinByBin);

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldBinByBin& rhs)
  : RooUnfoldBayes (rhs) {}

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit,
                                      const char* name, const char* title)
  : RooUnfoldBayes(res, meas, 1, smoothit, name, title) {}

// Override RooUnfoldBayes routines.
Int_t RooUnfoldBinByBin::unfold (vector<Double_t>& causes) { return _bayes->unfoldBinByBin (causes, _smoothit); }
Int_t RooUnfoldBinByBin::getCovariance() const             { return _bayes->getCovarianceBinByBin();            }

void
RooUnfoldBinByBin::GetSettings(){
	_minparm=0;
	_maxparm=0;
	_stepsizeparm=0;
	_defaultparm=0;
}
