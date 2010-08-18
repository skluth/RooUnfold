//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBinByBinOld.cxx,v 1.1 2010-08-18 12:58:04 fwx38934 Exp $
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

//____________________________________________________________
/* BEGIN_HTML
<p> Uses the correction factor method to unfold the distribution by looking at each bin individually.</p>
<p> This method cannot account for bin migration and as such cannot unfold reliably if a bias/smearing effects are applied.</p>
<p> The variables for unfolding are returned using the RooUnfoldBayes class, rather than being calculated independently</p>
END_HTML */

/////////////////////////////////////////////////////////////


#include "RooUnfoldBinByBinOld.h"

#include <vector>

#include "RooUnfoldBayes.h"
#include "RooUnfoldBayesImpl.h"

class TH1;
class RooUnfoldResponse;
using std::vector;

ClassImp (RooUnfoldBinByBinOld);

RooUnfoldBinByBinOld::RooUnfoldBinByBinOld (const RooUnfoldBinByBinOld& rhs)
  : RooUnfoldBayes (rhs) {Init();}

RooUnfoldBinByBinOld::RooUnfoldBinByBinOld (const RooUnfoldResponse* res, const TH1* meas, Bool_t smoothit,
                                      const char* name, const char* title)
  : RooUnfoldBayes(res, meas, 1, smoothit, name, title) {Init();}

// Override RooUnfoldBayes routines.
Int_t RooUnfoldBinByBinOld::unfold (vector<Double_t>& causes) { return _bayes->unfoldBinByBin (causes, _smoothit); }
Int_t RooUnfoldBinByBinOld::getCovariance() const             { return _bayes->getCovarianceBinByBin();            }
void
RooUnfoldBinByBinOld::GetSettings(){
	_minparm=0;
	_maxparm=0;
	_stepsizeparm=0;
	_defaultparm=0;
}
