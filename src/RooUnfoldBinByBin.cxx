//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Unfolding class using the bin by bin method of conversion factors. 
//
// Authors: Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p> Uses the correction factor method to unfold the distribution by looking at each bin individually.</p>
<p> This method cannot account for bin migration and as such cannot unfold reliably if a bias/smearing effects are applied.</p>
<p>Can only handle 1 dimensional distributions
<p>True and measured distributions must have the same binning
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldBinByBin.h"

#include "TH1.h"
#include "TH2.h"

#include "RooUnfoldResponse.h"

ClassImp (RooUnfoldBinByBin);

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldBinByBin& rhs)
  : RooUnfold (rhs)
{
  // Copy constructor.
  GetSettings();  
}

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas, 
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  GetSettings();
}

RooUnfoldBinByBin::~RooUnfoldBinByBin()
{
}

RooUnfoldBinByBin*
RooUnfoldBinByBin::Clone (const char* newname) const
{
    //Clones object
  RooUnfoldBinByBin* unfold= new RooUnfoldBinByBin(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}



TVectorD* RooUnfoldBinByBin::Impl()
{
  return &_factors;
}

void
RooUnfoldBinByBin::Unfold()
{
    const TVectorD& vmeas=  Vmeasured();
    const TVectorD& vtrain= _res->Vmeasured();
    const TVectorD& vtruth= _res->Vtruth();
    _rec.ResizeTo(_nt);
    _factors.ResizeTo(_nt);
    Int_t nb= _nm < _nt ? _nm : _nt;
    for (int i=0; i<nb; i++) {
      if (vtrain[i]==0.0) continue;
      Double_t c= vtruth[i]/vtrain[i];
      _factors[i]= c;
      _rec(i)= c*vmeas[i];
    }
    _unfolded= true;
}

void
RooUnfoldBinByBin::GetCov()
{
    const TMatrixD& covmeas= GetMeasuredCov();
    _cov.ResizeTo(_nt,_nt);
    Int_t nb= _nm < _nt ? _nm : _nt;
    for (int i=0; i<nb; i++)
      for (int j=0; j<nb; j++)
        _cov(i,j)= _factors[i]*_factors[j]*covmeas(i,j);
    _haveCov= true;
}

void
RooUnfoldBinByBin::GetSettings(){
    _minparm=0;
    _maxparm=0;
    _stepsizeparm=0;
    _defaultparm=0;
}
