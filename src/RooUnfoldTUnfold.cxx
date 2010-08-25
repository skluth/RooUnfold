//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTUnfold.cxx,v 1.16 2010-08-25 10:23:49 adye Exp $
//
// Description:
//      Unfolding class using TUnfold from ROOT to do the actual unfolding.
//
// Authors: Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>Uses the unfolding method implemented in ROOT's <a href="http://root.cern.ch/root/html/TUnfold.html">TUnfold</a> class
<p>Only included in ROOT versions 5.22 and higher
<p>Only able to reconstruct 1 dimensional distributions
<p>Can account for bin migration and smearing
<p>Errors come as a full covariance matrix. 
<p>Will sometimes warn of "unlinked" bins. These are bins with 0 entries and do not effect the results of the unfolding
<p>Regularisation parameter can be either optimised internally by plotting log10(chi2 squared) against log10(tau). The 'kink' in this curve is deemed the optimum tau value. This value can also be set manually (FixTau)
<p>The latest version (TUnfold 15 in ROOT 2.27.04) will not handle plots with an additional underflow bin. As a result overflows must be turned off
if v15 of TUnfold is used. ROOT versions 5.26 or below use v13 and so should be safe to use overflows.</ul>
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldTUnfold.h"

#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TUnfold.h"
#include "TGraph.h"

#include "RooUnfoldResponse.h"

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldTUnfold);

RooUnfoldTUnfold::RooUnfoldTUnfold (const RooUnfoldTUnfold& rhs)
  : RooUnfold (rhs)
{
  Init();
  CopyData (rhs);
}

RooUnfoldTUnfold::RooUnfoldTUnfold (const RooUnfoldResponse* res, const TH1* meas, TUnfold::ERegMode reg,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title),_reg_method(reg)
{
  Init();
}

void
RooUnfoldTUnfold::Destroy()
{
    delete _unf;
}

RooUnfoldTUnfold*
RooUnfoldTUnfold::Clone (const char* newname) const
{
    //Clones object
  RooUnfoldTUnfold* unfold= new RooUnfoldTUnfold(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void
RooUnfoldTUnfold::Reset()
{
    //Resets all values
  Destroy();
  Init();
  RooUnfold::Reset();
}


void
RooUnfoldTUnfold::Assign (const RooUnfoldTUnfold& rhs)
{
  RooUnfold::Assign (rhs);
  CopyData (rhs);
}

void
RooUnfoldTUnfold::CopyData (const RooUnfoldTUnfold& rhs)
{
  tau_set=rhs.tau_set;
  _tau=rhs._tau;
  _reg_method=rhs._reg_method;
}


void
RooUnfoldTUnfold::Init()
{
    //Sets error matrix
    tau_set=false;
    _tau=0;
    _unf=0;
  GetSettings();
}

TObject*
RooUnfoldTUnfold::Impl()
{
    return _unf;
}


void
RooUnfoldTUnfold::Unfold()
{
    /* Does the unfolding. Uses the optimal value of the unfolding parameter unless a value has already been set using FixTau*/
       
  const TH2D* Hres=_res->Hresponse();
  if (_fail) return;
  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);
  TH2D* Hresc=CopyOverflow(Hres);
  TH2D* Hres_flipped=TransposeHist(Hresc);
  TH1::AddDirectory (oldstat);

  _unf= new TUnfold(Hres_flipped,TUnfold::kHistMapOutputHoriz,_reg_method);
  Int_t nScan=30;
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,23,0)  /* TUnfold v6 (included in ROOT 5.22) didn't have setInput return value */
  if(_unf->SetInput(_meas)>=10000) {
    cerr<<"Unfolding result may be wrong\n";
  }
#else
  _unf->SetInput(_meas);
#endif
  //_unf->SetConstraint(TUnfold::kEConstraintArea);
  if (!tau_set){
    iBest=_unf->ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
    _tau=_unf->GetTau();  // save value, even if we don't use it unless tau_set
    cout <<"tau= "<<_tau<<endl;
  }
  else{
    _unf->DoUnfold(_tau);
  }
  TH1* reco=_unf->GetOutput("_rec","reconstructed dist",0,0);
  _rec.ResizeTo (reco->GetNbinsX());
  for (int i=0;i<reco->GetNbinsX();i++){
    _rec(i)=(reco->GetBinContent(i+1));
  }
  delete Hresc;
  delete reco;
  delete Hres_flipped;
  _unfolded= true;
  _haveCov=  false;
}

void
RooUnfoldTUnfold::GetCov()
{
    //Gets Covariance matrix
    Int_t nt=_rec.GetNrows();
    if (_overflow){nt+=2;}
    if (!_unfolded) Unfold();
    if (_fail) return;
    TH2D* ematrix=_unf->GetEmatrix("ematrix","error matrix",0,0);
    _cov.ResizeTo (nt,nt);
    for (Int_t i= 0; i<nt; i++) {
        for (Int_t j= 0; j<nt; j++) {
            _cov (i,j)= ematrix->GetBinContent(i+1,j+1);
        }
    }
    delete ematrix;
    _haveCov= true;
}


TH2D*
RooUnfoldTUnfold::TransposeHist(const TH2D* h)
{
    //Returns the transpose of a matrix expressed as a TH2D
    //Inefficiencies are returned in the underflow bin
  Int_t nx= h->GetNbinsX(), ny= h->GetNbinsY();
  Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax();
  Double_t ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax();
  TH2D* hx= new TH2D ("h_flipped", h->GetTitle(), ny, ylo, yhi, nx, xlo, xhi);
  if (nx<ny){
    cerr<<"Warning: fewer x bins than y bins. Unfolding may not work correctly"<<endl;
  } 
  for (Int_t i= 1; i <=ny; i++) {
    Double_t sumineff=0;
    for (Int_t j= 1; j <=nx; j++) {
      hx->SetBinContent (j, i, h->GetBinContent (i, j));
      hx->SetBinError   (j, i, h->GetBinError   (i, j));
      sumineff+=h->GetBinContent(i,j); 
    }
    hx->SetBinContent(i,0,_res->Htruth()->GetBinContent(i)-sumineff);
  }
  return hx;
}

TH2D*
RooUnfoldTUnfold::CopyOverflow (const TH2D* h) const
{
  if (!_overflow) return (TH2D*)h->Clone();
  Int_t nx= h->GetNbinsX(), ny= h->GetNbinsX();
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

void 
RooUnfoldTUnfold::FixTau(Double_t tau)
{
    _tau=tau;
    tau_set=true;
}

void
RooUnfoldTUnfold::SetRegMethod(TUnfold::ERegMode regmethod)
{
    /*
    Decides the regularisation method.

      regemthod setting             regularisation
      ===========================   ==============
      TUnfold::kRegModeNone         none
      TUnfold::kRegModeSize         minimize the size of (x-x0)
      TUnfold::kRegModeDerivative   minimize the 1st derivative of (x-x0)
      TUnfold::kRegModeCurvature    minimize the 2nd derivative of (x-x0)
      kRegModeDerivative and kRegModeCurvature are the optimal settings for a 1D input
      kRegModeSize is optimal for a 2D distribution (but still does not produce great results).
     */
    _reg_method=regmethod;
}

void
RooUnfoldTUnfold::OptimiseTau()
{
    tau_set=false;
}

void
RooUnfoldTUnfold::GetSettings()
{
    _minparm=0;
    _maxparm=1;
    _stepsizeparm=1e-2;
    _defaultparm=2;
}
