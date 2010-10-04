//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldInvert.cxx,v 1.6 2010-09-10 17:14:34 adye Exp $
//
// Description:
//      Unfolding class using inversion of the response matrix. This does not produce
//      good results and is designed to illustrate the need for more sophisticated
//      unfolding techniques
//
// Authors: Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML
<p>The simplest method of unfolding works by simply inverting the response matrix.</p> 
<p>This is not accurate for small matrices and produces inaccurate unfolded distributions.</p>
<p>The inversion method is included largely to illustrate the necessity of a more effective method of unfolding</p>
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldInvert.h"

#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"

#include "RooUnfoldResponse.h"

using std::cout;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldInvert);

RooUnfoldInvert::RooUnfoldInvert (const RooUnfoldInvert& rhs)
  : RooUnfold (rhs)
{
  Init();
}

RooUnfoldInvert::RooUnfoldInvert (const RooUnfoldResponse* res, const TH1* meas,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title)
{
  Init();
}

RooUnfoldInvert*
RooUnfoldInvert::Clone (const char* newname) const
{
  RooUnfoldInvert* unfold= new RooUnfoldInvert(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}


RooUnfoldInvert::~RooUnfoldInvert()
{
  delete _svd;
}

void
RooUnfoldInvert::Init()
{
  _svd= 0;
  GetSettings();
}

void
RooUnfoldInvert::Reset()
{
  delete _svd;
  Init();
  RooUnfold::Reset();
}

TDecompSVD*
RooUnfoldInvert::Impl()
{
  return _svd;
}

void
RooUnfoldInvert::Unfold()
{
  _svd= new TDecompSVD (_res->Mresponse());
  if (_svd->Condition()<0){
    cerr <<"Warning: response matrix bad condition= "<<_svd->Condition()<<endl;
  }

  _rec.ResizeTo(_nm);
  _rec= Vmeasured();
  Bool_t ok= _svd->Solve (_rec);
  if (!ok) {
    cerr << "Response matrix Solve failed" << endl;
    return;
  }

  _unfolded= true;
  _haveCov=  false;
}

void
RooUnfoldInvert::GetCov()
{
    Bool_t ok;
    TMatrixD resinv(_nt,_nm);
    resinv=_svd->Invert(ok);
    if (!ok) {
      cerr << "response matrix inversion failed" << endl;
      return;
    }

    const TVectorD& vmeasured= Vmeasured();
    _cov.ResizeTo(_nt,_nt);
    for (int i=0;i<_nt;i++){
        for (int j=0;j<=i;j++){
            Double_t c= 0;
            for (int k=0; k<_nm;k++){
                c += resinv(i,k) * resinv(j,k) * vmeasured(k);
            }
            _cov(i,j)= c;
            _cov(j,i)= c;
        }
    }
    _haveCov= true;
}

void
RooUnfoldInvert::GetSettings(){
    _minparm=0;
    _maxparm=0;
    _stepsizeparm=0;
    _defaultparm=0;
}
