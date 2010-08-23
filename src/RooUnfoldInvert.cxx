//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldInvert.cxx,v 1.4 2010-08-23 21:38:13 adye Exp $
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
}

RooUnfoldInvert::RooUnfoldInvert (const RooUnfoldResponse* res, const TH1* meas,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title)
{
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
}

void
RooUnfoldInvert::Unfold()
{
       
  const TH2D* Hres=_res->Hresponse();
  const TH1* Hmeas=_meas;
  int HresXbins=Hres->GetNbinsX();
  int HresYbins=Hres->GetNbinsY();
  int HmeasXbins=Hmeas->GetNbinsX();
  
  TMatrixD* Hres_m;
  if (_overflow){Hres_m=RooUnfoldResponse::H2M(Hres,HresXbins,HresYbins,_res->Htruth(),1);}
  else {Hres_m=RooUnfoldResponse::H2M(Hres,HresXbins,HresYbins,_res->Htruth(),0);}
  if (_overflow){
    HresXbins+=2;
    HresYbins+=2;
    HmeasXbins+=2;
  }
  TVectorD* Hmeas_mp=RooUnfoldResponse::H2V(Hmeas,HmeasXbins);
  Hmeas_m.ResizeTo(HmeasXbins);
  Hmeas_m=(*Hmeas_mp);
  TDecompSVD svd(*Hres_m);
  if (svd.Condition()<0){
    cerr <<"Error: bad condition= "<<svd.Condition()<<endl;
  }
  Hres_i.ResizeTo(HresYbins,HresXbins);
  Hres_i=svd.Invert();
  _rec.ResizeTo(Hres_i.GetNrows());
  _rec=Hres_i*Hmeas_m;
  _unfolded= true;
  _haveCov=  false;
}

void
RooUnfoldInvert::GetCov()
{
    if (!_unfolded){Unfold();}
    _cov.ResizeTo(Hres_i.GetNrows(),Hres_i.GetNcols());
    for (int i=0;i<Hres_i.GetNrows();i++){
        for (int j=0;j<Hres_i.GetNrows();j++){
            for (int k=0; k<Hres_i.GetNcols();k++){
                _cov(i,j)+=(Hres_i(i,k)*Hres_i(j,k)*Hmeas_m(k));
            }
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
