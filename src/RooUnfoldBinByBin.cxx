//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBinByBin.cxx,v 1.10 2010-08-19 12:31:07 fwx38934 Exp $
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
END_HTML */

/////////////////////////////////////////////////////////////

#include "RooUnfoldBinByBin.h"

#include <iostream>
#include <iomanip>
#include <cmath>

#include "TNamed.h"
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
using std::sqrt;

ClassImp (RooUnfoldBinByBin);

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldBinByBin& rhs)
  : RooUnfold (rhs)
{
  
}

RooUnfoldBinByBin::RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas, 
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title)
{
  
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



void
RooUnfoldBinByBin::Unfold()
{
	const TH2D* Hres=_res->Hresponse();
	HresXbins=Hres->GetNbinsX();
  	if (_overflow){
	  	HresXbins+=2;
    }
  	_rec.ResizeTo(HresXbins);
  	double c;
  	c_vector.ResizeTo(HresXbins);
  	for (int i=0; i<HresXbins;i++){ 	
		if(_overflow){
			c=(_res->Htruth()->GetBinContent(i))/(_res->Hmeasured()->GetBinContent(i));
  			if (_res->Hmeasured()->GetBinContent(i)==0){
  				_rec(i)=0;
  				c_vector(i)=0;
  			}
  			else{
  				c_vector(i)=c;
  				_rec(i)= _meas->GetBinContent(i)*c;
  			}
  		}
  		else{
  			//fiddling needed to exclude underflow bin (still included)
  			c=(_res->Htruth()->GetBinContent(i+1))/(_res->Hmeasured()->GetBinContent(i+1));
  			if (_res->Hmeasured()->GetBinContent(i+1)==0){
  				_rec(i)=0;
  				c_vector(i)=0;
  			}
  			else{
  				c_vector(i)=c;
  				_rec(i)= _meas->GetBinContent(i+1)*c;
  			}
  		}
  	}
  	_unfolded= true;
  	_haveCov=  false;
}

void
RooUnfoldBinByBin::GetCov()
{
	if(!_unfolded){Unfold();}
	_cov.ResizeTo(HresXbins,HresXbins);
	for (int i=0; i<HresXbins;i++){
		if(_overflow){
			double c=c_vector(i);
			if (_res->Hmeasured()->GetBinContent(i)==0){
				_cov(i,i)=0;
			}
			else{
				_cov(i,i)=c*c*_meas->GetBinContent(i);
			}
		}
		else{
			double c=c_vector(i);
			if (_res->Hmeasured()->GetBinContent(i+1)==0){
				_cov(i,i)=0;
			}
			else{
				_cov(i,i)=c*c*_meas->GetBinContent(i+1);
			}
		}
			
	}
	_haveCov= true;
}

void
RooUnfoldBinByBin::GetSettings(){
	_minparm=0;
	_maxparm=0;
	_stepsizeparm=0;
	_defaultparm=0;
}
