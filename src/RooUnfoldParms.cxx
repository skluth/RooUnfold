//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldParms.cxx,v 1.2 2010-08-04 14:53:04 fwx38934 Exp $
//
// Description:
//      Optimisation of regularisation parameter class
//
// Author: Richard Claridge <richard.claridge@stfc.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML

<p>Allows user to optimise the regularisation parameter of an unfolding technique.</p>
<p>If the true distribution is unknown, a plot of the rms of the errors is returned for each regularisation parameter (GetErr()).</p>
<p>If the true distribution is known, the following plots can be returned:
<ul>
<li> Chi squared values vs regularisation parameter (GetChi2())
<li> RMS of the residuals given by the true and the unfolded distrbutions vs regularisation parameter (GetRes())
<li> RMS spread of the residuals vs regularisation parameter (GetRMSSpread())
</ul>
<p>The numerical handling of the input data is done automatically on creating the RooUnfoldParms object. 
This object requires the maximum regularisation parameter and the unfolded distributions (RooUnfold Object). 
The method of getting errors, the number of iterations the error calculation is run over (default 500) and the truth distribution can also be set, 
these are only relevant if the true distribution is already known</p>
<p>For each regparm value, the measured distribution is unfolded. For each unfolded distribution residuals are plotted and rms found for the 
rms spread. The sum of the residuals over the whole distribution are calculated,divided by the number of bins and then rooted in order to 
return an rms. The chi squared values are calculated using the chi2() method in RooUnfold.</p>

 END_HTML */
////////////////////////////////////////////////////////////////

#include "RooUnfoldParms.h"

#include <cfloat>
#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TProfile.h"
#include "RooUnfold.h"
#include "TRandom.h"
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

ClassImp (RooUnfoldParms);

RooUnfoldParms::RooUnfoldParms(Int_t reg ,const RooUnfold* unfold_in,Int_t err,Int_t its,const TH1* truth)
:regparm(reg),unfold(unfold_in),doerror(err),Nits(its),hTrue(truth)
{
	herr=0;
	hch2=0;
	hres=0;  
	hrms=0;
	DoMath();
}

RooUnfoldParms::~RooUnfoldParms()
{
  delete herr;
  delete hch2;
  delete hres;  
  delete hrms;
}

TProfile*
RooUnfoldParms::GetChi2()
{
	/*Returns TProfile of Chi squared values for each regularisation parameter.
	 Requires a known truth distribution*/
	hch2->SetMarkerColor(2);
	hch2->SetMarkerStyle(4);
	hch2->SetMinimum(0);
	return dynamic_cast<TProfile*>(hch2->Clone());
}
TProfile*
RooUnfoldParms::GetErr()
{
	//Returns TProfile of errors for each regularisation parameter//
	herr->SetMarkerColor(4);
	herr->SetMarkerStyle(4);
	herr->SetMinimum(0);
	return dynamic_cast<TProfile*>(herr->Clone());
}

TProfile*
RooUnfoldParms::GetRes()
{
	/*Returns TProfile of RMS Residuals for each regularisation parameter.
	 Requires a known truth distribution*/
	hres->SetMarkerColor(8);
	hres->SetMarkerStyle(4);
	hres->SetMinimum(0);
	return dynamic_cast<TProfile*>(hres->Clone());
}

TH1*
RooUnfoldParms::GetRMSSpread()
{
	/*Returns TH1D of RMS spread of Residuals for each regularisation parameter.
	 Requires a known truth distribution*/
	hrms->SetMarkerColor(1);
	hrms->SetMarkerStyle(4);
	hrms->SetMinimum(0);
	return dynamic_cast<TH1D*>(hrms->Clone());
}

void
RooUnfoldParms::DoMath()
{
	/*Does all the data handling. Called in constructor.*/
	cout<<"Doing parms"<<endl;
	RooUnfold* u_temp = unfold->Clone("unfold_toy");
	u_temp->SetVerbose(unfold->verbose());
	u_temp->SetNits(Nits);
	u_temp->Setup(unfold->response(),unfold->Hmeasured());
	Double_t _maxparm=u_temp->Get_maxparm();
	Double_t _minparm=u_temp->Get_minparm();
	Double_t _stepsizeparm=u_temp->Get_stepsizeparm();
	Int_t nobins=Int_t((_maxparm-_minparm)/_stepsizeparm);
	
	vector<TH1D*> graph_vector;    
    Double_t graphage=_minparm;
    
    while(graphage<=_maxparm) {
	TString graph_title("Residuals at k= ");
    graph_title+=graphage;
    TH1D* graph_name = new TH1D (graph_title,graph_title, 200,0,1000);
    graph_vector.push_back(graph_name);
    graphage+=_stepsizeparm;
    }
    
    Double_t xlo=_minparm;
	Double_t xhi=_maxparm;
	hch2=new TProfile("hch2","chi^2 vs regparm",nobins,xlo,xhi);
	herr=new TProfile("herr","Error(squared) vs regparm",nobins,xlo,xhi);
	hres=new TProfile("hres","rms vs regparm",nobins,xlo,xhi);
	hrms=new TH1D("hrms","rms spread of residuals",nobins,xlo,xhi);
    Int_t bins=0;
    Int_t gvl=0;
    delete u_temp;
	for (Double_t k=_minparm;k<=_maxparm;k+=_stepsizeparm)
	{
		TH1* hMeas = dynamic_cast<TH1*>(unfold->Hmeasured()->Clone ("Measured"));
		RooUnfold* unf = unfold->Clone("unfold_toy");
		unf->Setup(unfold->response(),hMeas);
	    unf->SetVerbose(unfold->verbose());
	    unf->SetRegParm(k);
	    unf->SetNits(Nits);
		Double_t sq_err_tot=0;
		TH1* hReco=unf->Hreco(1);
		bins=hReco->GetXaxis()->GetNbins(); 
		for (int i=0;i<bins;i++)
		{
			sq_err_tot+=hReco->GetBinError(i);
		}
		herr->Fill(k,sqrt(sq_err_tot));
		if (hTrue)
		{	
			Double_t rsqt=0;	
			for (int i=0;i<bins;i++){
				if (hReco->GetBinContent(i)!=0.0 || (hReco->GetBinError(i)>0.0)) 
	    		{
					Double_t res=hReco->GetBinContent(i) - hTrue->GetBinContent(i);
					Double_t rsq=res*res;
					rsqt+=rsq;
					graph_vector[gvl]->Fill(sqrt(rsq));
	    		}
			}
			double chi2=unf->Chi2(hTrue,doerror);
			if (chi2<=1e10){
				hch2->Fill(k,chi2);
			}
		hres->Fill(k,(sqrt(rsqt/bins)));	
		}
		gvl++;
		delete unf;
	}
	Double_t bn=_minparm;
	for (unsigned int i=0; i<graph_vector.size(); i++){
		Double_t spr=(graph_vector[i]->GetRMS());
    	hrms->Fill(bn,spr/sqrt(bins));
    	bn+=_stepsizeparm;
	}
   	for (unsigned int i=0; i<graph_vector.size(); i++){
    	delete graph_vector[i];
    }
    
}

