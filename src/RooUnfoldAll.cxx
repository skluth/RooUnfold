//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldAll.cxx,v 1.6 2010-08-04 14:53:04 fwx38934 Exp $
//
// Description:
//      Unfolding errors class
//
// Author: Richard Claridge <richard.claridge@stfc.ac.uk>
//
//==============================================================================

//____________________________________________________________
/* BEGIN_HTML

<p> A graph drawing class to view the errors associated with an unfolding technique</p>
<p>Before these can be run, the RooUnfoldAll object must be created and the operation Plotting() run on the object in order to do the 
maths needed to plot these graphs. The object requires the number of iterations over which the errors are calculated and a RooUnfold object.</p>
<p>For each iteration each bin in the measured distribution is added to a random number from a gaussian with a width based on the error in that bin. This is then unfolded and the results plotted for each bin. The rms in each bin is then used as the spread of the values in every bin. This gives errors that are slightly larger than those returned by RooUnfold, but are a better representation of the spread in the data.</p> 
<p> If the true distribution is not known, the following can be returned:</p>
<ul>
<li> A graph of the errors from the unfolding (Unf_err())
<li> A graph of the errors due to the spread of the reconstructed points (Spread())
<li> An error matrix based on the spread of the reconstructed points (True_err())
</ul>
<p> If the true distribution is known then a plot of the chi squared values can also be returned (Chi2()).
 This requires the inclusion of the truth distribution and the error method on which the chi squared is based 
 (0 for a simple calculation, 1 or 2 for a method based on the covariance matrix, depending on the method used for calculation of errors.). </p>
<p>On some occasions the chi squared value can be very large. This is due to the covariance matrices being near singular and thus 
difficult to invert reliably. A warning will be displayed if this is the case. To plot the chi squared distribution use the option Draw("chi2"), to filter out the larger values use Draw("chi2","abs(chi2 < max") where max is the largest value to be included.</p> 
END_HTML */
/////////////////////////////////////////////////////////////////

#include "RooUnfoldAll.h"

#include <cfloat>
#include <iostream>
#include <cmath>

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TLine.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TAxis.h"

#include "RooUnfold.h"
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::fabs;

ClassImp (RooUnfoldAll);

RooUnfoldAll::RooUnfoldAll (int Nits,  const RooUnfold* unfold_in)
:iterations(Nits),unfold(unfold_in)
{
	h_err=0;
    h_err_res=0;
    h_err_res_sq=0;
    chi2=0; 
    All_hMeas(); 
}



RooUnfoldAll::~RooUnfoldAll()
{

  delete h_err;
  delete h_err_res;
  delete h_err_res_sq;
  delete chi2;  
}

void 
RooUnfoldAll::All_hMeas()
{
	//Gets graph size parameters//
	RooUnfold* u_c=unfold->Clone("clone");
	u_c->SetVerbose(unfold->verbose());
	u_c->SetNits(iterations);
	TH1* HR=u_c->Hreco();
	hMeas_const=unfold->Hmeasured();
	ntx=HR->GetNbinsX();
	xlo=HR->GetXaxis()->GetXmin();
	xhi=HR->GetXaxis()->GetXmax();
}


TNtuple*
RooUnfoldAll::Chi2(const TH1* hTrue,Int_t doerror)
{   
	if (!hTrue){
		cerr <<"Error: no truth distribution"<<endl;
	return 0;
	}
	/*Gets Chi squared values and returns an NTuple which can then be plotted. Prints error warning if chi squared 
values over 10^10 are returned. If Doerror=0 the chi squared value will come from a simple summation of the 
residuals. Doerror=1 will use the covariance matrix that comes from the unfolding. Doerror=2 will use the 
covariance matrix from True_err()*/
	double max=1e10;
	int odd_ch=0;
    chi2 = new TNtuple("chi2","chi2","chi2");
	for (int j=0; j<iterations;j++){
		TH1* hMeas_AR = dynamic_cast<TH1*>(hMeas_const->Clone ("Measured"));   hMeas_AR  ->SetTitle ("Measured");
	    TH1* hMeas=Add_Random(hMeas_AR);
	    RooUnfold* unfold_copy = unfold->Clone("unfold_toy");
		unfold_copy->Setup(unfold->response(),hMeas);
	    unfold_copy->SetVerbose(unfold->verbose());
	    unfold_copy->SetNits(iterations);
    	if (hTrue){
	    	double ch =unfold_copy->Chi2(hTrue,doerror);
	    	double f_ch = fabs(ch);
	    	chi2->Fill(ch);
	    	if (f_ch>=max){
	    		cerr<<"Large |chi^2| value: "<< ch << endl;
	    		odd_ch++;
	    	}
    	}
	}
    if (hTrue&&odd_ch!=0){
    	cout <<"There are " << odd_ch << " bins over "<<max <<endl;
    }
	chi2->SetFillColor(4);
	return chi2;
}

TMatrixD
RooUnfoldAll::True_err()
{
	//Returns a matrix of errors based on the spread of the reconstructed points//
	TMatrixD Error(ntx,ntx);
	for (int x=0;x<ntx;x++){
		Error(x,x)=h_err_res_sq->GetBinContent(x);
	}
	return Error;
}

TH1*
RooUnfoldAll::Spread(){
	//Returns a TH1D of the spread of the reconstructed points//
	h_err_res->SetMarkerColor(2);
	h_err_res->SetMarkerStyle(4);
	h_err_res->SetMinimum(0);
	return dynamic_cast<TH1D*>(h_err_res->Clone());
}

TH1* 
RooUnfoldAll::Unf_err(){
	//Returns a TH1D of the erros from the unfolding// 
	h_err = new TProfile ("h_err", "Unfolding errors",ntx,xlo,xhi);
	double dx=(xhi-xlo)/ntx;
	RooUnfold* unfold_copy = unfold->Clone("unfold_toy");
	unfold_copy->SetVerbose(unfold->verbose());
	unfold_copy->SetNits(iterations);
	TH1* hReco_= unfold_copy->Hreco(1);
	for (int i=0; i<ntx+2; i++) {    
    	if (hReco_->GetBinContent(i)!=0.0 || (hReco_->GetBinError(i)>0.0)) 
    	{
		Double_t u_error=hReco_->GetBinError(i); 
    	h_err->Fill(i*dx,u_error);
    	}
	}
	h_err->SetMarkerColor(4);
	h_err->SetMarkerStyle(24);
	h_err->SetMinimum(0);
	return dynamic_cast<TH1D*>(h_err->Clone());
}


void
RooUnfoldAll::Plotting()
{
	//Does the math//
	vector<TH1D*> graph_vector;    
    for (int i= 0; i<ntx+2; i++) {
	TString graph_title("Residuals at Bin ");
    graph_title+=i;
    TH1D* graph_name = new TH1D (graph_title,graph_title, 200,0,1000);
    graph_vector.push_back(graph_name);
    }
    h_err_res = new TH1D ("h_err_res", "Spread",ntx,xlo,xhi); 
    h_err_res_sq = new TH1D ("h_err_res_sq", "Spread^2",ntx,xlo,xhi);
	
	for (int j=0; j<iterations;j++){
		TH1* hMeas_AR = dynamic_cast<TH1*>(hMeas_const->Clone ("Measured"));   hMeas_AR  ->SetTitle ("Measured");
	    TH1* hMeas=Add_Random(hMeas_AR);
	    RooUnfold* unfold_copy = unfold->Clone("unfold_toy");
		unfold_copy->Setup(unfold->response(),hMeas);
	    unfold_copy->SetVerbose(unfold->verbose());
	    unfold_copy->SetNits(iterations);
	    hReco= unfold_copy->Hreco();
    	for (int i=0; i<ntx+2; i++) {    
    		if (hReco->GetBinContent(i)!=0.0 || (hReco->GetBinError(i)>0.0)) 
    		{
      		Double_t res= hReco->GetBinContent(i);
     		graph_vector[i]->Fill(res);  
        	}
    	}
	}
	for (unsigned int i=0; i<graph_vector.size(); i++){
		Double_t spr=(graph_vector[i]->GetRMS());
    	h_err_res_sq->SetBinContent(i,spr*spr);
    	h_err_res->SetBinContent(i,spr);
    }
    
    for (unsigned int i=0; i<graph_vector.size(); i++){
    	delete graph_vector[i];
    }
    
}



TH1*
RooUnfoldAll::Add_Random(TH1* hMeas_AR)
{
	//Adds a random number to the measured distribution before the unfolding//
	TH1* hMeas2=   dynamic_cast<TH1*>(hMeas_AR->Clone ("Measured"));   hMeas2  ->SetTitle ("Measured");
	for (Int_t i=0; i<hMeas_AR->GetNbinsX()+2 ; i++){
      		Double_t err=hMeas_AR->GetBinError(i);
      		Double_t new_x = hMeas_AR->GetBinContent(i) + gRandom->Gaus(0,err);
      		hMeas2->SetBinContent(i,new_x);
      		hMeas2->SetBinError(i,err);
    }
	return hMeas2;
}
