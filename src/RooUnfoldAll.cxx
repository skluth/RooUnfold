//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldAll.cxx,v 1.4 2010-07-22 16:15:30 fwx38934 Exp $
//
// Description:
//       Graph Drawing Class for use with RooUnfold.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Richard Claridge <richard.claridge@stfc.ac.uk>
//
//==============================================================================


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
	hMeas_const=unfold->Hmeasured();
	ntx=hMeas_const->GetNbinsX();
	xlo=hMeas_const->GetXaxis()->GetXmin();
	xhi=hMeas_const->GetXaxis()->GetXmax();
}


TNtuple*
RooUnfoldAll::Chi2()
{
	chi2->SetFillColor(4);
	return chi2;
}

TMatrixD
RooUnfoldAll::True_err()
{
	TMatrixD Error(ntx,ntx);
	for (int x=0;x<ntx;x++){
		Error(x,x)=h_err_res_sq->GetBinContent(x);
	}
	return Error;
}

TH1*
RooUnfoldAll::Spread(){
	h_err_res->SetMarkerColor(2);
	h_err_res->SetMarkerStyle(4);
	h_err_res->SetMinimum(0);
	return dynamic_cast<TH1D*>(h_err_res->Clone());
}

TH1* 
RooUnfoldAll::Unf_err(){
	h_err->SetMarkerColor(4);
	h_err->SetMarkerStyle(24);
	h_err->SetMinimum(0);
	return dynamic_cast<TH1D*>(h_err->Clone());
}


void
RooUnfoldAll::Plotting(const TH1 *hTrue)
{
	vector<TH1D*> graph_vector;    
    
    int odd_ch=0;
    double max=1e10;
    double dx=(xhi-xlo)/ntx;
    
    for (int i= 0; i<ntx+2; i++) {
	TString graph_title("Residuals at Bin ");
    graph_title+=i;
    TH1D* graph_name = new TH1D (graph_title,graph_title, 200,0,1000);
    graph_vector.push_back(graph_name);
    }
    
    chi2 = new TNtuple("chi2","chi2","chi2");
    h_err = new TProfile ("h_err", "Unfolding errors",ntx,xlo,xhi);
    h_err_res = new TH1D ("h_err_res", "Spread",ntx,xlo,xhi); 
    h_err_res_sq = new TH1D ("h_err_res_sq", "Spread^2",ntx,xlo,xhi);

	for (int j=0; j<iterations;j++){
		TH1* hMeas_AR = dynamic_cast<TH1*>(hMeas_const->Clone ("Measured"));   hMeas_AR  ->SetTitle ("Measured");
	    TH1* hMeas=Add_Random(hMeas_AR);
	    RooUnfold* unfold_copy = unfold->Clone("unfold_toy");
		unfold_copy->Setup(unfold->response(),hMeas);
	    unfold_copy->SetVerbose(unfold->verbose());
	    hReco= unfold_copy->Hreco();
    	for (int i=0; i<ntx+2; i++) {    
    		if (hReco->GetBinContent(i)!=0.0 || (hReco->GetBinError(i)>0.0)) 
    		{
      		Double_t res= hReco->GetBinContent(i);
     		graph_vector[i]->Fill(res);
     		Double_t u_error=hReco->GetBinError(i); 
    		h_err->Fill(i*dx,u_error);  
        	}
    	}
    	if (hTrue){
	    	double ch =unfold_copy->Chi2(hTrue);
	    	double f_ch = fabs(ch);
	    	chi2->Fill(ch);
	    	if (f_ch>=max){
	    		cerr<<"Large |chi^2| value: "<< ch << endl;
	    		odd_ch++;
	    	}
    	}
    	unfold_copy->Reset();
    	delete hMeas_AR;
    	delete hMeas;
    	
	}
	for (unsigned int i=0; i<graph_vector.size(); i++){
		Double_t spr=(graph_vector[i]->GetRMS());
    	h_err_res_sq->SetBinContent(i,spr*spr);
    	h_err_res->SetBinContent(i,spr);
    }
    for (unsigned int i=0; i<graph_vector.size(); i++){
    	delete graph_vector[i];
    }
    if (hTrue&&odd_ch!=0){
    	cout <<"There are " << odd_ch << " bins over "<<max <<endl;
    }
    
}



TH1*
RooUnfoldAll::Add_Random(TH1* hMeas_AR)
{
	TH1* hMeas2=   dynamic_cast<TH1*>(hMeas_AR->Clone ("Measured"));   hMeas2  ->SetTitle ("Measured");
	for (Int_t i=0; i<hMeas_AR->GetNbinsX()+2 ; i++){
      		Double_t err=hMeas_AR->GetBinError(i);
      		Double_t new_x = hMeas_AR->GetBinContent(i) + gRandom->Gaus(0,err);
      		hMeas2->SetBinContent(i,new_x);
      		hMeas2->SetBinError(i,err);
    }
	return hMeas2;
}
