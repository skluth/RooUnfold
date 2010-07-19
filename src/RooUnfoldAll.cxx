//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldAll.cxx,v 1.2 2010-07-19 21:45:10 adye Exp $
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

ClassImp (RooUnfoldAll);

RooUnfoldAll::RooUnfoldAll (int Nits,const TH1* h_in,  const RooUnfold* unfold_in)
:iterations(Nits), hTrue(h_in),unfold(unfold_in)
{
	All_hTrue();
	All_hMeas();
	Plotting();
}


RooUnfoldAll::~RooUnfoldAll()
{
}

void 
RooUnfoldAll::All_hTrue()
{
	ntx=hTrue->GetNbinsX();
	xlo=hTrue->GetXaxis()->GetXmin();
	xhi=hTrue->GetXaxis()->GetXmax();
}

void 
RooUnfoldAll::All_hMeas()
{
	hMeas_const=unfold->Hmeasured();
	nmx=hMeas_const->GetNbinsX();
}


TNtuple*
RooUnfoldAll::Chi2()
{
	chi2->SetFillColor(4);
	return chi2;
}

void
RooUnfoldAll::Plotting()
{
	TNtuple* chi2_wip = new TNtuple("chi2","chi2","chi2");
	vector<TH1D*> graph_vector;    
    
    for (int i= 0; i<ntx+2; i++) {
	TString graph_title("Residuals at Bin ");
    graph_title+=i;
    TH1D* graph_name = new TH1D (graph_title,graph_title, 200,-150,150);
    graph_vector.push_back(graph_name);
    }
    double dx=(xhi-xlo)/ntx;
    TProfile* h_err_wip = new TProfile ("h_err", "Unfolding errors",ntx,xlo,xhi);
    TH1D* h_err_res_wip = new TH1D ("h_err_res", "Spread",ntx,xlo,xhi); 
	
	for (int j=0; j<iterations;j++){
	TH1* hMeas_AR = dynamic_cast<TH1*>(hMeas_const->Clone ("Measured"));   hMeas_AR  ->SetTitle ("Measured");
    TH1* hMeas=Add_Random(hMeas_AR);
    RooUnfold* unfold_copy = unfold->Clone("unfold_toy");
	unfold_copy->Setup(unfold->response(),hMeas);
        unfold_copy->SetVerbose(unfold->verbose());
    hReco= unfold_copy->Hreco();
    	for (int i=0; i<ntx+2; i++) {    
    		if ((hReco->GetBinContent(i)!=0.0 || (hReco->GetBinError(i)>0.0)) &&
        	(hTrue->GetBinContent(i)!=0.0 || (hTrue->GetBinError(i)>0.0))) {
      		Double_t res= hReco->GetBinContent(i) - hTrue->GetBinContent(i);
      		TH1D* projection=graph_vector[i];
      		projection->Fill(res);
     		graph_vector[i]=projection; 
     		Double_t u_error=hReco->GetBinError(i); 
    		h_err_wip->Fill(i*dx,u_error);  	
        	}
    	}
    	chi2_wip->Fill(unfold_copy->Chi2(hTrue));
    	unfold_copy->Reset();
    	delete hMeas_AR;
    	delete hMeas;
    	delete hReco;
	}
	
	
	for (unsigned int i=0; i<graph_vector.size(); i++){
    	h_err_res_wip->SetBinContent(i,graph_vector[i]->GetRMS());
    }
    h_err_res = h_err_res_wip;
    chi2=chi2_wip;
    h_err=h_err_wip;
}

TH1*
RooUnfoldAll::Spread(){
	h_err_res->SetMarkerColor(2);
	h_err_res->SetMarkerStyle(4);
	return h_err_res;
}

TH1* 
RooUnfoldAll::Unf_err(){
	h_err->SetMarkerColor(4);
	h_err->SetMarkerStyle(24);
	return h_err;
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
