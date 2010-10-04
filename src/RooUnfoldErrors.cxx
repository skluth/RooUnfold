//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
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
<p>Before these can be run, the RooUnfoldErrors object must be created and the operation Plotting() run on the object in order to do the 
maths needed to plot these graphs. The object requires the number of toys over which the errors are calculated and a RooUnfold object.</p>
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

#include "RooUnfoldErrors.h"

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
#include "RooUnfoldResponse.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::fabs;

ClassImp (RooUnfoldErrors);

RooUnfoldErrors::RooUnfoldErrors (int NToys,  RooUnfold* unfold_in, const TH1* Truth)
:toys(NToys),unfold(unfold_in),hTrue(Truth)
{
    h_err=0;
    h_err_res=0;
    hchi2=0; 
    GraphParameters(); 
    CreatePlots();
}



RooUnfoldErrors::~RooUnfoldErrors()
{

  delete h_err;
  delete h_err_res;
  delete hchi2;  
}

void 
RooUnfoldErrors::GraphParameters()
{
    //Gets graph size parameters//
    const TH1* HR=unfold->response()->Htruth();
    ntx=HR->GetNbinsX();
    xlo=HR->GetXaxis()->GetXmin();
    xhi=HR->GetXaxis()->GetXmax();
}


TNtuple*
RooUnfoldErrors::Chi2()
{   
    //Returns TNtuple of chi squared values. 
    hchi2->SetFillColor(4);
    return hchi2;
}

TH1*
RooUnfoldErrors::RMSResiduals(){
    //Returns a TH1D of the spread of the reconstructed points//
    h_err_res->SetMarkerColor(2);
    h_err_res->SetMarkerStyle(4);
    h_err_res->SetMinimum(0);
    return dynamic_cast<TH1D*>(h_err_res->Clone());
}

TH1* 
RooUnfoldErrors::UnfoldingError(){
    //Returns a TH1D of the errors from the unfolding// 
    h_err->SetMarkerColor(4);
    h_err->SetMarkerStyle(24);
    h_err->SetMinimum(0);
    return dynamic_cast<TH1D*>(h_err->Clone());
}


void
RooUnfoldErrors::CreatePlots()
{
    /*Gets the values for plotting. Uses the Runtoy method from RooUnfold to get plots to analyse for
    spread and error on the unfolding. Can also give values for a chi squared plot if a truth distribution is known*/
    if (!hTrue){
        cerr <<"Error: no truth distribution"<<endl;
    }
    TH1* parms=unfold->Runtoy(RooUnfold::kCovariance);
    if (parms->GetDimension()!=1){
        parms=RooUnfoldResponse::H2H1D (parms, ntx);
        ntx=parms->GetNbinsX();
        xhi=1;
        xlo=0;
    }
    delete parms;
    
    double dx=(xhi-xlo)/ntx;
    h_err = new TProfile ("h_err", "Unfolding errors",ntx,xlo,xhi);
    double max=1e10;
    int odd_ch=0;
    hchi2 = new TNtuple("chi2","chi2","chi2");
    vector<TH1D*> graph_vector;
    for (int a=0; a<ntx; a++) {
    TString graph_title("Residuals at Bin ");
    graph_title+=a;
    TH1D* graph_name = new TH1D (graph_title,graph_title, 2000,0,10000);
    graph_vector.push_back(graph_name);
    }
    
    h_err_res = new TH1D ("h_err_res", "Spread",ntx,xlo,xhi); 
    
    for (int k=0; k<toys;k++){  
        double chi2;
        TH1* hReco_;
        TH1* hReco=unfold->Runtoy(RooUnfold::kCovariance,&chi2,hTrue);
        if (hReco->GetDimension()!=1){
            hReco_=RooUnfoldResponse::H2H1D (hReco, ntx);
        }
        else{
            hReco_=hReco;
        }
        for (int i=0; i<ntx; i++) {    
            Double_t res= hReco_->GetBinContent(i);
            graph_vector[i]->Fill(res);
            Double_t u_error=hReco_->GetBinError(i+1); 
            h_err->Fill(i*dx,u_error);
            } 
        if (hTrue){
            hchi2->Fill(chi2);
            if (TMath::Abs(chi2)>=max && unfold->verbose()>=1){
                cerr<<"Large |chi^2| value: "<< chi2 << endl;
                odd_ch++;
            }
        }
    delete hReco_;
    }     
    for (unsigned int i=0; i<graph_vector.size(); i++){
        Double_t spr=(graph_vector[i]->GetRMS());
        h_err_res->SetBinContent(i,spr);
    }
    
    for (unsigned int i=0; i<graph_vector.size(); i++){
        delete graph_vector[i];
    }
    
    if (odd_ch){
        cout <<"There are " << odd_ch << " bins over outside the range of 0 to "<<max <<endl;
    }

}

