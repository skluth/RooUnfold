#ifndef ROOUNFOLDPARMS_H_
#define ROOUNFOLDPARMS_H_

#include "TNamed.h"

class TH1;
class RooUnfold;
class TProfile;

class RooUnfoldParms : public TNamed {
	public:
	RooUnfoldParms(Int_t reg=1,const RooUnfold* unfold_in=0,Int_t err=1,Int_t its=500,const TH1* truth=0);
	virtual ~RooUnfoldParms();
	TProfile* GetChi2();
	TProfile* GetErr();
	TProfile* GetRes();
	TH1* GetRMSSpread();
	Int_t regparm; // Regularisation Parameter
	const RooUnfold* unfold; // Input object from RooUnfold
	Int_t doerror; // Set error calculation method
	Int_t Nits; // Number of iterations
	const TH1* hTrue; // Truth Distribution
	
	protected:
	TH1* hrms; // Output plot
	TProfile* hch2; // Output plot
	TProfile* herr; // Output plot
	TProfile* hres; // Output plot
	TH1* htau;
	void DoMath();
};
#endif /*ROOUNFOLDPARMS_H_*/
