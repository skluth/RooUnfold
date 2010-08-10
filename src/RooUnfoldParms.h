#ifndef ROOUNFOLDPARMS_H_
#define ROOUNFOLDPARMS_H_

#include "TNamed.h"

class TH1;
class RooUnfold;
class TProfile;

class RooUnfoldParms : public TNamed {
	public:
	RooUnfoldParms(const RooUnfold* unfold_in=0,Int_t err=1,Int_t its=500,const TH1* truth=0);
	virtual ~RooUnfoldParms();
	TProfile* GetChi2();
	TProfile* GetErr();
	TProfile* GetRes();
	TH1* GetRMSSpread();
	const RooUnfold* unfold; // Input object from RooUnfold
	Int_t doerror; // Set error calculation method
	Int_t Nits; // Number of iterations
	const TH1* hTrue; // Truth Distribution
	void Set_Min(double min);
	void Set_Max(double max);
	void Set_Stepsize(double size);
	
	protected:
	bool _done_math;
	TH1* hrms; // Output plot
	TProfile* hch2; // Output plot
	TProfile* herr; // Output plot
	TProfile* hres; // Output plot
	void DoMath();
	void Init();
	Double_t _maxparm; //Maximum parameter
	Double_t _minparm; //Minimum parameter
	Double_t _stepsizeparm; //Step size
public:
	ClassDef (RooUnfoldParms, 0)
};
#endif /*ROOUNFOLDPARMS_H_*/
