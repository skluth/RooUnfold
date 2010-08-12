#ifndef ROOUNFOLDTUNFOLD_H_
#define ROOUNFOLDTUNFOLD_H_

#include "RooUnfold.h"
#include "TUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;

class RooUnfoldTUnfold : public RooUnfold {

public:

  RooUnfoldTUnfold(); // default constructor
  RooUnfoldTUnfold (const char*    name, const char*    title); // named constructor
  RooUnfoldTUnfold (const TString& name, const TString& title); // named constructor
  RooUnfoldTUnfold (const RooUnfoldTUnfold& rhs); // copy constructor
  virtual ~RooUnfoldTUnfold(); // destructor
  RooUnfoldTUnfold& operator= (const RooUnfoldTUnfold& rhs); // assignment operator
  virtual RooUnfoldTUnfold* Clone (const char* newname= 0) const;
  RooUnfoldTUnfold (const RooUnfoldResponse* res, const TH1* meas,Int_t reg=2, 
                const char* name= 0, const char* title= 0);
	void Reset();
	
	TObject* Impl();
	void FixTau(Double_t tau);
	void OptimiseTau();
	virtual void SetRegParm(Double_t parm){FixTau(parm);}
	Double_t GetTau() const { return _tau;    }
	virtual Double_t GetRegParm() const { return GetTau(); }
	void SetRegMethod (Int_t regmethod);
	Int_t GetRegMethod() const {return _reg_method;}

protected:
  void Init();
  void Destroy();
  virtual void Unfold();
  virtual void GetCov();
  void Assign   (const RooUnfoldTUnfold& rhs); // implementation of assignment operator
  void CopyData (const RooUnfoldTUnfold& rhs);
  virtual void GetSettings();

	
private:
  Int_t _reg_method; //Regularisation method 
  TUnfold* _unf; //TUnfold object
  TH2D* TransposeHist(const TH2D* Hres);
  TH2D* CopyOverflow(const TH2D* h)const;
  Bool_t tau_set;
  Double_t _tau;
  
public:

  ClassDef (RooUnfoldTUnfold, 0) 
};

// Inline method definitions

inline RooUnfoldTUnfold::RooUnfoldTUnfold()                                           : RooUnfold()           {Init();}
inline RooUnfoldTUnfold::RooUnfoldTUnfold (const char* name, const char* title)       : RooUnfold(name,title) {Init();}
inline RooUnfoldTUnfold::RooUnfoldTUnfold (const TString& name, const TString& title) : RooUnfold(name,title) {Init();}
inline RooUnfoldTUnfold& RooUnfoldTUnfold::operator= (const RooUnfoldTUnfold& rhs) {Assign(rhs); return *this;}
inline RooUnfoldTUnfold::~RooUnfoldTUnfold() {Destroy();}

#endif /*ROOUNFOLDTUNFOLD_H_*/
