#ifndef ROOUNFOLDDAGOSTINI_H_
#define ROOUNFOLDDAGOSTINI_H_

#include "RooUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;

class RooUnfoldDagostini : public RooUnfold {

public:
  RooUnfoldDagostini(); // default constructor
  RooUnfoldDagostini (const char*    name, const char*    title); // named constructor
  RooUnfoldDagostini (const TString& name, const TString& title); // named constructor
  RooUnfoldDagostini (const RooUnfoldDagostini& rhs); // copy constructor
  virtual ~RooUnfoldDagostini(); // destructor
  RooUnfoldDagostini& operator= (const RooUnfoldDagostini& rhs); // assignment operator
  virtual RooUnfoldDagostini* Clone (const char* newname= 0) const;
  RooUnfoldDagostini (const RooUnfoldResponse* res, const TH1* meas, Int_t niter= 4, const char* name=0, const char* title=0);

  void SetIterations (Int_t niter= 4)  { _niter= niter; }
  Int_t GetIterations() const { return _niter; }
  virtual void  SetRegParm (Double_t parm) { SetIterations(Int_t(parm+0.5)); }
  virtual Double_t GetRegParm() const { return GetIterations(); }

  virtual void Reset();

protected:
  void Assign (const RooUnfoldDagostini& rhs); // implementation of assignment operator
  virtual void Unfold();
  virtual void GetCov();
  virtual void GetSettings();

private:
  void Init();

protected:
  // instance variables
  Int_t _niter;
  
public:
  ClassDef (RooUnfoldDagostini, 0) 
};

inline RooUnfoldDagostini::RooUnfoldDagostini()                                           : RooUnfold()           {Init();}
inline RooUnfoldDagostini::RooUnfoldDagostini (const char* name, const char* title)       : RooUnfold(name,title) {Init();}
inline RooUnfoldDagostini::RooUnfoldDagostini (const TString& name, const TString& title) : RooUnfold(name,title) {Init();}
inline RooUnfoldDagostini& RooUnfoldDagostini::operator= (const RooUnfoldDagostini& rhs) {Assign(rhs); return *this;}

#endif /*ROOUNFOLDDAGOSTINI_H_*/
