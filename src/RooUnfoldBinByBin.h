#ifndef ROOUNFOLDBINBYBIN_H_
#define ROOUNFOLDBINBYBIN_H_

#include "RooUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;

class RooUnfoldBinByBin : public RooUnfold {

public:


 RooUnfoldBinByBin(); // default constructor
  RooUnfoldBinByBin (const char*    name, const char*    title); // named constructor
  RooUnfoldBinByBin (const TString& name, const TString& title); // named constructor
  RooUnfoldBinByBin (const RooUnfoldBinByBin& rhs); // copy constructor
  virtual ~RooUnfoldBinByBin(); // destructor
  RooUnfoldBinByBin& operator= (const RooUnfoldBinByBin& rhs); // assignment operator
  virtual RooUnfoldBinByBin* Clone (const char* newname= 0) const;
RooUnfoldBinByBin (const RooUnfoldResponse* res, const TH1* meas, const char* name=0, const char* title=0);
  

protected:
virtual void Unfold();
virtual void GetCov();
virtual void GetSettings();

private:
  	int HresXbins;
  	TVectorD c_vector;
public:

  ClassDef (RooUnfoldBinByBin, 0) 
};

inline RooUnfoldBinByBin::RooUnfoldBinByBin()                                           : RooUnfold()           {Init();}
inline RooUnfoldBinByBin::RooUnfoldBinByBin (const char* name, const char* title)       : RooUnfold(name,title) {Init();}
inline RooUnfoldBinByBin::RooUnfoldBinByBin (const TString& name, const TString& title) : RooUnfold(name,title) {Init();}
inline RooUnfoldBinByBin& RooUnfoldBinByBin::operator= (const RooUnfoldBinByBin& rhs) {Assign(rhs); return *this;}


#endif /*ROOUNFOLDBINBYBIN_H_*/
