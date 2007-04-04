//--------------------------------------------------------------------
//
// $Id: RooUnfoldBayesImpl.h,v 1.1.1.1 2007-04-04 21:27:02 adye Exp $
//
//--------------------------------------------------------------------
#ifndef ROOUNFOLDBAYESIMPL_HH
#define ROOUNFOLDBAYESIMPL_HH

#include "TNamed.h"
#include <vector>
using std::vector;

class Array2D;
class TH1;

class RooUnfoldBayesImpl : public TNamed {

 public:
  enum UnfoldType {BinByBin=0, Inversion=1, Iterative=2, Blobel=3, SVD=4, MEM=5};

private:
//Prevent default methods being used (they won't work).
  RooUnfoldBayesImpl& operator=(const RooUnfoldBayesImpl& rhs);
  RooUnfoldBayesImpl(const RooUnfoldBayesImpl &unfold);

  void init();

  // return a unique bin number from 3 dimensional space.
  Int_t getTruthBin(Double_t x, Double_t y=0.0, Double_t z=0.0) const;
  Int_t getRecoBin (Double_t x, Double_t y=0.0, Double_t z=0.0) const;
  Int_t getBin     (Double_t x, Double_t xmin, Double_t xmax, Int_t nx) const;
  Double_t deltaM  (Int_t k, Int_t i, Int_t r, Int_t u) const;
  Int_t getCovariance(const vector<Double_t>& effects);
  Int_t getCovarianceBinByBin(vector<Double_t> effects);
  Double_t sum(const vector<Double_t>& x) const;
  Int_t smooth(const vector<Double_t>& PbarCi) const;
  Bool_t proceed() const;
  Double_t getnbarCi(const vector<Double_t>& effects,
                     vector<Double_t> &causes) const;
  Int_t getBin(Int_t index, Bool_t truth, vector<Int_t>& coords) const;

  // Truth variables
  vector<Int_t> _nt;      // number of bins per dimension
  vector<Double_t> _tmin; // minimum allowed value per dimension
  vector<Double_t> _tmax; // maximum allowed value per dimension
  Int_t _nc;              // number of causes (product of bins per dimension)
  Int_t _ndims;           // number/dimension of truth variables
  Int_t _nOutGeneratedRange; // number of truth training events outside dimension range

  // reconstructed variables
  vector<Int_t> _nr;      // number of bins per dimension
  vector<Double_t> _rmin; // minimum allowed value per dimension
  vector<Double_t> _rmax; // maximum allowed value per dimension
  Int_t _ne;              // number of effects (product of bins per dimension)
  Int_t _mdims;           // number/dimension of measured variables
  Int_t _nOutMeasuredRange; // number of measured training events outside dimension range

  Int_t _nOutEffectsRange; // number of measured events outside dimension range
  Double_t _neffects;      // number of measured events inside dimension range
  Double_t _ncauses;       // Estimated number of unfolded events
  Double_t _nCausesError;  // Estimated error on number of unfolded events

  Double_t _nbartrue; // best estimate of number of true events

  //  TVectorD *_PEj; // Probability of Effect E_j
  vector<Double_t> _nEj;   // Number of measured events from Effect E_j
  vector<Double_t> _nEstj; // Number of measured events from Effect E_j
  vector<Double_t> _PCi; // Probability of cause C_i
  vector<Double_t> _nCi; // Number of true events from cause C_i
  vector<Double_t> _efficiencyCi; // efficiency for detecting cause C_i
  vector<Double_t> _causes; // vector of measured variables to be unfolded
  vector<Double_t> _causesError; // vector of error on measured variables to be unfolded

  Array2D *_Sij; // Smearing matrix
  Array2D *_Nij; // mapping of causes to effects
  Array2D *_Mij; // unfolding matrix
  Array2D *_Vij; // covariance matrix

  Int_t _fDebug; // debug level
  Bool_t _initialised_vectors; // true after calling build
  Bool_t _trained; // true after training
  UnfoldType _UnfoldType; // which unfolding method to use

 public:

  RooUnfoldBayesImpl();

  RooUnfoldBayesImpl(const char* name, const char* title);

  ~RooUnfoldBayesImpl();

  // Initialisation functions
  Int_t build(Int_t ndims, const vector<Int_t>& nx, const vector<Double_t>& xmin, const vector<Double_t>& xmax,
              Int_t mdims, const vector<Int_t>& mx, const vector<Double_t>& ymin, const vector<Double_t>& ymax);

  // mapping of measured to generated events used for training
  Int_t accumulate(const vector<Double_t>& xtrue, const vector<Double_t>& xmeas);

  Int_t accumulate(const vector<Double_t>& xmeas);

  // fill arrays directly instead of using accumulate()
  Int_t setupTrain(const vector<Double_t>& nCi, const vector<Double_t>& nEj, const Array2D& Nij);
  Int_t setupUnfold(const vector<Double_t>& nEstj);

  // clear arrays
  void reset();
  Int_t clear();

  // set the debug level
  Int_t setDebug(Int_t Db);

  //
  Int_t train(Int_t iterations=3, Bool_t smoothit = false);
  Int_t trainBinByBin(Bool_t smoothit = false);

  //
  Int_t info(Int_t level=0) const; // print some useful info
  //
  Int_t unfold(vector<Double_t>& causes); // unfold vector of causes after training
  Int_t unfoldBinByBin(vector<Double_t>& causes); // unfold vector of causes after training

  Double_t getError(); // Estimated error on number of unfolded events

  // returns histogram of unfolded causes
  TH1* histCauses (const char* name= "causes", const char* title= "Unfolded results") const;

  const vector<Double_t>& truth()  const { return _nCi;    } // vector of training causes
  const vector<Double_t>& reco()   const { return _nEj;    } // vector of training effects
  const vector<Double_t>& input()  const { return _nEstj;  } // vector of measured effects
  const vector<Double_t>& output() const { return _causes; } // return unfold vector of causes
  const Array2D&      covariance() const { return *_Vij;   } // return unfold covariance matrix
  Double_t                error()  const { return _nCausesError; }  // return estimated error on number of unfolded events
  
  ClassDef(RooUnfoldBayesImpl,1) // Bayes Unfolding Algorithms
};

#endif
