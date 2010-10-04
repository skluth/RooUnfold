//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//   A class for unfolding 1, 2 or 3 dimensions of data using the
//   Bayesian Unfolding algorithm.
//
// Example (1-D):
//
//   // Training
//   // --------
//   // Create a new object
//   RooUnfoldBayesImpl* mytest= new RooUnfoldBayesImpl("unfold","test 1-D Unfolding");
//   // initialise the generated truth arrays and the reconstructed/measured arrays.
//   // in this example both have 30 bins and are in the range -12.5 to 10.0
//   ntdims = 1 ; ntrue.push_back(30) ; tmin.push_back(-12.5) ; tmax.push_back(10.0);
//   nrdims = 1 ; rtrue.push_back(30) ; rmin.push_back(-12.5) ; rmax.push_back(10.0);
//   mytest.build(ntdims, ntrue, tmin, tmax, nrdims, rtrue, rmin, rmax) ;
//   // Now train. For each event, pass the true generated value and the measured value. If the
//   // the measured value does not exist (because it failed a trigger cut say),
//   // then the measured array should have zero length.
//   for (int ievt=0; ievt < ngenerated; ievt++) {
//      // example event where measured value is accepted
//     xtrue.clear() ; xtrue.push_back(3.56); // generated value
//     xmeas.clear() ; xmeas.push_back(2.75); // measured value
//     mytest.accumulate(xtrue,xmeas);
//      // example event where measured value is rejected
//     xtrue.clear() ; xtrue.push_back(3.56); // generated value
//     xmeas.clear() ;                        // measured value failed cut
//     mytest.accumulate(xtrue,xmeas);
//   }
//   // etc...
//   // You can now save the response matrix if you like.
//
//   // Unfolding
//   // ---------
//   // Accumulate the measured variables you wish to unfold
//   for (int ievt=0; ievt < nmeasured; ievt++) {
//     xmeas.clear() ; xmeas.push_back(4.15); // measured value
//     mytest.accumulate(xmeas); // etc..
//   }
//   //when finished unfold and return the unfolded vector of values
//   vector<Double_t> causes;
//   mytest.unfold(causes);
//
// Authors: Fergus Wilson <fwilson@slac.stanford.edu> and Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldBayesImpl.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "Array2D.h"

using std::vector;
using std::fabs;
using std::sqrt;
using std::min;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::exit;

ClassImp(RooUnfoldBayesImpl);

//-------------------------------------------------------------------------
// default constructor
RooUnfoldBayesImpl::RooUnfoldBayesImpl(): TNamed()
{
  init();
}

//-------------------------------------------------------------------------
RooUnfoldBayesImpl::RooUnfoldBayesImpl(const char* name, const char* title): TNamed(name, title)
{
  init();
}

//-------------------------------------------------------------------------
RooUnfoldBayesImpl::~RooUnfoldBayesImpl()
{
  reset();
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::build(Int_t ndims, const vector<Int_t>& nx, const vector<Double_t>& xmin, const vector<Double_t>& xmax,
              Int_t mdims, const vector<Int_t>& mx, const vector<Double_t>& ymin, const vector<Double_t>& ymax)
{
  // Initialise the unfolding class.
  // ndims = number of truth variables to be unfolded (1 to 3 at the moment).
  // nx   = number of bins for each truth variable
  // xmin = minimum value for each truth variable
  // xmax = maximum value for each truth variable
  // mdims = number of measured variables to be unfolded (probably the same as ndims).
  // ny   = number of bins for each measured variable
  // ymin = minimum value for each measured variable
  // ymax = maximum value for each measured variable

  _ndims = ndims;
  _mdims = mdims;
  if ( (_ndims > 3) || (_ndims < 1)) {
    cerr << "Can only do 3 or less truth dimensions ! Exiting." << endl;
    return(0);
  }

  if ( (_mdims > 3) || (_mdims < 1)) {
    cerr << "Can only do 3 or less measured dimensions ! Exiting." << endl;
    return(0);
  }

  UInt_t nd = ndims;
  if ((nx.size() != nd) || (xmin.size() != nd) || (xmax.size() != nd)) {
    cerr << "Expected truth vector of length " << ndims << "  ! Exiting." << endl;
  }

  UInt_t md = mdims;
  if ((mx.size() != md) || (ymin.size() != md) || (ymax.size() != md)) {
    cerr << "Expected measured vector of length " << mdims << "  ! Exiting." << endl;
  }

  reset();

  _nt.clear();
  _tmin.clear();
  _tmax.clear();
  _nr.clear();
  _rmin.clear();
  _rmax.clear();
  _nOutEffectsRange = _nOutMeasuredRange = _nOutGeneratedRange = 0;
  _nCausesError = _neffects = _ncauses = 0.0;

  _nc = _ne = 1;
  for (Int_t i=0 ; i < _ndims ; i++) {
    // truth dimensions
    _nt.push_back(nx[i]);
    _tmin.push_back(xmin[i]);
    _tmax.push_back(xmax[i]);
    _nc *= nx[i];
  }

  for (Int_t i=_ndims ; i < 3 ; i++) {
    // truth dimensions
    _nt.push_back(0);
  }

  for (Int_t i=0 ; i < _mdims ; i++) {
    // measured dimensions
    _nr.push_back(mx[i]);
    _rmin.push_back(ymin[i]);
    _rmax.push_back(ymax[i]);
    _ne *= mx[i];
  }

  for (Int_t i=_mdims ; i < 3 ; i++) {
    // truth dimensions
    _nr.push_back(0);
  }

  // initialise 2D arrays
  _Sij = new Array2D(_nc,_ne);
  _Nij = new Array2D(_nc,_ne);
  _Mij = new Array2D(_nc,_ne);
  _Vij = new Array2D(_nc,_nc);

  // initialise vectors
  for (Int_t i=0; i < _nc ; i++) {
    _nCi.push_back(0.0);
    _causes.push_back(0.0);
    _causesError.push_back(0.0);
  }
  for (Int_t i=0; i < _ne ; i++) {
    _nEj.push_back(0.0);
    _nEstj.push_back(0.0);
  }

  _initialised_vectors = true;
  _trained = false;

  return(1);
}

//-------------------------------------------------------------------------
TH1*
RooUnfoldBayesImpl::histCauses(const char* name, const char* title) const
{
  // Returns pointer histogram of unfolded causes
  // name = name of histogram
  // title = Histogram title
  //
  Bool_t goterr= (_nCausesError != 0.0);
  if (!goterr) cerr << "Covariance matrix not calculated - fill histogram errors with sqrt(N)" << endl;

  if        (_ndims == 1) {
    TH1D* h = new TH1D(name, title, _nt[0], _tmin[0], _tmax[0]);
    for (Int_t i = 0 ; i < _nc ; i++) {
      h->SetBinContent (i+1, _causes[i]);
      h->SetBinError   (i+1, goterr ? _causesError[i] : sqrt (fabs (_causes[i])));
    }
    return h;
  } else if (_ndims == 2) {
    TH2D* h = new TH2D(name, title,
                       _nt[0], _tmin[0], _tmax[0],
                       _nt[1], _tmin[1], _tmax[1]);
    Int_t ix= 1, iy= 1, nbx= _nt[0];
    for (Int_t i= 0; i < _nc; i++) {
      h->SetBinContent (ix, iy, _causes[i]);
      h->SetBinError   (ix, iy, goterr ? _causesError[i] : sqrt (fabs (_causes[i])));
      if (++ix>nbx) { ix= 1; iy++; }
    }
    return h;
  } else if (_ndims == 3) {
    TH3D* h = new TH3D(name, title,
                       _nt[0], _tmin[0], _tmax[0],
                       _nt[1], _tmin[1], _tmax[1],
                       _nt[2], _tmin[2], _tmax[2]);
    Int_t ix= 1, iy= 1, iz= 1, nbx= _nt[0], nby= _nt[1];
    for (Int_t i= 0; i < _nc; i++) {
      h->SetBinContent (ix, iy, iz, _causes[i]);
      h->SetBinError   (ix, iy, iz, goterr ? _causesError[i] : sqrt (fabs (_causes[i])));
      if (++ix>nbx) {
        ix= 1;
        if (++iy>nby) { iy= 1; iz++; }
      }
    }
    return h;
  } else {
    return 0;
  }
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::accumulate(const vector<Double_t>& xtrue, const vector<Double_t>& xmeas)
{
  // Mapping of measured to generated events for training. Called once per generated event
  // even if the measured value is rejected by efficiency cuts etc...
  //
  // xtrue = vector of ndims length containing the generated truth values
  // xmeas = vector of mdims length containing the corresponding measured values. If the
  //         the measured point is rejected (e.g. by the trigger, failed to reconstruct it)
  //         this vector should have zero length.

  Double_t xm[3], xt[3];
  for (Int_t i=0; i < 3 ; i++) {xm[i]=0;xt[i]=0;}

  Int_t it(-1), jr(-1);

  // measured
  if (xmeas.size() > 0) {
    if ((UInt_t)_mdims != xmeas.size()) {
      cerr << "Unexpected length for reconstructed variables. Expected "
           << _mdims <<"; Found "<< xmeas.size() << endl;
      return(0);
    }

    for (Int_t i=0; i < _mdims ; i++) {xm[i] = xmeas[i];}

    jr = getRecoBin (xm[0], xm[1], xm[2]);
    if (jr==-1) {
      if (debug()) cout << "Outside Measured Range: " << xm[0] << endl;
      _nOutMeasuredRange++;
      return(1);
    }

    _nEj[jr] += 1.0;  // P(E_j) * Nobs
  }

  // truth
  if (xtrue.size() > 0) {
    if ((UInt_t)_ndims != xtrue.size()) {
      cerr << "Unexpected length for generated variables. Expected "
           << _ndims <<"; Found "<< xtrue.size() << endl;
      return(0);
    }

    for (Int_t i=0; i < _ndims ; i++) { xt[i] = xtrue[i];}

    it = getTruthBin(xt[0], xt[1], xt[2]);
    if (it==-1) {
      if (debug()) cout << "Outside Generated Range: " << xt[0] << endl;
      _nOutGeneratedRange++;
      return(1);
    }

    _nCi[it] += 1.0;  // P_0(C_i) * Ntrue
  }

  //Mapping of Causes to Effects
  if (debug()) cout << _Nij->GetNrows() << " " <<_Nij->GetNcols() << endl;
  if ((it!=-1) && (jr!=-1)) { _Nij->Add(it,jr,1.0);}

  return(1);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::accumulate(const vector<Double_t>& xmeas)
{
  // Accumulate the measured events you want to unfold.
  // xmeas = vector of mdims length containing the measured values.

  if ((UInt_t)_mdims != xmeas.size()) {
    cerr << "Unexpected length for reconstructed variables. Expected "
         << _mdims <<"; Found "<< xmeas.size() << endl;
    return(0);
  }

  Double_t xm[3];
  for (Int_t i=0; i < 3 ; i++) {xm[i]=0;}

  for (Int_t i=0; i < _mdims ; i++) { xm[i] = xmeas[i]; }

  Int_t jr = getRecoBin (xm[0], xm[1], xm[2]);
  if (jr==-1) {
    if (warn()) cout << "Outside Measured Range: " << xm[0] << " " << xm[1] << endl;
    _nOutEffectsRange++;
    return(1);
  }

  // measured
  _nEstj[jr] += 1.0;
  _neffects++;
  return(1);
}

Int_t
RooUnfoldBayesImpl::setupTrain(const vector<Double_t>& nCi, const vector<Double_t>& nEj, const Array2D& Nij)
{
  if ((Int_t)nCi.size() != _nc || (Int_t)nEj.size() != _ne || Nij.GetNrows() != _nc || Nij.GetNcols() != _ne) {
    cerr << "setupTrain: wrong vector size" << endl;
    return(0);
  }
  _nCi= nCi;
  _nEj= nEj;
  *_Nij= Nij;
  return(1);
}

Int_t RooUnfoldBayesImpl::setupUnfold(const vector<Double_t>& nEstj)
{
  if ((Int_t)nEstj.size() != _ne) {
    cerr << "setupUnfold: wrong vector size" << endl;
    return(0);
  }
  _nEstj= nEstj;
  for (Int_t i= 0; i < _ne; i++)
    _neffects += _nEstj[i];
  return(1);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getTruthBin(Double_t x, Double_t y, Double_t z) const
{
  // return the index of the bin corresponding to truth value (x,y,z)
  // where y,z may be optional.
  Int_t i(0), j(0), k(0);
  if (_ndims > 0) { i = getBin(x, _tmin[0], _tmax[0], _nt[0]);}
  if (_ndims > 1) { j = getBin(y, _tmin[1], _tmax[1], _nt[1]);}
  if (_ndims > 2) { k = getBin(z, _tmin[2], _tmax[2], _nt[2]);}

  Int_t ioffset = k*(_nt[0]*_nt[1]) + j*_nt[0] + i;


  if (ioffset<0||ioffset>_nc) { ioffset = -1; }

  if (ioffset==-1) {
    if (debug()) cout << "gTB : nt " << _nt[0] << " "<< _nt[1] << " "<< _nt[2] << endl
                      << "gTB : " << x << " "<< y << " "<< z << " "
                      << i << " " << j << " "<< k << " " << ioffset << endl;
  }
  return(ioffset);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getRecoBin(Double_t x, Double_t y, Double_t z) const
{
  // Return the index of the bin corresponding to measured value (x,y,z)
  // where y,z may be optional.
  Int_t i(0), j(0), k(0);
  if (_mdims > 0) { i = getBin(x, _rmin[0], _rmax[0], _nr[0]);}
  if (_mdims > 1) { j = getBin(y, _rmin[1], _rmax[1], _nr[1]);}
  if (_mdims > 2) { k = getBin(z, _rmin[2], _rmax[2], _nr[2]);}

  Int_t ioffset = k*(_nr[0]*_nr[1]) + j*_nr[0] + i;

  if (ioffset<0||ioffset>_ne) { ioffset = -1; }

  return(ioffset);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getBin(Double_t x, Double_t xmin, Double_t xmax, Int_t nx) const
{
  // Return the index to the element corresponding to value x

  if (debug()) cout << x << " " << xmin << " " << xmax << " " << nx << " " << endl;
  Double_t temp = (x - xmin) / (xmax - xmin);
  temp *= nx;
  Int_t bin = (Int_t) floor(temp);

  return (bin);
}

//-------------------------------------------------------------------------
void
RooUnfoldBayesImpl::init()
{
  _Sij= _Nij= _Mij= _Vij= 0;
  _initialised_vectors = false;
  _trained = false;
  _fDebug = 1;
  _UnfoldType = Iterative;
}

//-------------------------------------------------------------------------
void
RooUnfoldBayesImpl::reset()
{
  // Delete the internal matrices
  delete _Sij;
  delete _Nij;
  delete _Mij;
  delete _Vij;
  init();
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::clear()
{
  // Reset the unfolding. clears out training and measured samples
  vector<Int_t> nx(_nt), ny(_nr);
  nx.resize(_ndims);
  ny.resize(_mdims);
  // Use copies of vectors because _tmin etc will be modified inside build
  return build (_ndims, nx, vector<Double_t>(_tmin), vector<Double_t>(_tmax),
                _mdims, ny, vector<Double_t>(_rmin), vector<Double_t>(_rmax));
}


//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::train(Int_t iterations, Bool_t smoothit)
{
  // After accumulating the training sample, calculate the unfolding matrix.
  // iterations = number of iterations to perform (3 by default).
  // smoothit = smooth the matrix in between iterations (false, not implemented yet).

  if (!proceed()) {return(0);}

  // zero efficiency
  _efficiencyCi.clear();
  for (Int_t i=0 ; i< _nc ; i++) {
    _efficiencyCi.push_back(0);
  }
  Double_t ntrue = sum(_nCi);
  if (debug()) cout << " efficiency "<< _efficiencyCi.size() << " " << ntrue << endl;

  // Initial distribution
  vector<Double_t> P0C(_nCi);
  for (UInt_t i = 0 ; i <_nCi.size(); i++) {
    P0C[i] /= ntrue;
    //P0C[i] = 1.0/_nCi.size();
    if (debug()) cout << "i _nCi P0C " << i << "\t" << _nCi[i] << "\t" << P0C[i] << endl;
  }

  for (Int_t kiter = 0 ; kiter < iterations; kiter++) {

    if (verbose()) cout << "Iteration : " << kiter << endl;
    // Smearing matrix S
    for (Int_t i = 0 ; i < _nc ; i++) {
      _efficiencyCi[i] = 0.0;
      if (debug()) cout << "smearing i " << i << endl;
      for (Int_t j = 0 ; j < _ne ; j++) {
        _Sij->Set(i,j,0.0);
        if (_nCi[i] <= 0.0) {if (debug()) cout << "_nCi[i="<<i<<"] is " << _nCi[i] << endl; continue;}
        if (_Nij->Get(i,j) <= 0.0) {if (debug()) cout << "_Nij("<<i<<","<<j<<") is " << _Nij->Get(i,j) << endl; continue;}
        if (debug()) cout << "  smearing j " << j << endl;

        Double_t PEjCi = _Nij->Get(i,j) / _nCi[i];
        if (debug()) cout << " _Nij("<<i<<","<<j<<") "<< _Nij->Get(i,j) << "PEjCi " << PEjCi << " " << _nEj[j] << " "<< _nCi[i] << endl;
        // efficiency of detecting the cause Ci in Effect Ej
        _efficiencyCi[i] +=  PEjCi ;

        Double_t numer = PEjCi * P0C[i];
        Double_t denom(0.0);
        for (Int_t l = 0 ; l < _nc ; l++) {
          if (_nCi[l] <= 0.0) {if (debug()) cout << "_nCi[l="<<l<<"] is " << _nCi[l] << endl; continue;}
          Double_t PEjCl = _Nij->Get(l,j) / _nCi[l];
          if (debug()) cout << "l " << l << " " << PEjCl << " " << P0C[l] << endl;
          denom += (PEjCl * P0C[l]) ;
        }
        if (denom <= 0.0) {if (debug()) cout << "would set Sij("<<i<<","<<j<<") to " << numer << "/" << denom << endl; continue;}
        _Sij->Set(i,j,numer / denom) ;
        if (debug()) cout << "Sij " << _Sij->Get(i,j) << endl;
      }
    }
    if (debug()) {
      for (Int_t i=0 ; i< _nc ; i++) {
        cout << "efficiency i " << i << " :  " << _efficiencyCi[i] << endl;;
      }
    }

    vector<Double_t> nbarCi(_nc);

    _nbartrue = 0.0;

    // best estimate of true number of events
    for (Int_t i = 0 ; i < _nc ; i++) {
      nbarCi[i] = 0.0;
      for (Int_t j = 0 ; j < _ne ; j++) {
        _Mij->Set(i,j,0.0);   // unfolding matrix
        if (_efficiencyCi[i] <= 0.0) {continue;}
        Double_t x = _Sij->Get(i,j) / _efficiencyCi[i];
        _Mij->Set(i,j,x);
        nbarCi[i] += (_nEstj[j] * x);
        if (debug() && _Mij->Get(i,j)>0) {cout << "Mij " << i << " " << j << " " << _Mij->Get(i,j) << endl;}
      }
      if (verbose()>=2) cout << "nbarCi " << i << "\t" << nbarCi[i] << endl;
      _nbartrue += nbarCi[i];
    }

    // new estimate of true distribution
    vector<Double_t> PbarCi(nbarCi);
    if (verbose()>=2) cout << "nbartrue " << _nbartrue << endl;
    for (UInt_t i = 0 ; i < nbarCi.size(); i++) {
      PbarCi[i] /= _nbartrue;
      if (verbose()>=2) cout << "i PbarCi P0C " << i
                          << "\t" << PbarCi[i]*_nbartrue
                          << "\t" << P0C[i]*_nbartrue << endl;
    }

    // chi^2 (but need to calculate error)
    // takes a verrrrrrrry long time
    // only do when unfolded?
    //if (_nc*_ne < 50000) {getCovariance();}

    // so not smooth the last iteraction
    if (kiter < (iterations-1)) {
      if (smoothit) {smooth(PbarCi);}
      //smooth(PbarCi);
    }

    Double_t chi2 = getChi2(PbarCi, P0C, _nbartrue);
    if (verbose()) cout << "Chi^2 of change " << chi2 << endl;

    // replace P0C
    P0C = PbarCi;

    // and repeat
  }

  _trained = true;

  return(1);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::smooth(vector<Double_t>& PbarCi, Double_t nevts) const
{
  // Smooth unfolding distribution. PbarCi is the array of proababilities
  // to be smoothed PbarCi; nevts is the numbers of events
  // (needed to calculate suitable errors for the smearing).
  // PbarCi is returned with the smoothed distribution.

  if (_ndims != 1) {
    cerr << "Smoothing only implemented for 1-D distributions" << endl;
    return (0);
  } else {
    if (verbose()) cout << "Smoothing." << endl;
    TH1D* h = new TH1D("hsmooth", "Smoothed Causes", _nt[0], _tmin[0], _tmax[0]);
    // pack histogram
    for (Int_t i = 0 ; i < _nc ; i++) {
      h->SetBinContent(i+1, PbarCi[i]);
      // calculate error on probability
      Double_t error = sqrt (fabs (PbarCi[i]/nevts));
      h->SetBinError(i+1, error);
    }
    h->Smooth(); // smooth the histogram
    // unpack histogram into array
    for (Int_t i = 0 ; i < _nc ; i++) {
      PbarCi[i] = h->GetBinContent(i+1);
    }
    delete h;
  }

  return(1);
}

//-------------------------------------------------------------------------
Double_t
RooUnfoldBayesImpl::getChi2(const vector<Double_t> prob1,
                            const vector<Double_t> prob2,
                            Double_t nevents) const
{
  // calculate the chi^2. prob1 and prob2 are the probabilities
  // and nevents is the number of events used to calculate the probabilities
  Double_t chi2(0);
  if (debug()) cout << "chi2 " << prob1.size() << " " << prob2.size() << " " << nevents << endl;
  for (UInt_t i = 0 ; i < prob1.size() ; i++) {
    Double_t psum  = (prob1[i] + prob2[i])*nevents;
    Double_t pdiff = (prob1[i] - prob2[i])*nevents;
    if (psum > 1) {
      chi2 = chi2 + (pdiff*pdiff)/psum;
    } else {
      chi2 = chi2 + (pdiff*pdiff);
    }
  }
  return(chi2);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::info(Int_t level) const
{
  // Print out some useful info of progress so far

  cout << "-------------------------------------------" << endl;
  cout << "Unfolding Algorithm" << endl;
  cout << "Generated (Training):" << endl;
  cout << "  Number of Generated Dimensions : " << _ndims << endl;
  for (Int_t i=0 ; i < _ndims; i++) {
    cout << "  Dimension " << i+1 << ": Range " << _tmin[i] << " to "
         << _tmax[i] << " in " << _nt[i] << " bins." << endl;
  }
  cout << "  Total Number of bins   : " << _nc << endl;
  Double_t ntrue = sum(_nCi);
  cout << "  Total Number of events : " << ntrue << endl;
  cout << "  Events outside Range   : " << _nOutGeneratedRange << endl;

  cout << "Measured (Training):" << endl;
  cout << "  Number of Reconstructed Dimensions : " << _mdims << endl;
  for (Int_t i=0 ; i < _mdims; i++) {
    cout << "  Dimension " << i+1 << ": Range " << _rmin[i] << " to "
         << _rmax[i] << " in " << _nr[i] << " bins." << endl;
  }
  cout << "  Total Number of bins   : " << _ne << endl;
  Double_t nobs = sum(_nEj);
  cout << "  Total Number of events : " << nobs << endl;
  cout << "  Events outside Range   : " << _nOutMeasuredRange << endl;

  cout << "Input (for unfolding):" << endl;
  cout << "  Total Number of events : " << _neffects << endl;
  cout << "  Events outside Range   : " << _nOutEffectsRange << endl;

  cout << "Output (unfolded):" << endl;
  cout << "  Total Number of events : " << _ncauses
       << " +/- " << _nCausesError <<endl;

  if (level == 0) {}
  if (level == 1) {}
  if (level == 2) {}

  cout << "-------------------------------------------\n" << endl;

  if ((sum(_nEj)!=0) || (sum(_nCi)!=0)) {
    Int_t iend = min(_nCi.size(),_nEj.size());
    if (_ndims==2) {
      cout << "    \t     \t     \t Train \t Train\t Test\t Unfolded"<< endl;
      cout << "Bin \tBin x\tBin y\t Truth \t Reco\t Input\t Output"<< endl;
    } else {
      cout << "    \tTrain \tTrain\tTest\tUnfolded"<< endl;
      cout << "Bin \tTruth \tReco\tInput\tOutput"<< endl;
    }
    Int_t ir(0), ic(0);
    for (Int_t i=0; i < iend ; i++) {
      ic = i / _nr[0];
      ir = i - (ic*_nr[0]);
      if ((_nCi[i] == 0) && (_nEj[i] == 0) &&
          (_nEstj[i] == 0) && (_causes[i]==0)) {continue;}
      if (_ndims==2) {
        cout << i << "\t " << ir << "\t " << ic << "\t" << _nCi[i] \
             << "\t " << _nEj[i] << "\t " << _nEstj[i] << "\t " << _causes[i] << endl;
      } else {
        cout << i << "\t" << _nCi[i] \
             << "\t " << _nEj[i] << "\t " << _nEstj[i] << "\t " << _causes[i] << endl;
      }
    }

    // if the number of bins is different
    if (_nCi.size() > _nEj.size() ) {
      for (UInt_t i=iend; i < _nCi.size() ; i++) {
        cout << i << "\t " << _nCi[i] << endl;
      }
    } else {
      for (UInt_t i=iend; i < _nEj.size() ; i++) {
        cout << i << "\t \t " << _nEj[i] << endl;
      }
    }

    cout << "--------------------------------------------------------" << endl;
    if (_ndims==2) {
      cout << " \t\t\t" << sum(_nCi) << "\t " << sum(_nEj)
           << "\t " << sum(_nEstj) << "\t " << sum(_causes) << endl;
    } else {
      cout << " \t" << sum(_nCi) << "\t " << sum(_nEj)
           << "\t " << sum(_nEstj) << "\t " << sum(_causes) << endl;
    }
    cout << "--------------------------------------------------------\n" << endl;
  }


  return(1);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::trainBinByBin(Bool_t smoothit)
{
  // After accumulating the training sample, calculate the unfolding matrix.

  if (!proceed()) {return(0);}

  // calculate efficiency
  _efficiencyCi.clear();
  for (Int_t i=0 ; i< _nc ; i++) {
    Double_t eff(0.0);
    if (_nCi[i]>0) { eff = _nEj[i] / _nCi[i];}
    _efficiencyCi.push_back(eff);
  }

  //Double_t ngen  = sum(_nCi);
  //Double_t nmeas = sum(_nEj);

  // Unfolding matrix (inverse of efficiency)
  for (Int_t i = 0 ; i < _nc ; i++) {
    for (Int_t j = 0 ; j < _ne ; j++) {
      _Mij->Set(i,j,0.0);
      if ((i==j) && (_efficiencyCi[i]>0)) {_Mij->Set(i,j, 1.0/_efficiencyCi[i]);}
    }
  }

  _trained = true;

  return(1);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getCovarianceBinByBin()
{
  const vector<Double_t>& effects= _nEstj;
  // Calculate error

  // Create the covariance matrix for Bin by Bin

  // Unfolded = effects * Truth / Generated
  // => variance = sum in quadrature of errors

  for (Int_t k = 0 ; k < _nc ; k++) {
    for (Int_t l = 0 ; l < _nc ; l++) {
      Double_t variance(0);
      if ((k==l) && (_nEj[k]>0) && (_nCi[k]>0) && (effects[k]>0)) {
        Double_t f = effects[k] * _nCi[k] / _nEj[k];
        variance = (1.0/effects[k] + 1.0/_nCi[k] + 1.0/_nEj[k]) * (f*f);
      }
      _Vij->Set(k, l, variance);
    }
  }
  _nCausesError = getError();
  return(1);
}
//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getVariance()
{
    const vector<Double_t>& effects= _nEstj;
  vector<Double_t> dummy;
  Double_t nbartrue = getnbarCi(effects,dummy);
  for (Int_t k = 0 ; k < _nc ; k++) {
      Double_t temp=0; 
      Double_t temp2=0;
      for (Int_t j = 0 ; j < _ne ; j++) {
        Double_t Mkj = _Mij->Get(k,j);
        Double_t ratio = effects[j]/nbartrue;
        temp += (Mkj*Mkj*effects[j]*(1-ratio));

        for (Int_t i = 0 ; i <j-1 ; i++) {
          Double_t Mki = _Mij->Get(k,i);
          temp2 += 2*(Mki*Mkj*effects[i]*ratio);
        }
      }

      _Vij->Set(k,k,(temp-temp2));
  }
  _nCausesError = getError();
    return 1;   
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getCovariance(Bool_t doUnfoldSystematic)
{
  const vector<Double_t>& effects= _nEstj;
  // Create the covariance matrix
  if (warn() && _nc*_ne >= 50625)
    cout << "getCovariance (this takes some time with " << _nc << " x " << _ne << " bins)." << endl;
  vector<Double_t> dummy;
  Double_t nbartrue = getnbarCi(effects,dummy);

  // error  error from data
  if (verbose()) cout << "Calculating covariance due to number of measured events..." << endl;
  Array2D *Vc0 = new Array2D(_nc,_nc);
  for (Int_t k = 0 ; k < _nc ; k++) {
    for (Int_t l = k ; l < _nc ; l++) {
      Vc0->Set(k,l,0.0);
      Double_t temp(0), temp2(0);
      for (Int_t j = 0 ; j < _ne ; j++) {
        Double_t Mlj = _Mij->Get(l,j);
        if (Mlj == 0) {continue;}  // skip zero elements
        Double_t Mkj = _Mij->Get(k,j);

        Double_t ratio = effects[j]/nbartrue;
        temp += (Mkj*Mlj*effects[j]*(1-ratio));

        for (Int_t i = 0 ; i < _ne ; i++) {
          if (i==j) {continue;}
          Double_t Mki = _Mij->Get(k,i);
          temp2 += (Mki*Mlj*effects[i]*ratio);
        }
      }

      Vc0->Set(k,l,(temp-temp2));
      Vc0->Set(l,k,(temp-temp2)); // symmetric matrix
    }
  }

  // error due to uncertainty on unfolding matrix M
  // This is disabled by default: I'm not sure it is correct, it is very slow, and
  // the effect should be small with good MC statistics.
  Array2D *Vc1 = new Array2D(_nc,_nc);  // automatically zeroed
  if (doUnfoldSystematic) {
    if (verbose()) cout << "Calculating covariance due to unfolding matrix..." << endl;

    // Pre-compute some numbers
    vector<Double_t> inv_nCi(_nCi);
    Array2D *inv_npec = new Array2D(_nc,_ne);  // automatically zeroed
    for (Int_t k = 0 ; k < _nc ; k++) {
      if (inv_nCi[k] != 0) {inv_nCi[k] = 1.0 / inv_nCi[k];}
      for (Int_t i = 0 ; i < _ne ; i++) {
        if (inv_nCi[k] == 0) {continue;}
        Double_t pec  = _Nij->Get(k,i) / _nCi[k];
        Double_t temp = inv_nCi[k] / pec;
        if (pec !=0) {inv_npec->Set(k,i,temp); }
      }
    }
    //
    Array2D *M_tmp = new Array2D(_nc,_ne);  // automatically zeroed
    for (Int_t i = 0 ; i < _ne ; i++) {
      Double_t temp(0);
      // diagonal element
      for (Int_t u = 0 ; u < _nc ; u++) {
        temp = _Mij->Get(i,u) * _Mij->Get(i,u) * _efficiencyCi[u] * _efficiencyCi[u]
          * (inv_npec->Get(u,i) - inv_nCi[u]);
        M_tmp->Add(i, i, temp);
      }

      // off-diagonal element
      for (Int_t j = i+1 ; j < _ne ; j++) {
        for (Int_t u = 0 ; u < _nc ; u++) {
          temp = -1.0 * _Mij->Get(i,u) * _Mij->Get(j,u) * _efficiencyCi[u] * _efficiencyCi[u]
            * inv_nCi[u];
          M_tmp->Add(j, i, temp);
        }
        M_tmp->Set(i, j, M_tmp->Get(j, i)); // symmetric matrix
      }
    }

    // now calculate covariance
    Double_t neff_inv(0);
    for (Int_t k = 0 ; k < _nc ; k++) {
      (_efficiencyCi[k] != 0) ? neff_inv = inv_nCi[k]  / _efficiencyCi[k] : neff_inv = 0;
      for (Int_t l = k ; l < _nc ; l++) {
        for (Int_t i = 0 ; i < _ne ; i++) {
          for (Int_t j = 0 ; j < _ne ; j++) {
            Double_t covM = _Mij->Get(i,l) * inv_nCi[l] +
                            _Mij->Get(j,k) * inv_nCi[k] + M_tmp->Get(j,i);
            if (k==l) {
              covM -= neff_inv;
              if (i==j) {covM += inv_npec->Get(k,i);}
            }
            if (i==j) {
              covM -=  (_Mij->Get(i,l) * _efficiencyCi[l] * inv_npec->Get(l,i) +
                        _Mij->Get(i,k) * _efficiencyCi[k] * inv_npec->Get(k,i) );
            }
            covM +=  _Mij->Get(i,k) * _Mij->Get(j,l);

            Double_t temp = effects[i] * effects[j] * covM;
            Vc1->Add(l, k, temp);
          } // j...
        } // i...
        Vc1->Set(k, l, Vc1->Get(l,k));
      } // l...
    } // k...
    delete inv_npec;
  } // if (doUnfodsystematic

  // to get complete covariance add together
  Double_t nbar2 = nbartrue*nbartrue;
  for (Int_t k = 0 ; k < _nc ; k++) {
    for (Int_t l = k ; l < _nc ; l++) {
      Double_t temp = Vc0->Get(k,l);
      if (doUnfoldSystematic) temp += Vc1->Get(k,l) / nbar2;  // divide by nbartrue*nbartrue to get probability covariance matrix
      _Vij->Set(k,l,temp);
      _Vij->Set(l,k,temp);
    }
  }

  if (debug(3)) {
    cout << " k   l  Vc0(k,l)   Vc1(k,l)  _Vij(k,l)" << endl;
    for (Int_t k = 0 ; k < _nc ; k++) {
      for (Int_t l = 0 ; l < _nc ; l++) {
        if (_Vij->Get(k,l) == 0) {continue;}
        cout << left << setw(4) << k << setw(4) << l
             << setw(11) << Vc0->Get(k,l) << setw(10) << Vc1->Get(k,l)
             << setw(10)<< _Vij->Get(k,l) << right << endl;
      }
    }
  }
  delete Vc0;

  _nCausesError = getError();
  return(1);
}

//-------------------------------------------------------------------------
Double_t
RooUnfoldBayesImpl::deltaM(Int_t k, Int_t i, Int_t r, Int_t u) const
{
  // Helper function for getCovariance. Calculates derivative
  // d M_{ki} / d P(E_r | C_u)
  //
  Double_t temp(0);

  if (_Mij->Get(k,i) == 0.0) {return(temp);}

  if ((k==u) && (r==i)) {
    Double_t PErCu = _Nij->Get(u,r)  / _nCi[u];
    temp = 1.0 / PErCu;
  }
  if (k==u) {
    temp -= (1.0 / _efficiencyCi[u]);
  }
  if (r==i) {
    Double_t PEiCu = _Nij->Get(u,i) / _nCi[u];
    temp -= (_Mij->Get(u,i) * _efficiencyCi[u] / PEiCu);
  }

  temp = temp * _Mij->Get(k,i);
  return (temp);
}

//-------------------------------------------------------------------------
Bool_t
RooUnfoldBayesImpl::proceed() const
{
  // Check we are in a fit state to calculate unfolding matrix

  Bool_t statusOK(true);

  if (!_initialised_vectors) {
    cerr << "build not properly called yet" << endl;
    statusOK = false;
  }

  if (!statusOK) {
    cerr << "Exiting" << endl;
    exit(1);
  }
  return (statusOK);

}

//----------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::setDebug(Int_t Db) {
  // Initialises copious debugging printout
  _fDebug = Db;
  return(1);
}

//----------------------------------------------------------------------
Double_t
RooUnfoldBayesImpl::sum(const vector<Double_t>& x) const
{
  // sum of array x
  Double_t total(0);
  for (UInt_t i = 0 ; i < x.size() ; i++) { total += x[i]; }
  return(total);
}

//----------------------------------------------------------------------
Double_t
RooUnfoldBayesImpl::getnbarCi(const vector<Double_t>& effects,
                        vector<Double_t>& causes) const
{
  // Unfold the measured effects and return the estimated true causes

  causes.clear();

  // unfold
  Double_t nbartrue(0);
  for (Int_t i = 0 ; i < _nc ; i++) {
    Double_t temp(0);
    for (Int_t j = 0 ; j < _ne ; j++) {
      temp += (_Mij->Get(i,j) * effects[j]);
    }
    causes.push_back(temp);
    nbartrue += temp;
  }

  if (debug()) cout << "getnbarCi nbartrue " << nbartrue << endl;
  return(nbartrue);
}

//----------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getBin(Int_t index, Bool_t truth, vector<Int_t>& coords) const
{
  // Return the bin index
  Int_t irow(-1), icol(-1), iheight(-1);
  Int_t nd(-1), nrows(-1), ncols(-1);
  if (truth) {
    nd = _mdims;
    if (nd>1) {nrows = _nt[0];}
    if (nd>2) {ncols = _nt[1];}
  } else {
    nd = _ndims;
    nrows  = _nr[0];
    if (nd>1) {nrows = _nr[1];}
    if (nd>2) {ncols = _nt[2];}
  }

  if (nd==1) {
    irow = index;
  } else if (nd==2) {
    icol = Int_t (index / nrows);
    irow = index - icol * nrows;
  } else if (nd==3) {
    // this is wrong
    iheight = Int_t (index / (nrows*ncols));
    icol = Int_t (index / nrows);
    irow = index - icol * nrows;
  }
  coords.clear();
  coords.push_back(irow);
  coords.push_back(icol);
  coords.push_back(iheight);

  return(1);
}

//----------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::unfold(vector<Double_t>& causes, Int_t iterations, Bool_t smoothit)
{
  // Unfold the accumulated measured values and return the true values in a matrix.
  // causes = vector of unfolded values.

  _UnfoldType = Iterative;
  if (verbose()) cout << "Now unfolding..." << endl;
  // unfold
  if (!train(iterations, smoothit)) return(0);
  _ncauses = getnbarCi(_nEstj, causes);
  _causes.clear();
  _causes = causes;

  return(1);
}

//----------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::unfoldBinByBin(vector<Double_t>& causes, Bool_t smoothit)
{
  // Unfold the accumulated measured values and return the true values in a matrix.
  // causes = vector of unfolded values.

  _UnfoldType = BinByBin;
  if (verbose()) cout << "Now unfolding (bin-by-bin)..." << endl;
  // unfold
  if (!trainBinByBin(smoothit)) return(0);
  _ncauses = getnbarCi(_nEstj, causes);
  _causes.clear();
  _causes = causes;

  return(1);
}

//----------------------------------------------------------------------
Double_t
RooUnfoldBayesImpl::getError()
{
  // Return the total statistical error on the estimated number of true events
  Double_t toterror(0);
  _causesError.clear();
  for (Int_t i = 0 ; i < _nc ; i++) {
    Double_t error = 0;
    for (Int_t j = 0 ; j < _nc ; j++) {
      error += _Vij->Get(i,j);
    }
    toterror += error;

    if (error>0) {error = sqrt(error);}
    if (debug()) cout << "bin " << i << " error " << error << endl;
    _causesError.push_back(error); // update array of errors
  }

  // total error
  if (toterror>0) {toterror = sqrt(toterror);}
  if (debug()) cout << "error " << toterror << endl;

  return(toterror);
}
