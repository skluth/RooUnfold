//-----------------------------------------------------------------
//
// File and Version Information:
//   $Id: RooUnfoldBayesImpl.cxx,v 1.1.1.1 2007-04-04 21:27:25 adye Exp $
//
// Description:
//   Bayesian Unfolding class 
//
// Author List:
//   F.Wilson, T.Adye
//
// Copyright Information:
//   Copyright (C) 2005-2006 Rutherford Appleton Laboratory
//
// Author: Fergus Wilson <mailto:fwilson@slac.stanford.edu>, Tim Adye <mailto:T.J.Adye@rl.ac.uk>
// @(#)
/*
 * Copyright (C) 2005-2006 Rutherford Appleton Laboratory
 */
//------------------------------------------------------------------
#include "RooUnfoldBayesImpl.h"

#include <iostream>
#include <vector>
#include <math.h>

#include "TStopwatch.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "Array2D.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

ClassImp(RooUnfoldBayesImpl);

//-------------------------------------------------------------------------
// default constructor
RooUnfoldBayesImpl::RooUnfoldBayesImpl(): TNamed()
{
  // A class for unfolding 1, 2 or 3 dimensions of data using the Bayesian Unfolding algorithm.
  // Example (1-D):
  //
  // Training
  // --------
  // // Create a new object
  // RooUnfoldBayesImpl* mytest= new RooUnfoldBayesImpl("unfold","test 1-D Unfolding");
  // // initialise the generated truth arrays and the reconstructed/measured arrays.
  // // in this example both have 30 bins and are in the range -12.5 to 10.0
  // ntdims = 1 ; ntrue.push_back(30) ; tmin.push_back(-12.5) ; tmax.push_back(10.0);
  // nrdims = 1 ; rtrue.push_back(30) ; rmin.push_back(-12.5) ; rmax.push_back(10.0);
  // mytest.build(ntdims, ntrue, tmin, tmax, nrdims, rtrue, rmin, rmax) ;
  // // Now train. For each event, pass the true generated value and the measured value. If the
  // // the measured value does not exist (because it failed a trigger cut say), 
  // // then the measured array should have zero length.
  // for (int ievt=0; ievt < ngenerated; ievt++) {
  //    // example event where measured value is accepted
  //   xtrue.clear() ; xtrue.push_back(3.56); // generated value
  //   xmeas.clear() ; xmeas.push_back(2.75); // measured value
  //   mytest.accumulate(xtrue,xmeas);
  //    // example event where measured value is rejected
  //   xtrue.clear() ; xtrue.push_back(3.56); // generated value
  //   xmeas.clear() ;                        // measured value failed cut
  //   mytest.accumulate(xtrue,xmeas);
  // }
  // // etc...
  // // After accumulating all the training events, you now calculated the unfolding matrix
  // mytest.train() ; 
  // // You can now save the unfolding matirx if you like.
  //
  // Unfolding
  // ---------
  // // Accumulate the measured variables you wish to unfold
  // for (int ievt=0; ievt < nmeasured; ievt++) {
  //   xmeas.clear() ; xmeas.push_back(4.15); // measured value
  //   mytest.accumulate(xmeas); // etc..
  // }
  // //when finished unfold and return the unfolded vector of values
  // vector<Double_t> causes;
  // mytest.unfold(causes);
  //

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
      //cout << "Outside Measured Range: " << xm[0] << endl;
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
      //cout << "Outside Generated Range: " << xt[0] << endl;
      _nOutGeneratedRange++;
      return(1);
    }

    _nCi[it] += 1.0;  // P_0(C_i) * Ntrue
  }

  //Mapping of Causes to Effects
  //cout << _Nij->GetNrows() << " " <<_Nij->GetNcols() << endl;
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
    cout << "Outside Measured Range: " << xm[0] << " " << xm[1] << endl;
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
  if (nCi.size() != _nc || nEj.size() != _ne || Nij.GetNrows() != _nc || Nij.GetNcols() != _ne) {
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
  if (nEstj.size() != _ne) {
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
    //cout << "gTB : nt " << _nt[0] << " "<< _nt[1] << " "<< _nt[2] << endl;
    //cout << "gTB : " << x << " "<< y << " "<< z << " "
    //   << i << " " << j << " "<< k << " " << ioffset << endl;
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

  //cout << x << " " << xmin << " " << xmax << " " << nx << " " << endl;
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
  Double_t nobs = sum(_nEj);
  //cout << " efficiency "<< _efficiencyCi.size() << " " << nobs << endl;

  // Initial distribution
  vector<Double_t> P0C(_nCi);
  for (UInt_t i = 0 ; i <_nCi.size(); i++) {
    P0C[i] /= nobs;
    //P0C[i] = 1.0/_nCi.size();
    //cout << "i _nCi P0C " << i << "\t" << _nCi[i] << "\t" << P0C[i] << endl;
  }

  for (Int_t kiter = 0 ; kiter < iterations; kiter++) {

    cout << "Iteration : " << kiter << endl;
    // Smearing matrix S
    for (Int_t i = 0 ; i < _nc ; i++) {
      _efficiencyCi[i] = 0.0;
      //cout << "smearing i " << i << endl;
      for (Int_t j = 0 ; j < _ne ; j++) {
        _Sij->Set(i,j,0.0);
        if (_nCi[i] <= 0.0) {continue;}
        if (_Nij->Get(i,j) <= 0.0) {continue;}
        //cout << "  smearing j " << j << endl;

        Double_t PEjCi = _Nij->Get(i,j) / _nCi[i];
        //cout << "PEjCi " << PEjCi << " _Nij "<< _Nij->Get(i,j) << " " << _nEj[j] << " "<< _nCi[i] << endl;
        // efficiency of detecting the cause Ci in Effect Ej
        _efficiencyCi[i] +=  PEjCi ;

        Double_t numer = PEjCi * P0C[i];
        Double_t denom(0.0);
        for (Int_t l = 0 ; l < _nc ; l++) {
          if (_nCi[l] <= 0.0) {continue;}
          Double_t PEjCl = _Nij->Get(l,j) / _nCi[l];
          //cout << "l " << l << " " << PEjCl << " " << P0C[l] << endl;
          denom += (PEjCl * P0C[l]) ;
        }
        _Sij->Set(i,j,numer / denom) ;
        //cout << "Sij " << _Sij->Get(i,j) << endl;
      }
    }
    for (Int_t i=0 ; i< _nc ; i++) {
      //cout << "efficiency i " << i << " :  " << _efficiencyCi[i] << endl;;
    }

    vector<Double_t> nbarCi(_nCi);

    _nbartrue = 0.0;

    // best estimate of true number of events
    for (Int_t i = 0 ; i < _nc ; i++) {
      nbarCi[i] = 0.0;
      for (Int_t j = 0 ; j < _ne ; j++) {
        _Mij->Set(i,j,0.0);   // unfolding matrix
        if (_efficiencyCi[i] <= 0.0) {continue;}
        Double_t x = _Sij->Get(i,j) / _efficiencyCi[i];
        _Mij->Set(i,j,x);
        nbarCi[i] += (_nEj[j] * _Mij->Get(i,j));
        //      if (_Mij[i][j]>0) {cout << "Mij " << i << " " << j << " " << _Mij[i][j] << endl;}
      }
      //      _nbartrue += nbarCi[i];
    }

    _nbartrue = sum(nbarCi);

    //cout << "nbartrue " << _nbartrue << endl;

    // new estimate of true distribution
    vector<Double_t> PbarCi(nbarCi);
    Double_t diff(0.);
    for (UInt_t i = 0 ; i < nbarCi.size(); i++) {
      PbarCi[i] /= _nbartrue;
      //Double_t delta = (PbarCi[i] - P0C[i])*_nbartrue;
      //cout << "i PbarCi P0C " << i << "\t" << PbarCi[i]*_nbartrue
      //     << "\t" << P0C[i]*_nbartrue << "\t" << _nCi[i] << "\t" << delta << endl;
      diff += fabs(PbarCi[i] - P0C[i]);
    }

    //cout << "Difference " << diff << endl;

    // chi^2 (but need to calculate error)
    // takes a verrrrrrrry long time
    if (_nc*_ne < 50000) {getCovariance(_nEj);}

    // replace P0C
    P0C = PbarCi;

    if (smoothit) {smooth(PbarCi);}

    // and repeat
  }

  _trained = true;

  return(1);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::smooth(const vector<Double_t>& PbarCi) const
{
  // Smooth unfolding distributions
  if (_ndims != 1) {
    cout << "Smoothing only implemented for 1-D distributions" << endl;
    return (0);
  } else {
    cout << "Smoothing not yet implemented ! continuing." << endl;
  }

  //TH1D *h = new TH1D(PbarCi);
  return(1);
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
    Int_t iend = std::min(_nCi.size(),_nEj.size());
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
RooUnfoldBayesImpl::getCovarianceBinByBin(const vector<Double_t> effects)
{
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
  return(1);
}

//-------------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::getCovariance(const vector<Double_t>& effects)
{
  // Create the covariance matrix
  if (_nc*_ne >= 50000) 
    cout << "getCovariance (this takes some time with many bins)." << endl;
  vector<Double_t> dummy;
  Double_t nbartrue = getnbarCi(effects,dummy);

  // error due n(Ej)
  // error from data
  TStopwatch clock;
  clock.Start();
  Int_t np(0);
  for (Int_t k = 0 ; k < _nc ; k++) {
    for (Int_t l = 0 ; l < _nc ; l++) {
      _Vij->Set(k,l,0.0);
      Double_t temp(0), temp2(0);
      for (Int_t j = 0 ; j < _ne ; j++) {
        Double_t Mlj = _Mij->Get(l,j);
        Double_t Mkj = _Mij->Get(k,j);

        temp += (Mkj*Mlj*effects[j]*(1-effects[j]/nbartrue));

        for (Int_t i = 0 ; i < _ne ; i++) {
          if (i==j) {continue;}
          Double_t Mki = _Mij->Get(k,i);
          temp2 += (Mki*Mlj*effects[i]*effects[j]/nbartrue);
        }
      }

      _Vij->Set(k,l,(temp-temp2));

      //cout << "Vij " << k << " " << l << " "
      //   << _Vij[k][l] << " " << _a2d->Get(k,l) << endl;

    }
    clock.Stop();
    if (k==1 && (_nc*_ne >= 50000)) {
      Double_t nsecs  = clock.RealTime() * _nc / 2.0;
      if (nsecs>600) {
        Double_t nhours = nsecs / 3600.0;
        cout << "Predicted time per iteration : " << nhours << " hours" << endl;
      } else {
        cout << "Predicted time per iteration : " << nsecs << " secs" << endl;
      }
    } else {
      cout << "."; np++; if(np==72) {cout << endl; np=0;}
    }
  }
  cout << endl;

  return(1);
  // return the systematic error. This still needs to be debugged.
  cout << "getCovariance" << endl;
  // error due P(Ej|Ci)
  // error from Monte Carlo
  for (Int_t k = 0 ; k < _nc ; k++) {
    for (Int_t l = 0 ; l < _nc ; l++) {
      Double_t temp(0);
      for (Int_t j = 0 ; j < _ne ; j++) {
        if (_nEj[j] <= 0) {continue;}
        for (Int_t i = 0 ; i < _ne ; i++) {
          if (_nEj[i] <= 0) {continue;}

          Double_t covM(0);
          for (Int_t r = 0 ; r < _nc ; r++) {
            for (Int_t s = 0 ; s < _nc ; s++) {

              for (Int_t u = 0 ; u < _nc ; u++) {
                // cov[P(E_r|C_u),P(E_s|C_u)]
                Double_t covP(0);
                Double_t PErCu =  _Nij->Get(u,r) / _nCi[u];
                if (r==s) {
                  covP = PErCu * (1.0 - PErCu) / _nCi[u];
                } else {
                  Double_t PEsCu = _Nij->Get(u,s) / _nCi[u];
                  covP = -1.0 * PErCu * PEsCu / _nCi[u];
                }

                // cov[M_ki,M_lj]
                covM = covM + (deltaM(k,i,r,u) * deltaM(l,j,s,u) * covP);
                //cout << k << l << j << i << r << s << u << endl;
              }  // u
            }  // s
          }  // r

          //_Vij[k][l] += temp;
          _Vij->Add(k,l,temp);
        }  // i
      } // j
    } // l
  } // k

  return(1);
}

//-------------------------------------------------------------------------
Double_t
RooUnfoldBayesImpl::deltaM(Int_t k, Int_t i, Int_t r, Int_t u) const
{
  // Helper function for getCovariance.

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

  //cout << "getnbarCi nbartrue " << nbartrue << endl;
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
RooUnfoldBayesImpl::unfold(vector<Double_t>& causes)
{
  // Unfold the accumulated measured values and return the true values in a matrix.
  // causes = vector of unfolded values.

  cout << "Now unfolding..." << endl;
  // unfold
  _ncauses = getnbarCi(_nEstj, causes);
  _causes.clear();
  _causes = causes;

  // Calculate error
  if  (_nc*_ne < 50000) {getCovariance(_nEstj);}
  _nCausesError = getError();

  return(1);
}

//----------------------------------------------------------------------
Int_t
RooUnfoldBayesImpl::unfoldBinByBin(vector<Double_t>& causes)
{
  // Unfold the accumulated measured values and return the true values in a matrix.
  // causes = vector of unfolded values.

  cout << "Now unfolding (bin-by-bin)..." << endl;
  // unfold
  _ncauses = getnbarCi(_nEstj, causes);
  _causes.clear();
  _causes = causes;

  // Calculate error
  getCovarianceBinByBin(_nEstj);
  _nCausesError = getError();

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
    //cout << "bin " << i << " error " << error << endl;
    _causesError.push_back(error); // update array of errors
  }

  // total error
  if (toterror>0) {toterror = sqrt(toterror);}
  //cout << "error " << toterror << endl;

  return(toterror);
}
