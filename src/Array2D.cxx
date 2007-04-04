//-------------------------------------------------------------------------
//
// File and Version Information:
//   $Id: Array2D.cxx,v 1.1.1.1 2007-04-04 21:27:25 adye Exp $
//
// Description:
//  A helper class designed to store an 2 dimensional array as a 1-D array. 
//
// Author List:
//   F.Wilson, T.Adye
//
// Copyright Information:
//   Copyright (C) 2005-2006 rutherford Appleton Laboratory
//
//-------------------------------------------------------------------------
#include "Array2D.h"

ClassImp(Array2D);

Array2D::Array2D()
{
  // default constructor
  _nRows = _nCols = _dataSize = 0;
  _data=0;
}

// constructor
Array2D::Array2D(Int_t R, Int_t C, Double_t initval)
{
  // construct the array with R rows, C columns and initialise
  // each element with initval
  _nRows = R;
  _nCols = C;
  _dataSize = _nRows*_nCols;
  _data = new Double_t[_dataSize];
  if (_data == 0) {
    _nRows = _nCols = _dataSize = 0;
  } else {
    for (Int_t i =0 ; i < _dataSize; i++) {_data[i] = initval;}
  }
}

Array2D::Array2D(const Array2D &rhs)
{
  // copy constructor
  _nRows = rhs._nRows;
  _nCols = rhs._nCols;
  _dataSize = _nRows*_nCols;
  _data  = new Double_t[_dataSize];
  for (Int_t i =0 ; i < _dataSize; i++) {
    _data[i] = rhs._data[i];
  }
}

Bool_t
Array2D::operator==(const Array2D& rhs)
{
  // comparison operator
  if ((_nRows != rhs._nRows) || (_nCols != rhs._nCols)) {return(false);}
  Bool_t status(true);

  for (Int_t i =0 ; i < (_dataSize); i++) {
    if (_data[i] != rhs._data[i]) {status = false; break;}
  }
  return(status);
}

Array2D&
Array2D::operator=(const Array2D &rhs)
{
  // assignment operator
  if (this == &rhs) return (*this);

  delete [] _data;
  _nRows = rhs._nRows;
  _nCols = rhs._nCols;
  _dataSize = _nRows*_nCols;
  _data = new Double_t[_dataSize];

  for (Int_t i =0 ; i < _dataSize; i++) {
    _data[i] = rhs._data[i];
  }

  return (*this);
}

//
Array2D::~Array2D()
{
  delete [] _data;
  _data = 0;
}

Bool_t
Array2D::ValidRC(Int_t R, Int_t C) const
{
  // check if valid row and column number
  Bool_t status(false);
  if ((R>=0) && (R < _nRows) &&
      (C>=0) && (C < _nCols)) {status = true;}
  return(status);
}

Int_t
Array2D::Map(Int_t R, Int_t C) const
{
  // map row and column to 1D array
  return (R*_nCols+C);
}

Bool_t
Array2D::Set(Int_t R, Int_t C, Double_t value)
{
  // set the element in row R, column C to value initval
  Bool_t status(false);
  if (ValidRC(R,C)) {
    _data[Map(R,C)] = value;
    status = true;
  }
  return(status);
}

Bool_t
Array2D::Add(Int_t R, Int_t C, Double_t value)
{
  // Add value to element in row R, column C
  Bool_t status(false);
  if (ValidRC(R,C)) {
    _data[Map(R,C)] += value;
    status = true;
  }
  return(status);
}

Double_t
Array2D::Get(Int_t R, Int_t C) const
{
  // access the element in row R, column C
  Double_t value(0.0);
  if (ValidRC(R,C)) {
    value = _data[Map(R,C)];
  }
  return(value);
}
