//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: ArgVars.h,v 1.2 2010-01-19 15:33:45 adye Exp $
//
// Description:
//      Parse argument list for parameter settings
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ARGVARS_H
#define ARGVARS_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include <iostream>
#include "TObject.h"
#include "TList.h"
#endif

class TIter;

class ArgVar : public TObject {
public:
  const char *name;
  Int_t      *ivar;
  Double_t   *fvar;
  ArgVar (const char *n,    Int_t *i) : name(n), ivar(i), fvar(0) {}
  ArgVar (const char *n, Double_t *f) : name(n), ivar(0), fvar(f) {}
};

class ArgVars {
  TList lst;
public:
  ArgVars() {}
  ~ArgVars();
  void  Add (const char* name,    Int_t* ivar) { lst.Add (new ArgVar (name, ivar)); }
  void  Add (const char* name, Double_t* fvar) { lst.Add (new ArgVar (name, fvar)); }
  void  Add (const ArgVar&  arg)               { lst.Add (arg.Clone());             }
  void  Add (const ArgVars& args);
  Int_t SetArgs (int argc, const char* const* argv, bool split= false);
  void  Print (std::ostream& o, const char* sep= " ") const;
  void  Usage (const char* prog);
};

#endif
