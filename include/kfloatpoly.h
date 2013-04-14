#ifndef _KFLOATPOLY_H
#define _KFLOATPOLY_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>

using namespace std;

class K_FLOATPOLY
{
//  A class for representing polynomials with bigrational coefficients
  
  unsigned long num_vars;    //  Number of variables
  long*         deg;         //  Max. degree in each variable
  
  unsigned long num_coeffs;  //  Number of coefficients (terms)
  double*       coeffs;
  
  //  primitives
  
  long* index_to_powers(const unsigned long i) const;
  //  For an index, i, return the powers, p, of the term
  //  to which the i-th coefficient belongs to.
  
  unsigned long powers_to_index(const long* const p) const;
  //  Return the index of the term of the powers p.
  
  int reduce_deg();
  
  //  arithmetic
  
  K_FLOATPOLY add(const K_FLOATPOLY&) const;
  K_FLOATPOLY sub(const K_FLOATPOLY&) const;
  K_FLOATPOLY mul(const K_FLOATPOLY&) const;
  K_FLOATPOLY mul(const double) const;
  
	K_FLOATPOLY neg() const;
  
  //  stream
  
  ostream& output(ostream&) const;
  
public:
  
  //  constructors, assignment and destructor
  
	K_FLOATPOLY();
  K_FLOATPOLY(const unsigned long, const long* const);
  K_FLOATPOLY(const unsigned long, const long* const, const double* const);
  
  K_FLOATPOLY(const K_FLOATPOLY&);
  K_FLOATPOLY& operator =(const K_FLOATPOLY&);
  
  ~K_FLOATPOLY();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_FLOATPOLY&);
  
	//  arithmetic
  
  friend K_FLOATPOLY operator +(const K_FLOATPOLY&, const K_FLOATPOLY&);
  friend K_FLOATPOLY operator -(const K_FLOATPOLY&, const K_FLOATPOLY&);
  friend K_FLOATPOLY operator *(const K_FLOATPOLY&, const K_FLOATPOLY&);
  friend K_FLOATPOLY operator *(const K_FLOATPOLY&, const double);
  friend K_FLOATPOLY operator *(const double, const K_FLOATPOLY&);
  
  friend K_FLOATPOLY operator -(const K_FLOATPOLY&);
  
  //  others
  
  K_FLOATPOLY   derivative(const unsigned long) const;
  K_FLOATPOLY   subst_first_var(const double) const;
  double        evaluate(const double* const) const;
  K_FLOATPOLY   subst_val(const unsigned long, const double) const;
  K_FLOATPOLY   subst_expr(const unsigned long, const K_FLOATPOLY&) const;
  K_FLOATPOLY   subst_expr(const unsigned long,
                           const K_FLOATPOLY&, const K_FLOATPOLY&) const;
  
  //  other univariate functions
  
  //  unsigned long gen_fp_roots(const double l,
  //                             const doubld h,
  //                             double*& X) const
  //    PROVIDED *this is a univariate polynomial,
  //    compute roots X of *this in [l, h] by using only FP arithmetic.
  
  unsigned long gen_fp_roots(const double, const double, double*&) const;
};

#endif

