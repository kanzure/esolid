#ifndef _KINTPOLY_H
#define _KINTPOLY_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigint.h>

using namespace std;

class K_INTPOLY
{
//  A class for polynomials with bigint coefficients

  unsigned long num_vars;    //  Number of variables
  long*         deg;         //  Max. degree in each variable

  unsigned long num_coeffs;  //  Number of coefficients (terms)
  bigint*       coeffs;

  //  primitives

  long* index_to_powers(const unsigned long i) const;
  //  For an index, i, return the powers, p, of the term
  //  to which the i-th coefficient belongs to.

  unsigned long powers_to_index(const long* const p) const;
  //  Return the index of the term of the powers p.

  int reduce_deg();

  //  arithmetic

  K_INTPOLY add(const K_INTPOLY&) const;
  K_INTPOLY sub(const K_INTPOLY&) const;
  K_INTPOLY mul(const K_INTPOLY&) const;
  K_INTPOLY mul(const bigint&) const;

	K_INTPOLY neg() const;

  //  stream

  ostream& output(ostream&) const;

public:

  //  constructors, assignment and destructor

	K_INTPOLY();
  K_INTPOLY(const unsigned long, const long* const);
  K_INTPOLY(const unsigned long, const long* const, const bigint* const);

  K_INTPOLY(const K_INTPOLY&);
  K_INTPOLY& operator =(const K_INTPOLY&);

  ~K_INTPOLY();

  //  stream

  friend ostream& operator <<(ostream&, const K_INTPOLY&);

	//  arithmetic

  friend K_INTPOLY operator +(const K_INTPOLY&, const K_INTPOLY&);
  friend K_INTPOLY operator -(const K_INTPOLY&, const K_INTPOLY&);
  friend K_INTPOLY operator *(const K_INTPOLY&, const K_INTPOLY&);
  friend K_INTPOLY operator *(const K_INTPOLY&, const bigint&);
  friend K_INTPOLY operator *(const bigint&, const K_INTPOLY&);

  friend K_INTPOLY operator -(const K_INTPOLY&);

  //  other functions

  K_INTPOLY derivative(const unsigned long) const;
  K_INTPOLY subst_first_var(const bigint&) const;
  bigint    evaluate(const bigint* const) const;
  int       sgn_at(const bigint* const) const;
  K_INTPOLY subst_val(const unsigned long, const bigint&) const;
  K_INTPOLY subst_expr(const unsigned long, const K_INTPOLY&) const;
  K_INTPOLY subst_expr(const unsigned long,
                       const K_INTPOLY&, const K_INTPOLY&) const;

  //  other univariate functions

  int sgn_at(const bigint&) const;
};

#endif

