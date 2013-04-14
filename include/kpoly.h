#ifndef _KPOLY_H
#define _KPOLY_H

#include <cassert>
#include <cstdlib>
#include <iostream>

using namespace std;

class K_POLY
{
//  A base class for polynomials
  
protected:
  
  unsigned long num_vars;    //  Number of variables
  long*         deg;         //  Max. degree in each variable
  unsigned long num_coeffs;  //  Number of coefficients (terms)
  
public:
  
  //  constructor, assignment and destructor
  
	K_POLY();
  K_POLY(const unsigned long, const long* const);
  
  K_POLY(const K_POLY&);
  K_POLY& operator =(const K_POLY&);
  
  virtual ~K_POLY();
  
protected:
  
  long* index_to_powers(const unsigned long i) const;
  //  For an index, i, return the powers, p, of the term
  //  to which the i-th coefficient belongs to.
  
  long* index_to_powers_homog(const unsigned long t,
                              const unsigned long i) const;
  //  For an index, i, return the powers, p, of the term
  //  to which the i-th coefficient belongs to.
  
  unsigned long index_to_deg(const unsigned long i) const;
  //  Return the total degree of the i-th term.
  
  unsigned long powers_to_index(const long* const p) const;
  //  Return the index of the term of the powers p.
};

#endif

