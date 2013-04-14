//  file:    bigrational.h
//  update:  09/25/02

#ifndef _BIGRATIONAL_H
#define _BIGRATIONAL_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <gmp.h>
//#include <bigrand.h>
#include <bigint.h>

using namespace std;

class bigrational
{
friend class bigreal;
friend class bigcomplex;
  
  mpq_t rep;
  
  //  conversion
  
  bigint num() const;
  bigint den() const;
  
  bigint ceil() const;
  bigint floor() const;
  bigint trunc() const;
  
  long lgAbs() const;
  
  //  arithmetic
  
  bigrational add(const bigrational&) const;
  bigrational sub(const bigrational&) const;
  bigrational mul(const bigrational&) const;
  bigrational div(const bigrational&) const;
  
  bigrational neg() const;
  bigrational abs() const;
  bigrational inv() const;
  
  //  comparison
  
  int equal(const bigrational&) const;
  int cmp(const bigrational&) const;
  
  int sgn() const;
  
public:
  
  //  constructor, assignment and destructor
  
  bigrational();
  bigrational(const unsigned int);
  bigrational(const signed int);
  bigrational(const unsigned long);
  bigrational(const signed long);
  bigrational(const double);
  bigrational(const bigint&);
  bigrational(const bigint&, const bigint&);
  
  bigrational(const bigrational&);
  bigrational& operator =(const bigrational&);
  
  ~bigrational();
  
  //  random
  
//  friend bigrational DRand();
  
  //  conversion
  
  unsigned long len() const;
  
  double as_double() const;
  
  int    fits_bigint() const;
  bigint as_bigint() const;
  
  friend bigint num(const bigrational&);
  friend bigint den(const bigrational&);
  
  friend bigint ceil(const bigrational&);
  friend bigint floor(const bigrational&);
  friend bigint trunc(const bigrational&);
  
  friend long lb4lgAbs(const bigrational&);
  friend long ub4lgAbs(const bigrational&);
  
  //  arithmetic
  
  friend bigrational operator +(const bigrational&, const bigrational&);
  friend bigrational operator -(const bigrational&, const bigrational&);
  friend bigrational operator *(const bigrational&, const bigrational&);
  friend bigrational operator /(const bigrational&, const bigrational&);
  
  friend bigrational operator -(const bigrational&);
  friend bigrational abs(const bigrational&);
  friend bigrational inv(const bigrational&);
  
  //  comparison
  
  friend int operator ==(const bigrational&, const bigrational&);
  friend int operator !=(const bigrational&, const bigrational&);
  friend int operator <(const bigrational&, const bigrational&);
  friend int operator <=(const bigrational&, const bigrational&);
  friend int operator >(const bigrational&, const bigrational&);
  friend int operator >=(const bigrational&, const bigrational&);
  
  friend int sgn(const bigrational&);
  friend int operator &&(const bigrational&, const bigrational&);
  friend int operator ||(const bigrational&, const bigrational&);
  friend int operator !(const bigrational&);
  
  //  arithemtic and assignment
  
  bigrational& operator +=(const bigrational&);
  bigrational& operator -=(const bigrational&);
  bigrational& operator *=(const bigrational&);
  
  bigrational& operator /=(const bigrational&);
  
  //  stream
  
  friend ostream& operator <<(ostream&, const bigrational&);
  friend istream& operator >>(istream&, bigrational&);
};

#endif

