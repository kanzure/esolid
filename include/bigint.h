#ifndef _BIGINT_H
#define _BIGINT_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <gmp.h>
//#include <bigrand.h>

using namespace std;

class bigint
{
friend class bigrational;
friend class bigreal;
friend class bigcomplex;

  mpz_t rep;

  //  constructor, assignment and destructor

  bigint(const mpz_t);

  //  conversion

  double frexp(long*) const;

  //  arithmetic

  bigint add(const bigint&) const;
  bigint sub(const bigint&) const;
  bigint mul(const bigint&) const;
  bigint div(const bigint&) const;

  void inc();
  void dec();

  bigint lshift(const unsigned long) const;
  bigint rshift(const unsigned long) const;

  bigint neg() const;
  bigint abs() const;

  bigint root(const unsigned long) const;
  bigint sqrt() const;

  //  comparison

  int cmp(const bigint&) const;

  int sgn() const;

  //  bit

  bigint AND(const bigint&) const;
  bigint IOR(const bigint&) const;
  bigint XOR(const bigint&) const;

  //  algebra

  bigint gcd(const bigint&) const;
  bigint lcm(const bigint&) const;

public:

  //  constructor, assignment and destructor

  bigint();
  bigint(const unsigned int);
  bigint(const signed int);
  bigint(const unsigned long);
  bigint(const signed long);
  bigint(const double);
  bigint(const char*);

  bigint(const bigint&);
  bigint& operator =(const bigint&);

  ~bigint();

  //  random

//  friend bigint rand();

  //  conversion

  unsigned long len() const;

  int           fits_long() const;
  long          as_long() const;
  int           fits_unsigned_long() const;
  unsigned long as_unsigned_long() const;

  friend double frexp(const bigint&, long*);
  friend long   lb4lgAbs(const bigint&);
  friend long   ub4lgAbs(const bigint&);

  //  arithmetic

  friend bigint operator +(const bigint&, const bigint&);
  friend bigint operator -(const bigint&, const bigint&);
  friend bigint operator *(const bigint&, const bigint&);
  friend bigint operator /(const bigint&, const bigint&);

  bigint& operator ++();     //  prefix
  bigint  operator ++(int);  //  postfix
  bigint& operator --();     //  prefix
  bigint  operator --(int);  //  postfix

  friend bigint operator <<(const bigint&, const unsigned long);
  friend bigint operator >>(const bigint&, const unsigned long);

  friend bigint operator -(const bigint&);
  friend bigint abs(const bigint&);

  friend bigint root(const bigint&, const unsigned long);
  friend bigint sqrt(const bigint&);

  //  comparison

  friend int operator ==(const bigint&, const bigint&);
  friend int operator !=(const bigint&, const bigint&);
  friend int operator <(const bigint&, const bigint&);
  friend int operator <=(const bigint&, const bigint&);
  friend int operator >(const bigint&, const bigint&);
  friend int operator >=(const bigint&, const bigint&);

  friend int sgn(const bigint&);
  friend int operator &&(const bigint&, const bigint&);
  friend int operator ||(const bigint&, const bigint&);
  friend int operator !(const bigint&);

  //  bit

  friend bigint operator &(const bigint&, const bigint&);
  friend bigint operator |(const bigint&, const bigint&);
  friend bigint operator ^(const bigint&, const bigint&);

  int bit(unsigned long) const;

  //  algebra

  friend bigint gcd(const bigint&, const bigint&);
  friend bigint lcm(const bigint&, const bigint&);

  //  arithmetic and assignment

  bigint& operator +=(const bigint&);
  bigint& operator -=(const bigint&);
  bigint& operator *=(const bigint&);
  bigint& operator /=(const bigint&);

  bigint& operator <<=(const unsigned long);
  bigint& operator >>=(const unsigned long);

  bigint& operator &=(const bigint&);
  bigint& operator |=(const bigint&);
  bigint& operator ^=(const bigint&);

  //  stream

  friend ostream& operator <<(ostream&, const bigint&);
  friend istream& operator >>(istream&, bigint&);
};

#endif

