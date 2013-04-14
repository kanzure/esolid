#ifndef _BIGRATIONAL_VECTOR_H
#define _BIGRATIONAL_VECTOR_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>

using namespace std;

class bigrational_vector
{
private:
  
  unsigned long dim;
  bigrational*  rep;
  
  //  arithmetic
  
  bigrational_vector add(const bigrational_vector&) const;
  bigrational_vector sub(const bigrational_vector&) const;
  
  bigrational_vector neg() const;
  
  bigrational        in_prod(const bigrational_vector&) const;
  bigrational_vector scalar_mul(const bigrational&)     const;
  
  //  comparison
  
  int cmp(const bigrational_vector&) const;
  
  int is_zero() const;
  
public:
  
  //  constructor, assignment and destructor
  
  bigrational_vector(const unsigned long = 1);
  bigrational_vector(const unsigned long, const bigrational*);
//  bigrational_vector(const bigint_vector&);
  
  bigrational_vector(const bigrational_vector&);
  bigrational_vector& operator =(const bigrational_vector&);
  
  ~bigrational_vector();
  
  //  size
  
  unsigned long get_dim() const;
  
  //  element
  
  bigrational& operator [](const unsigned long) const;
  
  //  arithmetic
  
  friend bigrational_vector operator +(const bigrational_vector&,
                                       const bigrational_vector&);
  friend bigrational_vector operator -(const bigrational_vector&,
                                       const bigrational_vector&);
  
  friend bigrational_vector operator -(const bigrational_vector&);
  
  friend bigrational        in_prod(const bigrational_vector&,
                                    const bigrational_vector&);
  friend bigrational_vector scalar_mul(const bigrational&,
                                       const bigrational_vector&);
  
  //  comparison
  
  friend int operator ==(const bigrational_vector&, const bigrational_vector&);
  friend int operator !=(const bigrational_vector&, const bigrational_vector&);
  
  friend int is_zero(const bigrational_vector&);
  
  //  arithemtic and assignment
  
  bigrational_vector& operator +=(const bigrational_vector&);
  bigrational_vector& operator -=(const bigrational_vector&);
  
//  bigrational_vector& scalar_mul_assign(const bigrational&);
  
  //  stream
  
  friend ostream& operator <<(ostream&, const bigrational_vector&);
};

inline unsigned long bigrational_vector :: get_dim() const
{
  return dim;
}

inline
bigrational& bigrational_vector :: operator [](const unsigned long n) const
{
  assert(n < dim);
  
  return rep[n];
}

#endif

