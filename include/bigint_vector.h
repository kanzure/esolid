#ifndef _BIGINT_VECTOR_H
#define _BIGINT_VECTOR_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigint.h>

using namespace std;

class bigint_vector
{
friend class bigrational_vector;
  
private:
  
  unsigned long dim;
  bigint*       rep;
  
  //  arithmetic
  
  bigint_vector add(const bigint_vector&) const;
  bigint_vector sub(const bigint_vector&) const;
  
  bigint_vector neg() const;
  
  bigint        in_prod(const bigint_vector&) const;
  bigint_vector scalar_mul(const bigint&)     const;
  
  //  comparison
  
  int cmp(const bigint_vector&) const;
  
  int is_zero() const;
  
public:
  
  //  constructor, assignment and destructor
  
  bigint_vector(const unsigned long = 1);
  bigint_vector(const unsigned long, const bigint*);
  
  bigint_vector(const bigint_vector&);
  bigint_vector& operator =(const bigint_vector&);
  
  ~bigint_vector();
  
  //  size
  
  unsigned long get_dim() const;
  
  //  element
  
  bigint& operator [](const unsigned long) const;
  
  //  arithmetic
  
  friend bigint_vector operator +(const bigint_vector&, const bigint_vector&);
  friend bigint_vector operator -(const bigint_vector&, const bigint_vector&);
  
  friend bigint_vector operator -(const bigint_vector&);
  
  friend bigint        in_prod(const bigint_vector&, const bigint_vector&);
  friend bigint_vector scalar_mul(const bigint&, const bigint_vector&);
  
  //  comparison
  
  friend int operator ==(const bigint_vector&, const bigint_vector&);
  friend int operator !=(const bigint_vector&, const bigint_vector&);
  
  friend int is_zero(const bigint_vector&);
  
  //  arithemtic and assignment
  
  bigint_vector& operator +=(const bigint_vector&);
  bigint_vector& operator -=(const bigint_vector&);
  
//  bigint_vector& scalar_mul_assign(const bigint&);
  
  //  stream
  
  friend ostream& operator <<(ostream&, const bigint_vector&);
};

inline unsigned long bigint_vector :: get_dim() const
{
  return dim;
}

inline bigint& bigint_vector :: operator [](const unsigned long n) const
{
  assert(n < dim);
  
  return rep[n];
}

#endif

