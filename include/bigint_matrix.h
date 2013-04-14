//  file:    bigint_matrix.h
//  update:  09/25/02

#ifndef _BIGINT_MATRIX_H
#define _BIGINT_MATRIX_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigint.h>

using namespace std;

class bigint_matrix
{
friend class bigrational_matrix;
  
private:
  
  unsigned long num_row;
  unsigned long num_col;
  bigint*       rep;
  
  //  arithmetic
  
  bigint_matrix add(const bigint_matrix&) const;
  bigint_matrix sub(const bigint_matrix&) const;
  bigint_matrix mul(const bigint_matrix&) const;
  
  bigint_matrix neg() const;
  
//  bigint_matrix scalar_mul(const bigint&) const;
  
  //  comparison
  
  int cmp(const bigint_matrix&) const;
  
  //  square matrix
  
  int    is_identity() const;
  bigint det() const;
  
public:
  
  //  constructor, assignment and destructor
  
  bigint_matrix(const unsigned long = 1, const unsigned long = 1);
  bigint_matrix(const unsigned long, const unsigned long, const bigint*);
  
  bigint_matrix(const bigint_matrix&);
  bigint_matrix& operator =(const bigint_matrix&);
  
  ~bigint_matrix();
  
  //  size
  
  unsigned long get_num_row() const;
  unsigned long get_num_col() const;
  
  //  element
  
  bigint& operator ()(const unsigned long, const unsigned long) const;
  
  //  arithmetic
  
  friend bigint_matrix operator +(const bigint_matrix&, const bigint_matrix&);
  friend bigint_matrix operator -(const bigint_matrix&, const bigint_matrix&);
  friend bigint_matrix operator *(const bigint_matrix&, const bigint_matrix&);
  
  friend bigint_matrix operator -(const bigint_matrix&);
  
//  friend bigint_matrix scalar_mul(const bigint&, const bigint_matrix&);
  
  //  comparison
  
  friend int operator ==(const bigint_matrix&, const bigint_matrix&);
  friend int operator !=(const bigint_matrix&, const bigint_matrix&);
  
  //  arithemtic and assignment
  
  bigint_matrix& operator +=(const bigint_matrix&);
  bigint_matrix& operator -=(const bigint_matrix&);
  
//  bigint_matrix& scalar_mul_assign(const bigint&);
  
  //  square matrix
  
  friend bigint det(const bigint_matrix&);
  
  //  stream
  
  friend ostream& operator <<(ostream&, const bigint_matrix&);
};

inline unsigned long bigint_matrix :: get_num_row() const
{
  return num_row;
}

inline unsigned long bigint_matrix :: get_num_col() const
{
  return num_col;
}

inline bigint& bigint_matrix :: operator ()(const unsigned long r,
                                            const unsigned long c) const
{
  return rep[r * num_col + c];
}

#endif

