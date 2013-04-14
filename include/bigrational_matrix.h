//  file:    bigrational_matrix.h
//  update:  09/25/02

#ifndef _BIGRATIONAL_MATRIX_H
#define _BIGRATIONAL_MATRIX_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>

using namespace std;

class bigrational_matrix
{
private:
  
  unsigned long num_row;
  unsigned long num_col;
  bigrational*  rep;
  
  //  arithmetic
  
  bigrational_matrix add(const bigrational_matrix&) const;
  bigrational_matrix sub(const bigrational_matrix&) const;
  bigrational_matrix mul(const bigrational_matrix&) const;
  
  bigrational_matrix neg() const;
  
//  bigrational_matrix scalar_mul(const bigrational&) const;
  
  //  comparison
  
  int cmp(const bigrational_matrix&) const;
  
  //  linear algebra
  
  unsigned long rank() const;
  
  //  square matrix
  
  bigrational det() const;
  
public:
  
  //  constructor, assignment and destructor
  
  bigrational_matrix(const unsigned long = 1, const unsigned long = 1);
  bigrational_matrix(const unsigned long, const unsigned long,
                     const bigrational*);
//  bigrational_matrix(const bigint_matrix&);
  
  bigrational_matrix(const bigrational_matrix&);
  bigrational_matrix& operator =(const bigrational_matrix&);
  
  ~bigrational_matrix();
  
  //  size
  
  unsigned long get_num_row() const;
  unsigned long get_num_col() const;
  
  //  element
  
  bigrational& operator ()(const unsigned long, const unsigned long) const;
  
  //  arithmetic
  
  friend bigrational_matrix operator +(const bigrational_matrix&,
                                       const bigrational_matrix&);
  friend bigrational_matrix operator -(const bigrational_matrix&,
                                       const bigrational_matrix&);
  friend bigrational_matrix operator *(const bigrational_matrix&,
                                       const bigrational_matrix&);
  
  friend bigrational_matrix operator -(const bigrational_matrix&);
  
//  friend bigrational_matrix scalar_mul(const bigrational&,
//                                       const bigrational_matrix&);
  
  //  comparison
  
  friend int operator ==(const bigrational_matrix&, const bigrational_matrix&);
  friend int operator !=(const bigrational_matrix&, const bigrational_matrix&);
  
  //  arithemtic and assignment
  
  bigrational_matrix& operator +=(const bigrational_matrix&);
  bigrational_matrix& operator -=(const bigrational_matrix&);
  
//  bigrational_matrix& scalar_mul_assign(const bigrational&);
  
  //  linear algebra
  
  friend unsigned long rank(const bigrational_matrix&);
  
  //  square matrix
  
  int                is_identity() const;
  friend bigrational det(const bigrational_matrix&);
  bigrational        Bareiss(bigrational_matrix&) const;
  bigrational_matrix inverse() const;
  
  //  stream
  
  friend ostream& operator <<(ostream&, const bigrational_matrix&);
};

inline unsigned long bigrational_matrix :: get_num_row() const
{
  return num_row;
}

inline unsigned long bigrational_matrix :: get_num_col() const
{
  return num_col;
}

inline
bigrational& bigrational_matrix :: operator ()(const unsigned long r,
                                               const unsigned long c) const
{
  return rep[r * num_col + c];
}

#endif

