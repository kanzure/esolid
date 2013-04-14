#ifndef _KSOLID_H
#define _KSOLID_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <kpatch.h>
#include <kpartition.h>
#include <kgraph.h>
#include <bigrational_matrix.h>

using namespace std;

class K_SOLID
{
  unsigned long num_patches;
  K_PATCH**     patches;
  
  //  stream
  
  ostream& output(ostream&) const;
  
public:
  
  //  constructors, assignment and destructor
  
  K_SOLID();
  K_SOLID(K_PATCH* const [], const unsigned long);
  
  K_SOLID(const K_SOLID&);
  K_SOLID& operator =(const K_SOLID&);
  
  ~K_SOLID();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_SOLID&);
  
  K_SOLID transform(const bigrational_matrix&) const;
  int     Bezier_output(ostream&) const;
  
  //  primitives
  
  //  int classify_pt(const bigrational& x,
  //                  const bigrational& y,
  //                  const bigrational& z) const
  //    returns IN  if (x, y, z) lies inside *this and
  //            OUT if (x, y, z) lies outside *this.
  
  int classify_pt(const bigrational&,
                  const bigrational&,
                  const bigrational&) const;
  
  //  int classify_pt(const bigrational_vector& pt) const
  //    returns IN  if pt lies inside *this and
  //            OUT if pt lies outside *this.
  
  int classify_pt(const bigrational_vector&) const;
  
  //  arithmetic
  
  K_SOLID boolean(K_SOLID&, const char);
  K_SOLID merge(const K_SOLID&) const;
  
  //  other functions
  
  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
                               K_PARTITION**, const unsigned long);
};

K_SOLID read_solid(istream&, const bigrational& = 0);
K_SOLID read_CSG(const char*, const char*, const bigrational& = 0);

bigrational_matrix read_BRLCAD_matrix(istream&);

#endif

