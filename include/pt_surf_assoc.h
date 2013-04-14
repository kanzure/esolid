//  file:    pt_surf_assoc.h
//  update:  12/05/02

#ifndef PT_SURF_ASSOC_H
#define PT_SURF_ASSOC_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <kpoint2d.h>
#include <ksurf.h>

using namespace std;

class PT_SURF_ASSOC
{
  unsigned long len;
  K_POINT2D**   pts;
  K_SURF**      surfs;
  long*         rep;
  
  unsigned long find(const K_POINT2D*, const K_SURF*);
  unsigned long locate(K_POINT2D* const, K_SURF* const);
  
public:
  
  //  constructors and the destructor
  
  PT_SURF_ASSOC();
  ~PT_SURF_ASSOC();
  
  //  other functions
  
  unsigned long record_as_assoc(K_POINT2D* const x, K_SURF* const y,
                                K_POINT2D* const z, K_SURF* const w);
  //  associate (x, y) and (z, w)
  
  K_POINT2D* find_assoc_pt(const K_POINT2D* x, const K_SURF*y,
                           const K_SURF* w);
  //  return z where (x, y) and (z, w) are associated
};

#endif

