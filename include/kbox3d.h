#ifndef _KBOX3D_H
#define _KBOX3D_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>

using namespace std;

class K_BOX3D
{
friend class K_SURF;
friend class K_SOLID;

  bigrational low[3];
  bigrational high[3];

  int low_infty[3];   //  - 1 iff low[i] or high[i] is - infty
  int high_infty[3];  //    1 iff low[i] or high[i] is + infty
                      //    0 otherwise, i.e., low[i] or high[i] is bigrational

  //  stream

  ostream& output(ostream&) const;

  //  comparisons
public:
  int overlap(const K_BOX3D&) const;
  int contains(const K_BOX3D&) const;

//public:

  //  constructors, assignment and destructor

  K_BOX3D();
  K_BOX3D(const bigrational&, const bigrational&,
          const bigrational&, const bigrational&,
          const bigrational&, const bigrational&);
  K_BOX3D(const bigrational* const, const bigrational* const);

  K_BOX3D(const K_BOX3D&);
  K_BOX3D& operator =(const K_BOX3D&);

  ~K_BOX3D();

  //  stream

  friend ostream& operator <<(ostream&, const K_BOX3D&);
};

#endif

