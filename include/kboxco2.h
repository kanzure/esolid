#ifndef _KBOXCO2_H
#define _KBOXCO2_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>

using namespace std;

class K_POINT2D;
class K_CURVE;

class K_BOXCO2
{
friend class K_SEGMENT;
friend class K_CURVE;
friend class K_PATCH;
friend class K_PARTITION;

  bigrational low[2];
  bigrational high[2];

  int low_open[2];
  int high_open[2];

  //  stream

  ostream& output(ostream&) const;

  //  comparison

  int overlap(const K_BOXCO2&) const;
  int contains(const K_BOXCO2&) const;

  //  other

  K_BOXCO2 merge(const K_BOXCO2&) const;

public:

  //  constructors, assignment and destructor

  K_BOXCO2();
  K_BOXCO2(const bigrational&, const bigrational&,
           const bigrational&, const bigrational&,
           const int, const int, const int, const int);
  K_BOXCO2(const bigrational* const, const bigrational* const,
           const int* const, const int* const);

  K_BOXCO2(const K_BOXCO2&);
  K_BOXCO2& operator =(const K_BOXCO2&);

  ~K_BOXCO2();

  //  stream

  friend ostream& operator <<(ostream&, const K_BOXCO2&);

  //  primitive

  bigrational get_low_s() const;
  bigrational get_high_s() const;
  bigrational get_low_t() const;
  bigrational get_high_t() const;

  int is_low_s_open() const;
  int is_high_s_open() const;
  int is_low_t_open() const;
  int is_high_t_open() const;

  //  other

  friend int pt_inside_trim_curves(K_POINT2D&,
                                   K_CURVE** const, const unsigned long);
  friend int pt_in_on_out_trim_curves(K_POINT2D&,
                                      K_CURVE** const, const unsigned long);
};

#endif

