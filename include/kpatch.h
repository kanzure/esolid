#ifndef _KPATCH_H
#define _KPATCH_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <bigrational_vector.h>
#include <kratpoly.h>
#include <kcurve.h>
#include <ksurf.h>
#include <kpartition.h>

using namespace std;

class K_GRAPH;

class K_PATCH
{
friend class K_PARTITION;
friend class K_SOLID;

  unsigned long ID;

  K_SURF*     surf;     //  surface
  bool        is_head;  //  true if *this is the head of surf
  bigrational low_s;    //  lower bound for *this in s
  bigrational high_s;   //  upper bound for *this in s
  bigrational low_t;    //  lower bound for *this in t
  bigrational high_t;   //  upper bound for *this in t

  //  trimming curves

  unsigned long num_trim_curves;
  K_CURVE**     trim_curves;
  K_SURF**      adj_surfs;
  K_PATCH**     adj_patches;

  //  intersection curves

  unsigned long num_int_curves;
  K_CURVE**     int_curves;
  K_SURF**      adj_int_surfs;
  K_PATCH**     adj_int_patches;

  unsigned long  num_merged;
  K_CURVE***     ic;
  unsigned long* len_ic;
  int*           ic_closed;

  unsigned long ref_count;

  //  stream

  ostream& output(ostream&) const;

  //  primitives

  int set_range();

public:

  //  constructors, assignment and destructor

  K_PATCH();
//  K_PATCH(K_SURF* const, K_CURVE* const [], const unsigned long);
  K_PATCH(K_SURF* const, K_CURVE* const [], const unsigned long,
          const bool = true);
  K_PATCH(const K_PARTITION&);

  K_PATCH(const K_PATCH&);
  K_PATCH& operator =(const K_PATCH&);

  ~K_PATCH();

  //  stream

  friend ostream& operator <<(ostream&, const K_PATCH&);

  int Bezier_output(ostream&,
                    const unsigned int,
                    const unsigned int,
                    const unsigned int) const;

  //  primitives

  bigrational get_low_s() const;
  bigrational get_low_t() const;
  bigrational get_high_s() const;
  bigrational get_high_t() const;

  //  comparisons

  int in_dom(K_POINT2D&);
  int contains(K_POINT2D&);
  int in_on_out(K_POINT2D&);

  //  other functions

  int           intersect(K_PATCH&);
  int           split_trim_curves(const unsigned long, const unsigned long);
  unsigned long merge_curves();
  unsigned long split_loops(K_PATCH**&);

  K_POINT2D get_pt_in() const;

  friend unsigned long gen_partitions(K_PATCH* const, K_PARTITION**&);

  friend unsigned long gen_adjacency(K_PARTITION** const, const unsigned long,
                                     K_GRAPH*&,
                                     long*&, K_GRAPH*&);

//  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
//                               K_PARTITION**, const unsigned long);
  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
                               K_PARTITION**, const unsigned long,
                               const char);

  friend K_SOLID gen_box(const bigrational_vector [], const unsigned long);

  friend K_SOLID gen_cyl(const bigrational_vector&, const bigrational_vector&,
                         const bigrational_vector&, const bigrational_vector&,
                         const bigrational_vector&, const bigrational_vector&);

  friend unsigned long get_patches_ell(const bigrational_vector [4],
                                       K_PATCH**&);

  friend K_SOLID gen_ell(const bigrational_vector [4]);

  friend int get_patch_tor(K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_PATCH*&);

  friend K_SOLID gen_tor(const bigrational_vector&,
                         const bigrational_vector&,
                         const bigrational_vector&,
                         const bigrational_vector&,
                         const bigrational&, const bigrational&);
};

#endif

