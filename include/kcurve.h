#ifndef _KCURVE_H
#define _KCURVE_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <kratpoly.h>
#include <kpoint2d.h>
#include <ksegment.h>

using namespace std;

class K_GRAPH;
class K_SOLID;
class K_SEGMENT;

class K_CURVE
{
friend class K_PATCH;
friend class K_PARTITION;
friend class K_SOLID;
  
  K_RATPOLY*    poly;
  
  K_SEGMENT**   segments;
  unsigned long num_segments;
  
  K_CURVE*      curve_in_other_dom;
  int           dir_in_other_dom;
  
  unsigned long ref_count;
  
  //  stream
  
  ostream& output(ostream&) const;
  
  //  primitives
  
  int assoc(K_CURVE* const, const int);
  int reverse();
  
  int subdivide(const unsigned long);
  
public:
  
  //  constructors, assignment and destructor
  
  K_CURVE();
  K_CURVE(const K_RATPOLY&, K_SEGMENT* const [], const unsigned long);
  K_CURVE(K_RATPOLY* const, K_SEGMENT* const [], const unsigned long);
  K_CURVE(const K_CURVE&, const unsigned long, const unsigned long);
  
  K_CURVE(const K_CURVE&);
  K_CURVE& operator =(const K_CURVE&);
  
  ~K_CURVE();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_CURVE&);
  
  //  primitives
  
  K_POINT2D* start() const;
  K_POINT2D* end() const;
  K_BOXCO2   bbox() const;
  int        is_closed() const;
  int        rotate_closed_curve(const unsigned long);
  
  //  comparisons
  
  int sort_pts(K_POINT2D** const, const unsigned long) const;
  
  //  other functions
  
  unsigned long locate_pt_seg_start(K_POINT2D&);
  unsigned long locate_pt_seg_end(K_POINT2D&);
  int           is_start_or_end(K_POINT2D&);
  int           contains(K_POINT2D&, const int = 1);
  K_POINT2D     pt_on();
  int           add_pt(K_POINT2D* const, const int = 0);
  unsigned long find_intersections(const K_RATPOLY&, K_POINT2D**&, const int);
  int           add_on(const K_CURVE&);
  int           split(const unsigned long, K_CURVE&, K_CURVE&);
  
  int           mk_seg_fw(int* const, const int,
                          const long,
                          K_POINT2D** const, const unsigned long,
                          K_POINT2D** const, const unsigned long) const;
  int           mk_seg_bw(int* const, const int,
                          const long,
                          K_POINT2D** const, const unsigned long,
                          K_POINT2D** const, const unsigned long) const;
  
  friend int pt_inside_trim_curves(K_POINT2D&,
                                   K_CURVE** const, const unsigned long);
  friend int pt_in_on_out_trim_curves(K_POINT2D&,
                                      K_CURVE** const, const unsigned long);
  
  //  gen_curve_topo()
  
  friend unsigned long gen_curve_topo(const K_RATPOLY&,
                                      const bigrational&, const bigrational&,
                                      const bigrational&, const bigrational&,
                                      K_CURVE**&);
  
  friend unsigned long gen_curve_topo_proto(const K_RATPOLY&,
                                            const bigrational&,
                                            const bigrational&,
                                            const bigrational&,
                                            const bigrational&,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_CURVE**&);
  
  //  gen_box(), gen_cyl(), gen_ell(), gen_tor()
  
  friend int get_patch4(const bigrational_vector [4], K_PATCH*&);
  
  friend K_SOLID gen_box(const bigrational_vector [], const unsigned long);
  
  friend int get_cyl_cap(const bigrational_vector&,
                         const bigrational_vector&, const bigrational_vector&,
                         const bool,
                         K_PATCH*&);
  
  friend int get_cyl_side(K_RATPOLY* const,
                          K_RATPOLY* const, K_RATPOLY* const, K_RATPOLY* const,
                          K_RATPOLY* const,
                          K_PATCH*&);
  
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
  
  friend unsigned long gen_partitions(K_PATCH* const, K_PARTITION**&);
  
  friend unsigned long gen_adjacency(K_PARTITION** const, const unsigned long,
                                     K_GRAPH*&,
                                     long*&, K_GRAPH*&);
  
//  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
//                               K_PARTITION**, const unsigned long);
  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
                               K_PARTITION**, const unsigned long,
                               const char);
};

#endif

