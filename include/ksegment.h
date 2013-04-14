#ifndef _KSEGMENT_H
#define _KSEGMENT_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <kratpoly.h>
#include <kpoint2d.h>
#include <kboxco2.h>
#include <kpatch.h>

using namespace std;

class K_POINT2D;
class K_BOXCO2;

class K_SEGMENT
{
friend class K_CURVE;
friend class K_PATCH;
  
  K_POINT2D* start;
  K_POINT2D* end;
  
  unsigned long ref_count;
  
  //  stream
  
  ostream& output(ostream&) const;
  
//public:
  
  //  constructors, assignment and destructor
  
  //  K_SEGMENT()
  //    constructs a null segment.
  
  K_SEGMENT();
  
  //  K_SEGMENT(const K_POINT2D& x, const K_POINT2D& y)
  //    constructs a segment (x, y].
  
  K_SEGMENT(const K_POINT2D&, const K_POINT2D&);
  
  //  K_SEGMENT(K_POINT2D* const x, K_POINT2D* const y)
  //    constructs a segment (x, y].
  
  K_SEGMENT(K_POINT2D* const, K_POINT2D* const);
  
  K_SEGMENT(const K_SEGMENT&);
  K_SEGMENT& operator =(const K_SEGMENT&);
  
  ~K_SEGMENT();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_SEGMENT&);
  
  //  primitives
  
  K_SEGMENT reverse() const;
  
  //  other functions
  
  //  K_BOXCO2 outer_box() const
  //    returns the outer box for *this.
  
  K_BOXCO2  outer_box() const;
  
  //  K_BOXCO2* inner_box() const
  //    returns a pointer to the inner box for *this.
  //    returns 0 if the inner box for *this is not well-def'ed.
  
  K_BOXCO2* inner_box() const;
  
  //  int contains(K_POINT2D& x)
  //    returns 1 if x lies inside the inner box of *this,
  //            0 if x lies outside the outer box of *this.
  //    POSSIBLY DOES NOT TERMINATE!
  //  CAUTION!!!  The start point is OUT of the segment whereas
  //              the end point is IN the segment.
  
  int contains(K_POINT2D& x);
  
  //  gen_curve_topo()
  
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
  
  friend int get_cyl_cap(const bigrational_vector&,
                         const bigrational_vector&, const bigrational_vector&,
                         const bool,
                         K_PATCH*&);
  
  friend int get_cyl_side(K_RATPOLY* const,
                          K_RATPOLY* const,
                          K_RATPOLY* const,
                          K_RATPOLY* const,
                          K_RATPOLY* const,
                          K_PATCH*&);
  
  friend unsigned long get_patches_ell(const bigrational_vector [4],
                                       K_PATCH**&);
  
  friend int get_patch_tor(K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_RATPOLY* const,
                           K_PATCH*&);
};

#endif

