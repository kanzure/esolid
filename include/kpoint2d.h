#ifndef _KPOINT2D_H
#define _KPOINT2D_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <kratpoly.h>
#include <root1.h>
#include <kpoint1d.h>
#include <kboxco2.h>

using namespace std;

class K_BOXCO2;
class K_PARTITION;
class K_SURF;
class K_PATCH;

class K_POINT2D
{
friend class PTS_CACHE;
friend class K_SEGMENT;
friend class K_CURVE;
friend class PT_SURF_ASSOC;
friend class K_PATCH;
friend class K_PARTITION;
friend class K_SOLID;

  mutable unsigned int type;  //  type == 1 if ROOT1 in s and t
                              //          2 if ROOT1 in s and bigrational in t
                              //          3 if bigrational in s and ROOT1 in t
                              //          4 if bigrational in s and t

  mutable ROOT1*      PtRs;  //  the interval for the s-coordinate of *this,
                             //  0 if type == 3 or 4
  mutable ROOT1*      PtRt;  //  the interval for the t-coordinate of *this,
                             //  0 if type == 2 or 4
  mutable bigrational PtBs;  //  the point for the s-coordinate of *this,
                             //  undefined if type == 1 or 2
  mutable bigrational PtBt;  //  the point for the t-coordinate,
                             //  undefined if type == 1 or 3.

  K_RATPOLY* poly1;  //  polynomials
  K_RATPOLY* poly2;  //    that intersect with each other at *this

  mutable K_POINT2D* prev;  //  the pointers to the neighbors
  mutable K_POINT2D* next;  //    in the linked list of equivalent pts

  K_POINT2D* pt_in_other_dom;  //  the pointer to
                               //    the same pt in the other domain

  unsigned long ref_count;  //  reference counter

  //  constructors, assignment and destructor

  //  K_POINT2D(const bigrational& b_s, const K_POINT1D& x_t,
  //            const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    construct a type 3 or 4 instance from
  //      bigrational b_s, K_POINT1D x_t and K_RATPOLY P1 and P2.
  //    P1 must intersect with P2 at (b_s, x_t).

  K_POINT2D(const bigrational&, const K_POINT1D&,
            const K_RATPOLY&, const K_RATPOLY&);

  //  K_POINT2D(const K_POINT1D& x_s, const bigrational& b_t,
  //            const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    construct a type 2 or 4 instance from
  //      K_POINT1D x_s, bigrational b_t and K_RATPOLY P1 and P2.
  //    P1 must intersect with P2 at (x_s, b_t).

  K_POINT2D(const K_POINT1D&, const bigrational&,
            const K_RATPOLY&, const K_RATPOLY&);

  //  K_POINT2D(const K_POINT1D& x_s, const K_POINT1D& x_t,
  //            const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    construct a type 1 or 2 or 3 or 4 instance from
  //      K_POINT1D x_s and x_t and K_RATPOLY P1 and P2.
  //    P1 must intersect with P2 at (x_s, x_t).

  K_POINT2D(const K_POINT1D&, const K_POINT1D&,
            const K_RATPOLY&, const K_RATPOLY&);

  //  stream

  ostream& output(ostream&) const;

  //  primitive

  //  int cut_s(const bigrational& b_s) const
  //    cut *this by the line s = b_s, i.e.,
  //    refine the box for *this
  //      by setting get_low_s() or get_high_s() to be b_s.

  int cut_s(const bigrational&) const;

  //  int cut_t(const bigrational& b_t) const
  //    cut *this by the line t = b_t, i.e.,
  //    refine the box for *this
  //      by setting get_low_t() or get_high_t() to be b_t.

  int cut_t(const bigrational&) const;

  //  int halve_s() const
  //    cut *this by the line s = half that halves the box for *this, i.e.,
  //    refine the box for *this
  //      by setting get_low_s() or get_high_s() to be half.

  int halve_s() const;

  //  int halve_t() const
  //    cut *this by the line t = half that halves the box for *this, i.e.,
  //    refine the box for *this
  //      by setting get_low_t() or get_high_t() to be half.

  int halve_t() const;

  //  int reduce_s(const unsigned long num_bits) const
  //    reduce *this s.t.
  //      [get_low_s(), get_high_s()] will contain all the numbers
  //        approximated by PtRs->float_est to num_bits precision.
  //    return 1 if some reduction occurs and
  //           0 otherwise.

  int reduce_s(const unsigned long) const;

  //  int reduce_t(const unsigned long num_bits) const
  //    reduce *this s.t.
  //      [get_low_t(), get_high_t()] will contain all the numbers
  //        approximated by PtRt->float_est to num_bits precision.
  //    return 1 if some reduction occurs and
  //           0 otherwise.

  int reduce_t(const unsigned long) const;

  //  int reduce(const unsigned long num_bits_s,
  //             const unsigned long num_bits_t) const
  //    reduce *this s.t.
  //      [get_low_s(), get_high_s()] will contain all the numbers
  //        approximated by PtRs->float_est to num_bits_s precision and
  //      [get_low_t(), get_high_t()] will contain all the numbers
  //        approximated by PtRt->float_est to num_bits_t precision.
  //    return 1 if some reduction occurs and
  //           0 otherwise.

  int reduce(const unsigned long, const unsigned long) const;

  //  int contract_s(const bigrational& tol) const
  //    contract *this s.t.
  //      [get_low_s(), get_high_s()] will be no larger than tol.

  int contract_s(const bigrational&) const;

  //  int contract_t(const bigrational& tol) const
  //    contract *this s.t.
  //      [get_low_t(), get_high_t()] will be no larger than tol.

  int contract_t(const bigrational&) const;

  //  int contract(const bigrational& tol_s, const bigrational& tol_t) const
  //    contract *this s.t.
  //      [get_low_s(), get_high_s()] will be no larger than tol_s and
  //      [get_low_t(), get_high_t()] will be no larger than tol_t.

  int contract(const bigrational&, const bigrational&) const;

  //  int shrink(const bigrational& fac_s, const bigrational& fac_t) const
  //    shrink *this s.t.
  //      [get_low_s(), get_high_s()] will become smaller by fac_s and
  //      [get_low_t(), get_high_t()] will become smaller by fac_t.

  int shrink(const bigrational&, const bigrational&) const;

  //  int K_POINT2D :: separate(const K_POINT2D& x) const
  //    separate *this and x s.t. they will not overlap.
  //    POSSIBLY DOES NOT TERMINATE.

  int separate(const K_POINT2D&) const;

  //  K_BOXCO2 bbox() const
  //    return a bounding box for *this.

  K_BOXCO2 bbox() const;

public:

  //  constructors, assignment and destructor

  K_POINT2D();

  //  K_POINT2D(const bigrational& b_s, const bigrational& b_t)
  //    construct a type 4 instance from bigrational b_s and b_t.

  K_POINT2D(const bigrational&, const bigrational&);

  //  K_POINT2D(const ROOT1& r_s, const ROOT1& r_t,
  //            const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    construct a type 2 instance from
  //      ROOT1 r_s and r_t and K_RATPOLY P1 and P2.
  //    r_s.num_roots and r_t.num_roots must have been 1 and
  //    P1 must intersect with P2 at (r_s, r_t).

  K_POINT2D(const ROOT1&, const ROOT1&,
            const K_RATPOLY&, const K_RATPOLY&);

  //  K_POINT2D(const ROOT1& r_s, const bigrational& b_t,
  //            const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    construct a type 2 instance from
  //      bigrational b_s, ROOT1 r_t and K_RATPOLY P1 and P2.
  //    r_s.num_roots must have been 1 and
  //    P1 must intersect with P2 at (r_s, b_t).

  K_POINT2D(const ROOT1&, const bigrational&,
            const K_RATPOLY&, const K_RATPOLY&);

  //  K_POINT2D(const bigrational& b_s, const ROOT1& r_t,
  //            const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    construct a type 3 instance from
  //      bigrational b_s, ROOT1 r_t and K_RATPOLY P1 and P2.
  //    r_t.num_roots must have been 1 and
  //    P1 must intersect with P2 at (b_s, r_t).

  K_POINT2D(const bigrational&, const ROOT1&,
            const K_RATPOLY&, const K_RATPOLY&);

  //  K_POINT2D(const bigrational& b_s, const bigrational& b_t,
  //            const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    construct a type 4 instance from
  //      bigrational b_s and b_t and K_RATPOLY P1 and P2.
  //    P1 must intersect with P2 at (b_s, b_t).

  K_POINT2D(const bigrational&, const bigrational&,
            const K_RATPOLY&, const K_RATPOLY&);

	K_POINT2D(const K_POINT2D&);
  K_POINT2D& operator =(const K_POINT2D&);

  ~K_POINT2D();

  //  stream

  friend ostream& operator <<(ostream&, const K_POINT2D&);

  // primitive

  bigrational get_low_s() const;
  bigrational get_high_s() const;
  bigrational get_low_t() const;
  bigrational get_high_t() const;

  bigrational_vector get_low() const;
  bigrational_vector get_high() const;

  //  comparison

  //  int compare_s(const K_POINT2D& x) const
  //    return   1 if *this >  x,
  //             0 if *this == x, and
  //           - 1 if *this <  x  in the s-coordinate.

  int compare_s(const K_POINT2D&) const;

  //  int compare_t(const K_POINT2D& x) const
  //    return   1 if *this >  x,
  //             0 if *this == x, and
  //           - 1 if *this <  x  in the t-coordinate.

  int compare_t(const K_POINT2D&) const;

  //  int compare(const K_POINT2D& x, const unsigned long i) const
  //    return   1 if *this >  x,
  //             0 if *this == x, and
  //           - 1 if *this <  x  in the i-th coordinate.

//  int compare(const K_POINT2D&, const unsigned long) const;
//
//  //  int compare_s(const bigrational& b_s) const
//  //    return   1 if *this >  b_s,
//  //             0 if *this == b_s, and
//  //           - 1 if *this <  b_s  in the s-coordinate.

  int compare_s(const bigrational&) const;

  //  int compare_t(const bigrational& b_t) const
  //    return   1 if *this >  b_t,
  //             0 if *this == b_t, and
  //           - 1 if *this <  b_t  in the t-coordinate.

  int compare_t(const bigrational&) const;

//  //  int compare(const bigrational& b, const unsigned long i) const
//  //    return   1 if *this >  b,
//  //             0 if *this == b, and
//  //           - 1 if *this <  b  in the i-th coordinate.
//
//  int compare(const bigrational&, const unsigned long) const;

  //  int sort_s(K_POINT2D** const X, const unsigned long n)
  //    insertion-sort the array X of length n in the s-coordinate.
  //    return 1 if all the elements are distinct in the s-coordinate and
  //           0 otherwise.

  friend int sort_s(K_POINT2D** const, const unsigned long);

  //  int sort_t(K_POINT2D** const X, const unsigned long n)
  //    insertion-sort the array X of length n in the t-coordinate.
  //    return 1 if all the elements are distinct in the t-coordinate and
  //           0 otherwise.

  friend int sort_t(K_POINT2D** const, const unsigned long);

//  //  int sort(K_POINT2D** const X, const unsigned long n,
//  //           const unsigned long i)
//  //    insertion-sort the array X of length n in the i-th coordinate.
//  //    return 1 if all the elements are distinct in the i-th coordinate and
//  //           0 otherwise.
//
//  friend int sort(K_POINT2D** const, const unsigned long, const unsigned long);

  //  int overlap(const K_POINT2D& x) const
  //    return 1 if *this and x overlap, and
  //           0 otherwise.

  int overlap(const K_POINT2D&) const;

  //  int equiv(const K_POINT2D& x) const
  //    return 1 if *this and x are identical or
  //                            have already been identified, and
  //            0 otherwise

  int equiv(const K_POINT2D&) const;

  //  int equal(K_POINT2D& x)
  //    return 1 if the boxes for *this and x are equal, and
  //           0 otherwise

  int equal(K_POINT2D&);

  //  int equal_s(K_POINT2D& x)
  //    returns 1 if *this and x are equal in the s-coordinate, and
  //            0 otherwise.

  int equal_s(K_POINT2D& x);

  //  int equal_t(K_POINT2D& x)
  //    returns 1 if *this and x are equal in the t-coordinate, and
  //            0 otherwise.

  int equal_t(K_POINT2D& x);

//  //  int equal_dir(K_POINT2D& x, const unsigned long i)
//  //    returns 1 if *this and x are equal in the i-th coordinate, and
//  //            0 otherwise.
//
//  int equal_dir(K_POINT2D&, const unsigned long);

  //  get_pts(): get algebraic points

  //  unsigned long get_pts(const bigrational& l_s, const bigrational& h_s,
  //                        const bigrational& l_t, const bigrational& h_t,
  //                        const K_RATPOLY& P1, const K_RATPOLY& P2,
  //                        K_POINT2D**& pts,
  //                        const bigrational& tol, const int count_edges)
  //    computes the intersections "pts" of 2 bivariate polynomials P1 and P2
  //                              on/in the region [l_s, h_s] x [l_t, h_t].
  //    returns the number of intersections.
  //    if tol > 0 then
  //      a box for each intersection is no larger than tol in any direction.
  //    if tol = 0 then
  //      there is no limit on the size of boxes for intersections.
  //    if count_edges = 1 then
  //      the intersections on the boundary edges of the region are counted.
  //    if count_edges = 0 then
  //      the intersections on the boundary edges of the region are ignored.

  friend unsigned long get_pts(const bigrational&, const bigrational&,
                               const bigrational&, const bigrational&,
                               const K_RATPOLY&, const K_RATPOLY&,
                               K_POINT2D**&,
                               const bigrational&, const int);

  friend unsigned long get_pts_interior(const bigrational&, const bigrational&,
                                        const bigrational&, const bigrational&,
                                        const K_RATPOLY&, const K_RATPOLY&,
                                        K_POINT2D**&,
                                        const bigrational&);

  friend int refine_interior(K_POINT1D* const, K_POINT1D* const,
                             const K_RATPOLY&, const K_RATPOLY&,
                             K_POINT2D*&,
                             const bigrational&);

  friend unsigned long get_pts_proto(const bigrational&, const bigrational&,
                                     const K_RATPOLY&,
                                     const bigrational&,
                                     const K_RATPOLY&, const K_RATPOLY&,
                                     K_POINT2D**&,
                                     const bigrational&, const int);

  friend unsigned long get_pts_proto(const bigrational&,
                                     const bigrational&, const bigrational&,
                                     const K_RATPOLY&,
                                     const K_RATPOLY&, const K_RATPOLY&,
                                     K_POINT2D**&,
                                     const bigrational&, const int);

  //  unsigned long get_all_pts(const K_RATPOLY& P1, const K_RATPOLY& P2,
  //                            K_POINT2D**& pts,
  //                            const bigrational& tol)
  //    computes all the intersections "pts" of
  //      2 bivariate polynomials P1 and P2.
  //    returns the number of intersections.
  //    if tol >  0 then
  //      a box for each intersection is no larger than tol in any direction.
  //    if tol = 0 then
  //      there is no limit on the size of boxes for intersections.

  friend unsigned long get_all_pts(const K_RATPOLY&, const K_RATPOLY&,
                                   K_POINT2D**&,
                                   const bigrational&);

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

  //  other

  friend int match_pts(K_POINT2D** const, const unsigned long, K_SURF* const,
                       K_POINT2D** const, const unsigned long, K_SURF* const);

  friend int pt_inside_trim_curves(K_POINT2D&,
                                   K_CURVE** const, const unsigned long);

  friend int pt_in_on_out_trim_curves(K_POINT2D&,
                                      K_CURVE** const, const unsigned long);

  friend unsigned long gen_partitions(K_PATCH* const, K_PARTITION**&);

  int get_fp_approx(double*&) const;
};

#endif

