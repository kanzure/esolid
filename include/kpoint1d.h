#ifndef _KPOINT1D_H
#define _KPOINT1D_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <root1.h>
#include <kratpoly.h>

using namespace std;

class K_POINT1D
{
friend class K_POINT2D;
  
  mutable unsigned int type;  //  type == 1 if ROOT1
                              //          2 if bigrational
  
  mutable ROOT1*      PtR;  //  the interval for *this, 0 if type == 2
  mutable bigrational PtB;  //  the point for *this, undefined if type == 1
  
  K_RATPOLY*  poly;  //  a polynomial that vanished at *this
  
  unsigned long ref_count;  //  reference counter
  
  //  stream
  
  ostream& output(ostream&) const;
  
  //  primitive
  
  //  int cut(const bigrational& b) const
  //    cut *this at the point b, i.e.,
  //    refine the interval for *this
  //      by setting get_low() or get_high() to be b.
  
  int cut(const bigrational&) const;
  
  //  int halve() const
  //    cut *this at the point half that halves the interval for *this, i.e.,
  //    refine the interval for *this
  //      by setting get_low() or get_high() to be half.
  
  int halve() const;
  
  //  int reduce(const unsigned long num_bits) const
  //    reduce *this s.t.
  //      [get_low(), get_high()] will contain all the numbers
  //        approximated by float_est to num_bits precision.
  //    return 1 if some reduction occurs and
  //           0 otherwise.
  
  int reduce(const unsigned long) const;
  
  //  int contract(const bigrational& tol) const
  //    contract *this s.t. [get_low(), get_high()] will be no larger than tol.
  
  int contract(const bigrational&) const;
  
  //  int shrink(const bigrational& fac) const
  //    shrink *this s.t. [get_low(), get_high()] will become smaller by fac.
  //    e.g. fac == 1/10 => [get_low(), get_high()] becomes 1/10.
  
  int shrink(const bigrational&) const;
  
  //  int K_POINT1D :: separate(const K_POINT1D& x) const
  //    separate *this and x s.t. they will not overlap.
  //    POSSIBLY DOES NOT TERMINATE.
  
  int separate(const K_POINT1D&) const;
  
  //  arithmetic
  
  //  K_POINT1D add(const K_POINT1D& x) const
  //    return *this + x
  
  K_POINT1D add(const K_POINT1D&) const;
  
  //  K_POINT1D sub(const K_POINT1D& x) const
  //    return *this - x
  
  K_POINT1D sub(const K_POINT1D&) const;
  
  //  K_POINT1D mul(const K_POINT1D& x) const
  //    return *this * x
  
  K_POINT1D mul(const K_POINT1D&) const;
  
  //  K_POINT1D div(const K_POINT1D& x) const
  //    return *this / x
  
  K_POINT1D div(const K_POINT1D&) const;
  
public:
  
  //  constructors, assignment and destructor
  
  //  K_POINT1D()
  //    constructs a type 0 instance.
  
  K_POINT1D();
  
  //  K_POINT1D(const ROOT1& r)
  //    construct a type1 instance from ROOT1 r.
  //    r.num_root must have been 1.
  
  K_POINT1D(const ROOT1&);
  
  //  K_POINT1D(const ROOT1& r)
  //    construct a type1 instance from ROOT1 r of a root of K_RATPOLY P.
  //    r.num_roots must have been 1.
  
  K_POINT1D(const ROOT1&, const K_RATPOLY&);
  
  //  K_POINT1D(const bigrational& b)
  //    construct a type2 instance from bigrational b.
  
  K_POINT1D(const bigrational&);
  
  //  K_POINT1D(const bigrational& b, const K_RATPOLY& P)
  //    construct a type2 instance from bigrational b and K_RATPOLY P.
  //    P must vanish at b.
  
  K_POINT1D(const bigrational&, const K_RATPOLY&);
  
  K_POINT1D(const K_POINT1D&);
  K_POINT1D& operator =(const K_POINT1D&);
  
  ~K_POINT1D();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_POINT1D&);
  
  //  primitive
  
  bigrational get_low() const;
  bigrational get_high() const;
  
  //  arithmetic
  
  friend K_POINT1D operator +(const K_POINT1D&, const K_POINT1D&);
  friend K_POINT1D operator -(const K_POINT1D&, const K_POINT1D&);
  friend K_POINT1D operator *(const K_POINT1D&, const K_POINT1D&);
  friend K_POINT1D operator /(const K_POINT1D&, const K_POINT1D&);
  
  //  comparison
  
  //  int compare(const K_POINT1D& x) const
  //    return   1 if *this > x,
  //             0 if *this == x, and
  //           - 1 if *this < x.
  
  int compare(const K_POINT1D&) const;
  
  //  int compare(const bigrational& b) const
  //    return   1 if *this > b,
  //             0 if *this == b, and
  //           - 1 if *this < b.
  
  int compare(const bigrational&) const;
  
  //  int sort(K_POINT1D** const X, const unsigned long n)
  //    perform insertion sort on the array X of length n.
  //    return 1 if elements are distinct and
  //           0 otherwise.
  
  friend int sort(K_POINT1D** const, const unsigned long);
  
  //  int overlap(const K_POINT1D& x) const
  //    return 1 if *this and x overlap, and
  //           0 otherwise.
  
  int overlap(const K_POINT1D&) const;
  
  //  get_pts(): get algebraic points
  
  //  unsigned long get_pts(const bigrational& l, const bigrational& h,
  //                        const K_RATPOLY& P,
  //                        K_POINT1D**& pts,
  //                        const bigrational& tol, const int count_endpts)
  //    computes the roots "pts" of a univariate polynomial P
  //                      on/in the interval [l, h].
  //    returns the number of roots.
  //    if tol > 0 then an interval for each root is no larger than tol.
  //    if tol = 0 then there is no limit on the size of intervals for roots.
  //    if count_endpts = 1 then
  //      the roots at the endpoints of the interval [l, h] are counted.
  //    if count_endpts = 0 then
  //      the roots on the endpoints of the interval [l, h] are ignored.
  
  friend unsigned long get_pts(const bigrational&, const bigrational&,
                               const K_RATPOLY&,
                               K_POINT1D**&,
                               const bigrational&, const int);
  
  //  unsigned long get_all_pts(const K_RATPOLY& P,
  //                            K_POINT1D**& pts,
  //                            const bigrational& tol)
  //    computes all the intersections "pts" of a univariate polynomials P.
  //    returns the number of roots.
  //    if tol > 0 then an interval for each root is no larger than tol.
  //    if tol = 0 then there is no limit on the size of intervals for roots.
  
  friend unsigned long get_all_pts(const K_RATPOLY&,
                                   K_POINT1D**&,
                                   const bigrational);
  
  //  other
  
  friend unsigned long match_intervals(const bigrational* const,
                                       const bigrational* const,
                                       const unsigned long,
                                       K_POINT1D** const,
                                       long*&,
                                       const unsigned long);
  
  friend unsigned long get_pts_interior(const bigrational&, const bigrational&,
                                        const bigrational&, const bigrational&,
                                        const K_RATPOLY&, const K_RATPOLY&,
                                        K_POINT2D**&,
                                        const bigrational&);
  
  friend int refine_interior(K_POINT1D* const, K_POINT1D* const,
                             const K_RATPOLY&, const K_RATPOLY&,
                             K_POINT2D*&,
                             const bigrational&);
};

#endif

