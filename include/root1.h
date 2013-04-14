#ifndef _ROOT1_H
#define _ROOT1_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <kratpoly.h>

using namespace std;

class ROOT1
{
friend class K_POINT1D;
friend class K_POINT2D;

  K_RATPOLY*      poly;           //  a polynomial whose roots are *this

  bigrational     low;            //  an lower bound on the root
  bigrational     high;           //  an upper bound on the root

  long            num_perm_low;   //  the number of Sturm seq. permanencies
  long            num_perm_high;  //    at low and high

  unsigned long   num_roots;      //  the number of roots in the interval

  int             sig_low;        //  Set to be - 1 if poly(low) < 0,
                                  //              1 if poly(low) > 0,
                                  //  BUT         0 if
                                  //    num_roots > 1 or
                                  //    poly(low) * poly(high) > 0 or
                                  //    poly(low) * poly(high) == 0.

  double          float_est;      //  a floating-point estimate of the root
  int             ok_float;       //  Set ot be 1 if the estimate is OK and
                                  //            0 otherwise.

  unsigned long   ref_count;      //  reference counter

//  PseudorootList*	pseudo_roots;   //  Pseudoroots.

  //  stream

  ostream& output(ostream&) const;

  //  other

  //  int cut(const bigrational& b)
  //    cuts [low, high] at b s.t. [low, hgih] becomes
  //      either [low, b] or [b, high] that contains some roots of *poly.

  int cut(const bigrational&);

  //  int halve()
  //    cuts [low, high] at half = (low + high) / 2 s.t. [low, high] becomes
  //      either [low, half] or [half, high] that contains some roots of *poly.

  int halve();

  //  int reduce(const unsigned long num_bits)
  //    reduces *this s.t. [low, high] contains all the numbers
  //      approximated  by float_est to num_bits precision.
  //    returns 1 if *this is reduced and
  //            0 otherwise.
  //    works only when float_est has been computed.
  //    returns 0 if the float_est is not good.

  int reduce(const unsigned long);

  //  int contract(const bigrational& tol)
  //    contracts *this s.t. [low, high] is no larger than tol.

  int contract(const bigrational&);

  //  int shrink(const bigrational& fac)
  //    shrinks *this s.t. [low, high] becomes smaller by fac.
  //    e.g. fac == 1/10 => [low, high] becomes 1/10.

  int shrink(const bigrational&);

public:

  //  constructors, assignment and destructor

  ROOT1();

//  ROOT1(const K_RATPOLY&, const bigrational&, const bigrational&);

  //  ROOT1(const K_RATPOLY& P,
  //        const bigrational& l, const bigrational& h,
  //        const long n_l, const long n_h)
  //    constructs an instance for the set of the roots of P in [l, h].
  //    the num. of Sturm seq. permanencies at l and h are n_l and n_h
  //      unless n_l and n_h are - 1, resp.

  ROOT1(const K_RATPOLY&, const bigrational&, const bigrational&,
        const long = - 1, const long = - 1);

  //  ROOT1(const K_RATPOLY& P,
  //        const bigrational& l, const bigrational& h,
  //        const double d)
  //    constructs an instance for the set of the roots of P in [l, h].
  //    the floating-point estimate of the root is d if num_roots == 1.

  ROOT1(const K_RATPOLY&, const bigrational&, const bigrational&,
        const double);

  ROOT1(const ROOT1&);
  ROOT1& operator =(const ROOT1&);

  ~ROOT1();

  //  stream

  friend ostream& operator <<(ostream&, const ROOT1&);

  //  other

  //  unsigned long isolate_roots(ROOT1*& R, const bigrational& tol)
  //    isolates roots of *poly in [low, high] and store them in R.
  //    each R[i] contains 1 and only 1 root of *poly and
  //    (R[i].low, R[i].high) is no larger than tol.

  unsigned long isolate_roots(ROOT1*&, const bigrational& tol);

  //  get_pts(): get algebraic points

  friend unsigned long get_pts(const bigrational&, const bigrational&,
                               const K_RATPOLY&,
                               K_POINT1D**&,
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
};

#endif

