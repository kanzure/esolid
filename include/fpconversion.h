//  file:    fpconversion.h
//  update:  10/14/02

#ifndef _FPCONVERSION_H
#define _FPCONVERSION_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>

using namespace std;

//  int double2bigrational_interval(const double x,
//                                  const unisgned long n,
//                                  bigrational& l, bigrational& r)
//
//  Finds a closed interval `[l, r]' of rational endpoints
//  which contains an IEEE double precision floating-point number `x'
//  and is of width `2^(- n)'.
//
//  The denominators of `l' & `r' are always powers of 2.
//
//  If possible, the denominators of `l' & `r' will be of length
//  `num_bits + 1' bits or fewer.
//
//  Returns 0 if successful,
//          < 0 if unsuccessful (e.g. `x' is infinity, NaN, ...),
//          > 0 if `n' is too strict (see below).
//
//  `n' can be as large as 53 if `x' is normal or
//                         52 if `x' is subnormal.
//  If `n' is too large then the best approximation possible is made.
//  In this case, a positive number is returned as a warning.
//
//  Notes:
//
//  (1) This does not attempt to find the best approximation to `x';
//      it just makes a box of the specified size containing `x'.
//
//  (2) It assumes that doubles are 64 bits, unsigned ints are 32 bits,
//      and probably several other similar assumptions.
//
//  (3) It should run in constant time; there are almost no loops.
//
//  (4) It assumes that `x' is an exact binary number, whose unspecified
//      low-order bits are all 0. That is, `x' is treated as if computed
//      by truncation. It would have been better to assume that it was
//      computed by rounding.
//      example:  given the 3-decimal-digit number 3.14, it produces
//      the output [3.14, 3.15) when asked for 3 digits.
//      Perhaps you'd rather have [3.13, 3.15) since 3.14 may have been
//      rounded up from 3.136.

int double2bigrational_interval(const double,
                                const unsigned int,
                                bigrational&, bigrational&);

bigrational as_bigrational(const float);

#endif

