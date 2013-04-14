#ifndef _FPFILTER_H
#define _FPFILTER_H

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iostream>

#define RAT_TO_DBL_ERR 4 * DBL_EPSILON

using namespace std;

int round_nearest();
int round_down();
int round_up();
int round_zero();
//  Set rounding mode.

int
horner_with_err(
  const double* const poly,        //  coefficients of a univariate polynomial
                                   //  in decreasing degree order
  const unsigned long deg,         //  degree of poly
  const double        coeff_err,   //  error bound for coeffs
  const double        in_val,      //  evaluate poly here
  const double        in_err,      //  error bound for in_val
  double&             out_val,     //  value of poly
  double&             out_err      //  bound on absolute error
  );
//  Evaluate poly at in_val. Also, compute a bound on error.

int sign_given_err(const double val, const double err);
//  Returns the sign of val if val might be off from 0.0 by as much as err.
//  Returns 0 if we can't tell.

double newton_for_pseudoroot(const double* const poly,
                             const unsigned long deg,
                             const double coeff_err,
                             const double* const deriv,
                             const double low,
                             const double high);

#endif

