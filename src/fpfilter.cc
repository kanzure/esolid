//  file:    fpfilter.cc
//  update:  08/20/03

#include <config.h>

#include <fpfilter.h>

#if (defined(_Linux_i386_))
#include <fpu_control.h>
#elif ((defined(_SunOS_)) || (defined(_IRIX_)))
#include <ieeefp.h>
#endif

//#include <float.h>
//  Are u using Microsoft?!? Sorry I have no idea what a mess they do!

int round_nearest()
{
//#ifdef _Linux_i386_
  
  fpu_control_t w;
  
  _FPU_GETCW(w);
  _FPU_SETCW(w & 0xf3ff | _FPU_RC_NEAREST);
  
//#endif
//#ifdef _SunOS_
//  
//  fp_rnd m = FP_RN;
//  
//  fpsetround(m);
//  
//#endif
//#ifdef _IRIX_
//  
//  fp_rnd m = FP_RN;
//  
//  fpsetround(m);
//  
//#endif
  
  return 0;
}

int round_down()
{
//#ifdef _Linux_i386_
  
  fpu_control_t w;
  
  _FPU_GETCW(w);
  _FPU_SETCW(w & 0xf3ff | _FPU_RC_DOWN);
  
//#endif
//#ifdef _SunOS_
//  
//  fp_rnd m = FP_RM;
//  
//  fpsetround(m);
//  
//#endif
//#ifdef _IRIX_
//  
//  fp_rnd m = FP_RM;
//  
//  fpsetround(m);
//  
//#endif
  
  return 0;
}

int round_up()
{
//#ifdef _Linux_i386_
  
  fpu_control_t w;
  
  _FPU_GETCW(w);
  _FPU_SETCW(w & 0xf3ff | _FPU_RC_UP);
  
//#endif
//#ifdef _SunOS_
//  
//  fp_rnd m = FP_RP;
//  
//  fpsetround(m);
//  
//#endif
//#ifdef _IRIX_
//  
//  fp_rnd m = FP_RP;
//  
//  fpsetround(m);
//  
//#endif
  
  return 0;
}

int round_zero()
{
//#ifdef _Linux_i386_
  
  fpu_control_t w;
  
  _FPU_GETCW(w);
  _FPU_SETCW(w & 0xf3ff | _FPU_RC_ZERO);
  
//#endif
//#ifdef _SunOS_
//  
//  fp_rnd m = FP_RZ;
//  
//  fpsetround(m);
//  
//#endif
//#ifdef _IRIX_
//  
//  fp_rnd m = FP_RZ;
//  
//  fpsetround(m);
//  
//#endif
  
  return 0;
}

int horner_with_err(const double* const poly,
                    const unsigned long deg,
                    const double        coeff_err,
                    const double        in_val,
                    const double        in_err,
                    double&             out_val,
                    double&             out_err)
{
  unsigned long i;
  int           f;
  double        abs_in_val;
  double        rel_err;
  
  f = finite(in_val);
  
  //  1. Evaluate poly at in_val.
  
  if (f)
    f = finite(out_val = poly[0]);
  
  for (i = 1; f && i <= deg; i++)
    if (f = finite(out_val *= in_val))
      f = finite(out_val += poly[i]);
  
  //  2. Compute a bound on error.  We round towards +infinity
  //     to make sure the bound is sufficient. Since the errors
  //     on each term could conspire in this way, the error is
  //     bounded by
  //       (max per-term error) * (value of abs-val-coeff poly).
  
  //  2-1. Set MPU mode.
  
  round_up();
  
  //  2-2. Evaluate the polynomial whose coeffients are
  //       the absolute values of the coeffients of poly.
  
  if (f)
    if (f = finite(abs_in_val = fabs(in_val)))
      f = finite(out_err = fabs(poly[0]));
  
  for (i = 1; f && i <= deg; i++)
    if (f = finite(out_err *= abs_in_val))
      f = finite(out_err += fabs(poly[i]));
  
  //  2-3. Compute the maximum relative error for any single term.
  
  if (f)
    f = finite(rel_err = deg * (2 * DBL_EPSILON + in_err) + coeff_err);
  
  //  2-4. Compute a bound on error.
  
  if (f)
    f = finite(out_err *= rel_err);
  
  //  2-5. Retrieve the original rounding mode.
  
  round_nearest();
  
  return f;
}

//  int sign_given_err(double val, double err)
//    returns some positive value if val is certainly positive,
//            some negative value if val is certainly negative, and
//            0 otherwise.

int sign_given_err(const double val, const double err)
{
  assert(finite(val));
  assert(finite(err));
  
  double v_l, v_h;
  
  if (finite(v_l = val - err) && val > 0.0 && v_l > 0.0)
    return 1;
  else if (finite(v_h = val + err) && val < 0.0 && v_h < 0.0)
    return - 1;
  else
    return 0;
}

double newton_for_pseudoroot(const double* const poly,
                             const unsigned long deg,
                             const double coeff_err,
                             const double* const deriv,
                             const double l0,
                             const double h0)
{
  double l, h, l_val, l_err, h_val, h_err, m, p_val, p_err, d_val, d_err;
  int    l_sign, h_sign, p_sign, d_sign;
  
  l = l0;
  horner_with_err(poly, deg, coeff_err, l, 0.0, l_val, l_err);
  l_sign = sign_given_err(l_val, l_err);
  h = h0;
  horner_with_err(poly, deg, coeff_err, h, 0.0, h_val, h_err);
  h_sign = sign_given_err(h_val, h_err);
  assert(l_sign * h_sign < 0);
  
  m = (l + h) / 2.0;
  
  while (1)
  {
    horner_with_err(poly, deg, coeff_err, m, 0.0, p_val, p_err);
    p_sign = sign_given_err(p_val, p_err);
    
    if (p_sign == 0 || h - l < DBL_EPSILON * h)
      return m;
    else if (p_sign == l_sign)
      l = m;
    else
      h = m;
    
    horner_with_err(deriv, deg - 1, coeff_err + DBL_EPSILON,
                    m, 0.0,
                    d_val, d_err);
    m -= p_val / d_val;
    
    if (m <= l || m >= h)
      m = (l + h) / 2.0;
  }
  
  return 0.0;
}

