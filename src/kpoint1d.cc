#include <config.h>
#ifdef _EXPERIMENT
#include <counter.h>
#endif

#include <kpoint1d.h>

#include <fpconversion.h>

//  K_POINT1D :: K_POINT1D()
//    constructs a type 0 instance.

K_POINT1D :: K_POINT1D()
{
  type = 0;

  PtR = 0;
  PtB = 0;

  poly = 0;

  ref_count = 0;
}

//  K_POINT1D :: K_POINT1D(const ROOT1& r)
//    construct a type1 instance from ROOT1 r.
//    r.num_roots must have been 1.

K_POINT1D :: K_POINT1D(const ROOT1& r)
{
  assert(r.num_roots == 1);

  type = 1;

  PtR = new ROOT1(r);
  PtR->ref_count++;

  if (poly = PtR->poly)
    poly->ref_count++;

  ref_count = 0;
}

//  K_POINT1D :: K_POINT1D(const ROOT1& r)
//    construct a type1 instance from ROOT1 r of a root of K_RATPOLY P.
//    r.num_roots must have been 1.

K_POINT1D :: K_POINT1D(const ROOT1& r, const K_RATPOLY& P)
{
  assert(r.num_roots == 1);
//  assert(P.num_vars == 1);

//  K_RATPOLY R = rem(P, *r.poly);
//
//  assert(R.deg[0] == 0 && R.coeffs[0] == 0);

  type = 1;

  PtR = new ROOT1(r);
  PtR->ref_count++;

  poly = new K_RATPOLY(P);
  poly->ref_count++;

  ref_count = 0;
}

//  K_POINT1D :: K_POINT1D(const bigrational& b)
//    construct a type2 instance from bigrational b.

K_POINT1D :: K_POINT1D(const bigrational& b)
{
  type = 2;

  PtR = 0;
  PtB = b;

  poly = new K_RATPOLY(1, 0, PtB);
  poly->ref_count++;

  ref_count = 0;
}

//  K_POINT1D :: K_POINT1D(const bigrational& b, const K_RATPOLY& P)
//    construct a type2 instance from bigrational b and K_RATPOLY P.
//    P must vanish at b.

K_POINT1D :: K_POINT1D(const bigrational& b, const K_RATPOLY& P)
{
//  assert(P.num_vars == 1);
//  assert(!P.sgn_at(b));

  type = 2;

  PtR = 0;
  PtB = b;

  poly = new K_RATPOLY(P);
  poly->ref_count++;

  ref_count = 0;
}

//  K_POINT1D :: K_POINT1D(const K_POINT1D& x)
//    the copy constructor

K_POINT1D :: K_POINT1D(const K_POINT1D& x)
{
  type = x.type;

  if (PtR = x.PtR)
    PtR->ref_count++;

  PtB = x.PtB;

  if (poly = x.poly)
    poly->ref_count++;

  ref_count = 0;
}

//  K_POINT1D& K_POINT1D :: operator =(const K_POINT1D& x)
//    the assignment operator

K_POINT1D& K_POINT1D :: operator =(const K_POINT1D& x)
{
  if (this != &x)
  {
    if (PtR && !--PtR->ref_count)
      delete PtR;

    if (poly && !--poly->ref_count)
      delete poly;

    type = x.type;

    if (PtR = x.PtR)
      PtR->ref_count++;

    PtB = x.PtB;

    if (poly = x.poly)
      poly->ref_count++;
  }

  return *this;
}

//  K_POINT1D :: ~K_POINT1D()
//    the destructor

K_POINT1D :: ~K_POINT1D()
{
  if (PtR && !--PtR->ref_count)
    delete PtR;

  if (poly && !--poly->ref_count)
    delete poly;
}

ostream& K_POINT1D :: output(ostream& o) const
{
  if (type == 1)
    o << *PtR << flush;
  else if (type == 2)
    o << PtB << flush;
  else  //  if type == 0
    o << " NULL " << flush;

  return o;
}

ostream& operator <<(ostream& o, const K_POINT1D& x)
{
  return x.output(o);
}

bigrational K_POINT1D :: get_low() const
{
  assert(type);

  bigrational l;

  if (type == 1)
    l = PtR->low;
  else  //  if (type == 2)
    l = PtB;

  return l;
}

bigrational K_POINT1D :: get_high() const
{
  assert(type);

  bigrational h;

  if (type == 1)
    h = PtR->high;
  else  //  if (type == 2)
    h = PtB;

  return h;
}

//  int K_POINT1D :: cut(const bigrational& b) const
//    cuts *this at the point b, i.e.,
//    refines the interval for *this
//      by setting get_low() or get_high() to be b.

int K_POINT1D :: cut(const bigrational& b) const
{
  assert(type);

  if (type == 1 && b > PtR->low && b < PtR->high)
    if (!PtR->poly->sgn_at(b))
    {
      type = 2;

      if (!--PtR->ref_count)
        delete PtR;

      PtR = 0;
      PtB = b;
    }
    else  //  if (PtR->poly->sgn_at(b))
      PtR->cut(b);

  return 0;
}

//  int K_POINT1D :: halve() const
//    cuts *this at the point half that halves the interval for *this, i.e.,
//    refines the interval for *this
//      by setting get_low() or get_high() to be half.

int K_POINT1D :: halve() const
{
  assert(type);

  if (type == 1)
  {
    bigrational half;

    half = (PtR->low + PtR->high) / 2;

    if (!PtR->poly->sgn_at(half))
    {
      type = 2;

      if (!--PtR->ref_count)
        delete PtR;

      PtR = 0;
      PtB = half;
    }
    else  //  if (PtR->poly->sgn_at(half))
      PtR->cut(half);
  }

  return 0;
}

//  int K_POINT1D :: reduce(const unsigned long num_bits) const
//    reduces *this s.t.
//      [get_low(), get_high()] will contain all the numbers
//        approximated by float_est to num_bits precision.
//    returns 1 if some reduction occurs and
//            0 otherwise.

int K_POINT1D :: reduce(const unsigned long num_bits) const
{
  assert(type);

  int reduced;

  if (type == 1)
  {
    assert(PtR->ok_float);

    bigrational l, h;

    //  Compute l and h s.t. [l, h] contains all the numbers
    //    approximated by float_est to num_bits precision.

    double2bigrational_interval(PtR->float_est, num_bits, l, h);

    //  If [low, high] \sebseteq [l, h], or
    //     (h >=) l > high (>= low) or (l <=) h < low (<= high) then
    //    nothing to do.
    //  Otherwise, see if *this is reducible.

    if (l <= PtR->low && h >= PtR->high)
      reduced = 1;
    else if (l <= PtR->high && h >= PtR->low)
    //  if (PtR->low < l <= PtR->high && PtR->low <= h < PtR->high)
    {
      int sgn_l, sgn_h;

      reduced = 1;

      if (l > PtR->low)
        if (sgn_l = PtR->poly->sgn_at(l))
          if (PtR->sig_low < 0 && sgn_l < 0
              ||
              PtR->sig_low > 0 && sgn_l > 0
              ||
              PtR->sig_low == 0 &&
              PtR->poly->num_Sturm_seq_perm(l) == PtR->num_perm_low)
          //  Keep reduced == 1.
            PtR->low = l;
          else
            reduced = 0;
        else  //  if l is a root of PtR->poly
        //  Keep reduced == 1.
        {
          type = 2;

          if (!--PtR->ref_count)
            delete PtR;

          PtR = 0;
          PtB = l;
        }

      //  Do the following if PtR != 0 (and even if reduced == 0)
      //  since better high or low might be obtained.

      if (PtR && h < PtR->high)
        if (sgn_h = PtR->poly->sgn_at(h))
          if (PtR->sig_low < 0 && sgn_h > 0
              ||
              PtR->sig_low > 0 && sgn_h < 0
              ||
              PtR->sig_low == 0 &&
              PtR->poly->num_Sturm_seq_perm(h) == PtR->num_perm_high)
          //  Keep reduced as it is.
            PtR->high = h;
          else
            reduced = 0;
        else  //  if h is a root of PtR->poly
        {
          type = 2;

          if (!--PtR->ref_count)
            delete PtR;

          PtR = 0;
          PtB = h;

          reduced = 1;
        }
    }
    else  //  if (h >= l) > high (>= low) or (l <=) h < low (<= high)
      reduced = 0;
  }
  else  //  if (type == 2)
    reduced = 1;

  return reduced;
}

//  int K_POINT1D :: contract(const bigrational& tol) const
//    contracts *this s.t. [get_low(), get_high()] will be no larger than tol.

int K_POINT1D :: contract(const bigrational& tol) const
{
  assert(tol > 0);
  assert(type);

  if (type == 1 && PtR->high - PtR->low > tol)
  {
    if (PtR->ok_float)
    {
      double  d_tol;
      long    num_bits;

      d_tol    = tol.as_double();
      num_bits = 0;

      while (d_tol < 1.0 && num_bits < MAX_NUM_GOOD_FP_BITS)
      {
        d_tol *= 2.0;
        num_bits++;
      }

      while (num_bits > 0 && !reduce(num_bits))
        num_bits -= 3;
    }

    while (type == 1 && PtR->high - PtR->low > tol)
      halve();
  }

  return 0;
}

//  int K_POINT1D :: shrink(const bigrational& fac) const
//    shrinks *this s.t. [get_low(), get_high()] will become smaller by fac.
//    e.g. fac == 1/10 => [get_low(), get_high()] will become 1/10.

int K_POINT1D :: shrink(const bigrational& fac) const
{
  assert(type);

  if (type == 1)
    contract(fac * (PtR->high - PtR->low));

  return 0;
}

//  K_POINT1D K_POINT1D :: add(const K_POINT1D& x) const
//    returns *this + x

K_POINT1D K_POINT1D :: add(const K_POINT1D& x) const
{
  assert(type);
  assert(x.type);

  K_POINT1D y;

  if (type == 1 && x.type == 1)
  {
    y.type = 1;

    long d[2];
    long p[2];

    K_RATPOLY T = poly->add_var(1);
    //  T = T(V0, V1), V1 is dummy

    K_RATPOLY X = x.poly->add_var(0).add_var(2);
    //  X(V0, V1, V2), V0 and V2 are dummy

    //  Form V2 - V0.

    d[0] = 1;
    d[1] = 1;

    K_RATPOLY V2_minus_V0(2, d);

    p[0]                     = 1;
    p[1]                     = 0;
    V2_minus_V0.get_coeff(p) = - 1;
    p[0]                     = 0;
    p[1]                     = 1;
    V2_minus_V0.get_coeff(p) = 1;

    K_RATPOLY X_V1_minus_V0 = X.subst_expr(1, V2_minus_V0);
    //  X_V1_minus_V0(V0, V1) = X(V0, V2 - V0, V2) | V2 = V1 = *x.poly(V1 - V0)

    K_RATPOLY Y = T.GoodSylvester(X_V1_minus_V0, 0);
    //  Y(V0) = Res_V0 (*poly(V0), *x.poly(V1 - V0))

    y.PtR = new ROOT1(Y, PtR->low + x.PtR->low, PtR->high + x.PtR->high);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else if (type == 1 && x.type == 2)
  {
    y.type = 1;

    K_RATPOLY T = poly->add_var(0);
    //  T(V0, V1), V0 is dummy

    K_RATPOLY X = K_RATPOLY(1, 0, x.PtB);
    //  X(V0) = V0 - x.PtB

    K_RATPOLY Y = T.subst_expr(1, X);
    //  Y(V0) = *poly(V0 - x.PtB)

    y.PtR = new ROOT1(Y, PtR->low + x.PtB, PtR->high + x.PtB);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else if (type == 2 && x.type == 1)
  {
    y.type = 1;

    K_RATPOLY T = K_RATPOLY(1, 0, PtB);
    //  T(V0) = V0 - x.PtB

    K_RATPOLY X = x.poly->add_var(0);
    //  X(V0, V1), V0 is dummy

    K_RATPOLY Y = X.subst_expr(1, T);
    //  Y(V0) = *x.poly(V0 - PtB)

    y.PtR = new ROOT1(Y, PtB + x.PtR->low, PtB + x.PtR->high);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else  //  if (type == 2 && x.type == 2)
  {
    y.type = 2;

    y.PtR = 0;
    y.PtB = PtB + x.PtB;

    y.poly = new K_RATPOLY(1, 0, y.PtB);
    y.poly->ref_count++;
  }

  return y;
}

K_POINT1D operator +(const K_POINT1D& x, const K_POINT1D& y)
{
  return x.add(y);
}

//  K_POINT1D K_POINT1D :: sub(const K_POINT1D& x) const
//    return *this - x

K_POINT1D K_POINT1D :: sub(const K_POINT1D& x) const
{
  assert(type);
  assert(x.type);

  K_POINT1D y;

  if (type == 1 && x.type == 1)
  {
    y.type = 1;

    long d[2];
    long p[2];

    K_RATPOLY T = poly->add_var(1);
    //  T(V0, V1), V1 is dummy

    K_RATPOLY X = x.poly->add_var(0).add_var(2);
    //  X(V0, V1, V2), V0 and V2 are dummy

    //  Form V2 + V0.

    d[0] = 1;
    d[1] = 1;

    K_RATPOLY V2_plus_V0(2, d);

    p[0]                    = 1;
    p[1]                    = 0;
    V2_plus_V0.get_coeff(p) = 1;
    p[0]                    = 0;
    p[1]                    = 1;
    V2_plus_V0.get_coeff(p) = 1;

    K_RATPOLY X_V1_plus_V0 = X.subst_expr(1, V2_plus_V0);
    //  X_V1_plus_V0(V0, V1) = X(V0, V2 + V0, V2) | V2 = V1 = *x.poly(V1 + V0)

    K_RATPOLY Y = T.GoodSylvester(X_V1_plus_V0, 0);
    //  Y(V0) = Res_V0 (*poly(V0), *x.poly(V1 + V0))

    y.PtR = new ROOT1(Y, PtR->low - x.PtR->high, PtR->high - x.PtR->low);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else if (type == 1 && x.type == 2)
  {
    y.type = 1;

    K_RATPOLY T = poly->add_var(0);
    //  T(V0, V1), V0 is dummy

    K_RATPOLY X = K_RATPOLY(1, 0, - x.PtB);
    //  X(V0) = V0 + x.PtB

    K_RATPOLY Y = T.subst_expr(1, X);
    //  Y(V0) = *poly(V0 + x.PtB)

    y.PtR = new ROOT1(Y, PtR->low - x.PtB, PtR->high - x.PtB);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else if (type == 2 && x.type == 1)
  {
    y.type = 1;

    K_RATPOLY T = K_RATPOLY(1, 0, - PtB);
    //  T(V0) = V0 + x.PtB

    K_RATPOLY X = x.poly->add_var(0);
    //  X(V0, V1), V0 is dummy

    K_RATPOLY Y = X.subst_expr(1, T);
    //  Y(V0) = *x.poly(V0 + PtB)

    y.PtR = new ROOT1(Y, x.PtR->low - PtB, x.PtR->high - PtB);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else  //  if (type == 2 && x.type == 2)
  {
    y.type = 2;

    y.PtR  = 0;
    y.PtB  = PtB - x.PtB;

    y.poly = new K_RATPOLY(1, 0, y.PtB);
    y.poly->ref_count++;
  }

  return y;
}

K_POINT1D operator -(const K_POINT1D& x, const K_POINT1D& y)
{
  return x.sub(y);
}

//  K_POINT1D K_POINT1D :: mul(const K_POINT1D& x) const
//    return *this * x

K_POINT1D K_POINT1D :: mul(const K_POINT1D& x) const
{
  assert(type);
  assert(x.type);

  //  Cut *this and x at 0 s.t. their low and high have the same sign.

  cut(0);
  x.cut(0);

  K_POINT1D y;

  if (type == 1 && x.type == 1)
  {
    y.type = 1;

    long          d[2];
    long          p[2];
    bigrational   ly, hy;

    K_RATPOLY T = poly->add_var(1);
    //  T(V0, V1), V1 is dummy

    K_RATPOLY X = x.poly->add_var(0).add_var(2);
    //  X(V0, V1, V2), V0 and V2 are dummy

    //  Form V2 / V0.

    d[0] = 0;
    d[1] = 1;

    K_RATPOLY V2(2, d);

    p[0]            = 0;
    p[1]            = 1;
    V2.get_coeff(p) = 1;

    d[0] = 1;
    d[1] = 0;

    K_RATPOLY V0(2, d);

    p[0]            = 1;
    p[1]            = 0;
    V0.get_coeff(p) = 1;

    K_RATPOLY X_V1_over_V0 = X.subst_expr(1, V2, V0);
    //  X_V1_over_V0(V0, V1) = X(V0, V2 / V0, V2) | V2 = V1
    //                       = V0^(deg *x.poly) *x.poly(V1 / V0)

    K_RATPOLY Y = T.GoodSylvester(X_V1_over_V0, 0);
    //  Y(V0) = Res_V0 (*poly(V0), V0^(deg *x.poly) *x.poly(V1 / V0))

    assert(PtR->low * PtR->high >= 0);
    assert(x.PtR->low * x.PtR->high >= 0);

    if (PtR->high > 0)
      if (x.PtR->high > 0)  //  if *this > 0 && x > 0
      {
        ly = PtR->low  * x.PtR->low;
        hy = PtR->high * x.PtR->high;
      }
      else  //  if *this > 0 && x <= 0
      {
        ly = PtR->high * x.PtR->low;
        hy = PtR->low  * x.PtR->high;
      }
    else  //  if *this <= 0 && x > 0
      if (x.PtR->high > 0)
      {
        ly = PtR->low  * x.PtR->high;
        hy = PtR->high * x.PtR->low;
      }
      else  //  if *this <= 0 && x <= 0
      {
        ly = PtR->high * x.PtR->high;
        hy = PtR->low  * x.PtR->low;
      }

    y.PtR = new ROOT1(Y, ly, hy);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else if (type == 1 && x.type == 2)
    if (sgn(x.PtB))
    {
      y.type = 1;

      long        dx[1];
      bigrational ly, hy;

      K_RATPOLY T = poly->add_var(0);
      //  T(V0, V1), V0 is dummy

      //  Form X(V0) = V0 / x.PtB

      dx[0] = 1;

      K_RATPOLY X = K_RATPOLY(1, dx);

      X.coeffs[0] = 1 / x.PtB;

      K_RATPOLY Y = T.subst_expr(1, X);
      //  Y(V0) = *poly(V0 / x.PtB)

      assert(PtR->low * PtR->high >= 0);

      if (x.PtB > 0)
      {
        ly = PtR->low  * x.PtB;
        hy = PtR->high * x.PtB;
      }
      else  //  if (x.PtB < 0)
      {
        ly = PtR->high * x.PtB;
        hy = PtR->low  * x.PtB;
      }

      y.PtR = new ROOT1(Y, ly, hy);
      y.PtR->ref_count++;

      if (y.poly = y.PtR->poly)
        y.poly->ref_count++;
    }
    else  //  if (x.PtB == 0)
    {
      y.type = 2;

      y.PtR = 0;
      y.PtB = 0;

      y.poly = new K_RATPOLY(1, 0, 0);
      y.poly->ref_count++;
    }
  else if (type == 2 && x.type == 1)
    if (sgn(PtB))
    {
      y.type = 1;

      long        dt[1];
      bigrational ly, hy;

      //  Form T(V0) = V0 / PtB

      dt[0] = 1;

      K_RATPOLY T = K_RATPOLY(1, dt);

      T.coeffs[0] = 1 / PtB;

      K_RATPOLY X = x.poly->add_var(0);
      //  X(V0, V1), V0 is dummy

      K_RATPOLY Y = X.subst_expr(1, T);
      //  Y(V0) = *x.poly(V0 / PtB)

      assert(x.PtR->low * x.PtR->high >= 0);

      if (PtB > 0)
      {
        ly = PtB * x.PtR->low;
        hy = PtB * x.PtR->high;
      }
      else  //  if PtB < 0
      {
        ly = PtB * x.PtR->high;
        hy = PtB * x.PtR->low;
      }

      y.PtR = new ROOT1(Y, ly, hy);
      y.PtR->ref_count++;

      if (y.poly = y.PtR->poly)
        y.poly->ref_count++;
    }
    else  //  if (PtB == 0)
    {
      y.type = 2;

      y.PtR = 0;
      y.PtB = 0;

      y.poly = new K_RATPOLY(1, 0, 0);
      y.poly->ref_count++;
    }
  else  //  if (type == 2 && x.type == 2)
  {
    y.type = 2;

    y.PtB = PtB * x.PtB;
    y.PtR = 0;

    y.poly = new K_RATPOLY(1, 0, y.PtB);
    y.poly->ref_count++;
  }

  return y;
}

K_POINT1D operator *(const K_POINT1D& x, const K_POINT1D& y)
{
  return x.mul(y);
}

//  K_POINT1D K_POINT1D :: div(const K_POINT1D& x) const
//    return *this / x

K_POINT1D K_POINT1D :: div(const K_POINT1D& x) const
{
  assert(type);
  assert(x.type);

  //  Cut *this and x at 0 s.t. their low & high have the same sign.

  cut(0);
  x.cut(0);

  assert(x.type == 1 || x.type == 2 && sgn(x.PtB));

  while (x.type == 1 && (!sgn(x.PtR->low) || !sgn(x.PtR->high)))
    x.shrink(shrink_step);

  K_POINT1D y;

  if (type == 1 && x.type == 1)
  {
    y.type = 1;

    long          d[2];
    long          p[2];
    bigrational   ly, hy;

    K_RATPOLY T = poly->add_var(1);
    //  T(V0, V1), V1 is dummy

    K_RATPOLY X = x.poly->add_var(0).add_var(2);
    //  X(V0, V1, V2), V0 and V2 are dummy

    //  Form V0 / V2.

    d[0] = 1;
    d[1] = 0;

    K_RATPOLY V0(2, d);

    p[0]            = 1;
    p[1]            = 0;
    V0.get_coeff(p) = 1;

    d[0] = 0;
    d[1] = 1;

    K_RATPOLY V2(2, d);

    p[0]            = 0;
    p[1]            = 1;
    V2.get_coeff(p) = 1;

    K_RATPOLY X_V0_over_V1 = X.subst_expr(1, V0, V2);
    //  X_V0_over_V1(V0, V1) = X(V0, V0 / V2, V2) | V2 = V1
    //                       = V1^(deg *x.poly) *x.poly(V0 / V1)

    K_RATPOLY Y = T.GoodSylvester(X_V0_over_V1, 0);
    //  Y(V0) = Res_V0 (*poly(V0), V1^(deg *x.poly) *x.poly(V0 / V1))

    assert(PtR->low * PtR->high >= 0);
    assert(x.PtR->low * x.PtR->high >= 0);

    if (PtR->high > 0)
      if (x.PtR->high > 0)  //  if *this > 0 && x > 0
      {
        ly = PtR->low  / x.PtR->high;
        hy = PtR->high / x.PtR->low;
      }
      else  //  if *this > 0 && x < 0
      {
        ly = PtR->high / x.PtR->high;
        hy = PtR->low  / x.PtR->low;
      }
    else  //  if *this <= 0 && x > 0
      if (x.PtR->high > 0)
      {
        ly = PtR->low  / x.PtR->low;
        hy = PtR->high / x.PtR->high;
      }
      else  //  if *this <= 0 && x < 0
      {
        ly = PtR->high / x.PtR->low;
        hy = PtR->low  / x.PtR->high;
      }

    y.PtR = new ROOT1(Y, ly, hy);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else if (type == 1 && x.type == 2)
  {
    y.type = 1;

    long        dx[1];
    bigrational ly, hy;

    K_RATPOLY T = poly->add_var(0);
    //  T(V0, V1), V0 is dummy

    //  Form X(V0) = V0 * x.PtB

    dx[0] = 1;

    K_RATPOLY X = K_RATPOLY(1, dx);

    X.coeffs[0] = x.PtB;

    K_RATPOLY Y = T.subst_expr(1, X);
    //  Y(V0) = *poly(V0 * x.PtB)

    assert(PtR->low * PtR->high >= 0);

    if (x.PtB > 0)
    {
      ly = PtR->low  / x.PtB;
      hy = PtR->high / x.PtB;
    }
    else  //  if (x.PtB < 0)
    {
      ly = PtR->high / x.PtB;
      hy = PtR->low  / x.PtB;
    }

    y.PtR = new ROOT1(Y, ly, hy);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else if (type == 2 && x.type == 1)
  {
    y.type = 1;

    long        dt[1];
    bigrational ly, hy;

    //  Form T(V0) = V0 * PtB

    dt[0] = 1;

    K_RATPOLY T = K_RATPOLY(1, dt);

    T.coeffs[0] = PtB;

    K_RATPOLY X = x.poly->add_var(0);
    //  X(V0, V1), V0 is dummy

    K_RATPOLY Y = X.subst_expr(1, T);
    //  Y(V0) = *x.poly(V0 * PtB)

    assert(x.PtR->low * x.PtR->high >= 0);

    if (PtB > 0)
    {
      ly = x.PtR->low  / PtB;
      hy = x.PtR->high / PtB;
    }
    else  //  if (PtB < 0)
    {
      ly = x.PtR->high / PtB;
      hy = x.PtR->low  / PtB;
    }

    y.PtR = new ROOT1(Y, ly, hy);
    y.PtR->ref_count++;

    if (y.poly = y.PtR->poly)
      y.poly->ref_count++;
  }
  else  //  if (type == 2 && x.type == 2)
  {
    y.type = 2;

    y.PtR = 0;
    y.PtB = PtB / x.PtB;

    y.poly = new K_RATPOLY(1, 0, y.PtB);
    y.poly->ref_count++;
  }

  return y;
}

K_POINT1D operator /(const K_POINT1D& x, const K_POINT1D& y)
{
  return x.div(y);
}

//  int K_POINT1D :: compare(const K_POINT1D& x) const
//    return   1 if *this > x,
//             0 if *this == x, and
//           - 1 if *this < x.

int K_POINT1D :: compare(const K_POINT1D& x) const
{
  assert(type);
  assert(x.type);

  int c;

  if (type == 2)
    if (x.type == 2)
      if (PtB > x.PtB)
        c = 1;
      else if (PtB < x.PtB)
        c = - 1;
      else  //  if (PtB == x.PtB)
        c = 0;
    else  //  if (x.type == 1)
      if (PtB >= x.PtR->high)
        c = 1;
      else if (PtB <= x.PtR->low)
        c = - 1;
      else  //  if (x.PtR->low < PtB < x.PtR->high)
      {
        x.cut(PtB);
        c = compare(x);
      }
  else  //  if (type == 1)
    if (x.type == 2)
      if (PtR->low >= x.PtB)
        c = 1;
      else if (PtR->high <= x.PtB)
        c = - 1;
      else  //  if (PtR->low < x.PtB < PtR->high)
      {
        cut(x.PtB);
        c = compare(x);
      }
    else  //  if (x.type == 1)
      if (PtR->low > x.PtR->low && PtR->low < x.PtR->high)
      {
        x.cut(PtR->low);
        c = compare(x);
      }
      else if (PtR->high > x.PtR->low && PtR->high < x.PtR->high)
      {
        x.cut(PtR->high);
        c = compare(x);
      }
      else if (PtR->low < x.PtR->low && PtR->high > x.PtR->low)
      {
        cut(x.PtR->low);
        c = compare(x);
      }
      else if (PtR->low < x.PtR->high && PtR->high > x.PtR->high)
      {
        cut(x.PtR->high);
        c = compare(x);
      }
      else if (PtR->low >= x.PtR->high)
        c = 1;
      else if (PtR->high <= x.PtR->low)
        c = - 1;
      else  //  if (PtR->low == x.PtR->low && PtR->high == x.PtR->high)
        c = 0;

  return c;
}

//  int K_POINT1D :: compare(const bigrational& b) const
//    return   1 if *this >  b,
//             0 if *this == b, and
//           - 1 if *this <  b.

int K_POINT1D :: compare(const bigrational& b) const
{
  assert(type);

  int c;

  if (type == 2)
    if (PtB > b)
      c = 1;
    else if (PtB < b)
      c = - 1;
    else  //  if (PtB == b)
      c = 0;
  else  //  if (type == 1)
    if (PtR->low >= b)
      c = 1;
    else if (PtR->high <= b)
      c = - 1;
    else  //  if (PtR->low < b < PtR->high)
    {
      cut(b);
      c = compare(b);
    }

  return c;
}

//  int sort(K_POINT1D** const X, const unsigned long n)
//    insertion-sort the array X of length n.
//    return 1 if all the elements are distinct and
//           0 otherwise.

int sort(K_POINT1D** const X, const unsigned long n)
{
  long       i, j;
  K_POINT1D* x;
  int        c, distinct;

  distinct = 1;

  for (i = 1; i < n; i++)
  {
    x = X[i];

    for (j = i - 1; j >= 0 && (c = x->compare(*X[j])) < 0; j--)
      X[j + 1] = X[j];

    X[j + 1] = x;

    if (!c)
      distinct = 0;
  }

  return distinct;
}

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

unsigned long get_pts(const bigrational& l, const bigrational& h,
                      const K_RATPOLY& P,
                      K_POINT1D**& pts,
                      const bigrational& tol, const int count_endpts)
{
  assert(P.num_vars == 1);

#ifdef _EXPERIMENT
  num_kpoint1d_get_pts++;
#endif

  unsigned long num_pts;

  if (!P.deg[0])
  //  0.1.  when P is a constant polynomial
  {
    pts     = 0;
    num_pts = 0;
  }
  else if (l > h)
  //  0.2.  when the interval [l, h] is not well-defined
  {
    pts     = 0;
    num_pts = 0;
  }
  else
  {
    //  1.  See if 0 is a root of P.
    //      Obtain P1 by removing some power of s from P s.t.
    //        P1 has a non-zero constant term.

    long          i;
    unsigned long j;
    long          d[1];
    K_POINT1D**   pts1;
    unsigned long num_pts0, num_pts1;

    for (i = P.num_coeffs - 1; i >= 0 && !sgn(P.coeffs[i]); i--)
      ;

    d[0] = i;

    K_RATPOLY P1(1, d);

    for (j = 0; j < P1.num_coeffs; j++)
      P1.coeffs[j] = P.coeffs[j];

    if (i < P.num_coeffs - 1
        &&
        (sgn(l) < 0 && sgn(h) > 0
         ||
         count_endpts && sgn(l) <= 0 && sgn(h) >= 0))
    //  if 0 is a root of P in/on the interval [l, h]
      num_pts0 = 1;
    else
    //  if 0 is not a root of P or 0 is not in/on the interval [l. h]
      num_pts0 = 0;

    //  2.  Find roots of P1.

    if (!P1.deg[0])
    //  2.0.  when P1 is a constant polynomial
    {
      pts1     = 0;
      num_pts1 = 0;
    }
    else if (P1.deg[0] == 1)
    //  2.1.  when P1 is linear
    {
      bigrational b;

      b = - P1.coeffs[1] / P1.coeffs[0];
      //  P1.deg[0] == 1 => P1.coeffs[0] != 0

      if (b > l && b < h || count_endpts && b >= l && b <= h)
      {
        pts1    = new K_POINT1D* [num_pts1 = 1];
        pts1[0] = new K_POINT1D(b, P1);
        pts1[0]->ref_count++;
      }
      else  //  if b is not in/on the interval [l, s]
      {
        pts1     = 0;
        num_pts1 = 0;
      }
    }
    else
    //  2.2.  Otherwise, perform 1-D root finding.
    {
      unsigned long num_each, num_all;
      ROOT1*        each;
      bigrational   e_l;

      //  2.2.1  Find all roots of P1 in/on the half open interval [l, h) and
      //         isolate them.

      ROOT1 all(P1, l, h);

      num_all  = all.num_roots;
      num_each = all.isolate_roots(each, tol);

      assert(num_each == num_all);

      pts1     = new K_POINT1D* [num_all + 1];
      num_pts1 = 0;

      for (i = 0; i < num_all; i++)
      {
        e_l = each[i].low;

        if (!P1.sgn_at(e_l))
        //  2.2.1.1  A root e_l = each[i].low of P1 will belong to pts1
        //           iff (1) it is not zero (it was taken care of before) and
        //               (2) it is not l unless count_endpts == 1.
        {
          if (sgn(e_l) && (count_endpts || e_l != l))
          {
            pts1[num_pts1] = new K_POINT1D(e_l);
            pts1[num_pts1++]->ref_count++;
          }
        }
        else  //  if (P1.sgn_at(e_l))
        {
          pts1[num_pts1] = new K_POINT1D(each[i]);
          pts1[num_pts1++]->ref_count++;
        }
      }

      if (num_each > 0)
        delete [] each;

      //  2.2.2  See if h is a root of P1.

      if (count_endpts && !P1.sgn_at(h))
      {
        pts1[num_pts1] = new K_POINT1D(h);
        pts1[num_pts1++]->ref_count++;
      }
    }

    //  3.  Gather all together.

    pts = new K_POINT1D* [num_pts = num_pts0 + num_pts1];
    i   = 0;

    if (num_pts0 == 1)
    {
      pts[0] = new K_POINT1D(0);
      pts[0]->ref_count++;
      i      = 1;
    }

    for (j = 0; j < num_pts1; j++)
    {
      pts[i] = pts1[j];
      pts[i++]->ref_count++;
    }

    assert(i == num_pts);

    for (i = 0; i < num_pts1; i++)
      if (!--pts1[i]->ref_count)
        delete pts1[i];

    delete [] pts1;  //  pts1 != 0
  }

  return num_pts;
}

unsigned long match_intervals(const bigrational* const L,
                              const bigrational* const H,
                              const unsigned long num_intervals,
                              K_POINT1D** const pts,
                              long*& matches,
                              const unsigned long num_pts)
{
//  assert(sort(pts, num_pts));

  unsigned long i, j;
  int           c;

  matches = new long [num_pts];

  i = j = 0;

  while (i < num_intervals && j < num_pts)
  {
    while (i < num_intervals && pts[j]->get_low() >= H[i])
      i++;

    if (i < num_intervals)
    {
      c = 1;

      while (c && i < num_intervals)
      {
        pts[j]->cut(L[i]);
        pts[j]->cut(H[i]);

        if (pts[j]->get_low() >= H[i])
          i++;
        else if (pts[j]->get_high() <= L[i])
        {
          matches[j] = - 1;
          c          = 0;
        }
        else  //  if (L[i] < *pts[j] < H[i])
        {
          matches[j] = i;
          c          = 0;
        }
      }
    }

    if (i < num_intervals)
      j++;
  }

  if (i == num_intervals)
    while (j < num_pts)
      matches[j++] = - 1;

  return num_pts;
}

//  int K_POINT1D :: overlap(const K_POINT1D& x) const
//    return 1 if *this and x overlap, and
//           0 otherwise.

int K_POINT1D :: overlap(const K_POINT1D& x) const
{
  assert(type);
  assert(x.type);

  int o;

  if (type == 2)
    if (x.type == 2)
      if (PtB == x.PtB)
        o = 1;
      else  //  if (PtB != x.PtB)
        o = 0;
    else  //  if (x.type == 1)
      if (PtB > x.PtR->low && PtB < x.PtR->high)
        o = 1;
      else  //  if (PtB <= x.PtR->low || PtB >= x.PtR->high)
        o = 0;
  else  //  if (type == 1)
    if (x.type == 2)
      if (PtR->low < x.PtB && PtR->high > x.PtB)
        o = 1;
      else  //  if (PtR->low >= x.PtB || PtR->high <= x.PtB)
        o = 0;
    else  //  if (x.type == 1)
      if (PtR->low >= x.PtR->low && PtR->low < x.PtR->high
          ||
          PtR->high > x.PtR->low && PtR->high <= x.PtR->high
          ||
          PtR->low <= x.PtR->low && PtR->high > x.PtR->low
          ||
          PtR->low < x.PtR->high && PtR->high >= x.PtR->high)
        o = 1;
      else  //  if (PtR->low > x.PtR->high || PtR->high < x.PtR->low)
        o = 0;

  return o;
}

//  int K_POINT1D :: separate(const K_POINT1D& x) const
//    separate *this and x s.t. they will not overlap.
//    POSSIBLY DOES NOT TERMINATE!

int K_POINT1D :: separate(const K_POINT1D& x) const
{
  assert(type);
  assert(x.type);

  if (type == 2 && x.type == 1)
    x.cut(PtB);
  else if (type == 1 && x.type == 2)
    cut(x.PtB);
  else if (type == 1 && x.type == 1)
  {
    cut(x.PtR->low);
    cut(x.PtR->high);
    x.cut(PtR->low);
    x.cut(PtR->high);
  }

  while (overlap(x))
  {
    shrink(shrink_step);
    x.shrink(shrink_step);
    separate(x);
  }

  return 0;
}

//  unsigned long get_all_pts(const K_RATPOLY& P,
//                            K_POINT1D**& pts,
//                            const bigrational& tol)
//    computes all the intersections "pts" of a univariate polynomials P.
//    returns the number of roots.
//    if tol > 0 then an interval for each root is no larger than tol.
//    if tol = 0 then there is no limit on the size of intervals for roots.

unsigned long get_all_pts(const K_RATPOLY& P,
                          K_POINT1D**& pts,
                          const bigrational tol)
{
  assert(P.get_num_vars() == 1);

  unsigned long num_int_pts;
  bigrational   l, h;

  h = P.get_Mignotte_bd() + epsilon;
  l = - h - epsilon;

  return get_pts(l, h, P, pts, tol, 0);
}

