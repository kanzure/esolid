//  file:    root1.cc
//  update:  10/08/02

#include <config.h>
#ifdef _EXPERIMENT
#include <counter.h>
#endif

#include <root1.h>

#include <fpconversion.h>

//  ROOT1 :: ROOT1()

ROOT1 :: ROOT1()
  : poly(0),
    low(0), high(0),
    num_perm_low(- 1), num_perm_high(- 1),
    num_roots(0),
    sig_low(0),
    float_est(0.0), ok_float(0),
    ref_count(0)
{}

//  ROOT1 :: ROOT1(const K_RATPOLY& P,
//                 const bigrational& l, const bigrational& h)
//    constructs an instance for the set of the roots of P in [l, h].
//
//ROOT1 :: ROOT1(const K_RATPOLY& P, const bigrational& l, const bigrational& h)
//{
//  assert(P.num_vars == 1);
//  assert(P.deg[0] >= 1);
//  assert(l <= h);
//  
//  poly = new K_RATPOLY(P);
//  poly->ref_count++;
//  
//  low  = l;
//  high = h;
//  
//  poly->set_Sturm_seq();
//  
//  num_perm_low  = poly->num_Sturm_seq_perm(low);
//  num_perm_high = poly->num_Sturm_seq_perm(high);
//  
//  num_roots = num_perm_high - num_perm_low;
//  
//  if (num_roots == 1)
//  {
//    int sgn_low, sgn_high;
//    
//    sgn_low  = poly->sgn_at(low);
//    sgn_high = poly->sgn_at(high);
//    
//    if (sgn_low < 0 && sgn_high > 0)
//      sig_low = - 1;
//    else if (sgn_low > 0 && sgn_high < 0)
//      sig_low = 1;
//    else  //  if (sgn_low * sgn_high >= 0)
//      sig_low = 0;
//    
//    unsigned long num_fp_roots;
//    double*       fp_roots;
//    
//    num_fp_roots =
//      poly->as_FLOATPOLY().gen_fp_roots(low.as_double(), high.as_double(),
//                                        fp_roots);
//    
//    if (num_fp_roots == 1)
//    {
//      float_est = fp_roots[0];
//      ok_float  = 1;
//      delete [] fp_roots;
//    }
//    else  //  if (num_fp_roots != 1)
//    {
//      float_est = 0.0;
//      ok_float  = 0;
//    }
//  }
//  else  //  if (num_roots != 1)
//  {
//    sig_low   = 0;
//    
//    float_est = 0.0;
//    ok_float  = 0;
//  }
//  
//  ref_count = 0;
//}

//  ROOT1 :: ROOT1(const K_RATPOLY& P,
//                 const bigrational& l, const bigrational& h,
//                 const long n_l, const long n_h)
//    constructs an instance for the set of the roots of P in [l, h].
//    the num. of Sturm seq. permanencies at l and h are n_l and n_h
//      unless n_l and n_h are - 1, resp.

ROOT1 :: ROOT1(const K_RATPOLY& P,
               const bigrational& l, const bigrational& h,
               const long n_l, const long n_h)
{
  assert(P.num_vars == 1);
  assert(P.deg[0] >= 1);
  assert(l <= h);
  
  poly = new K_RATPOLY(P);
  poly->ref_count++;
  
  low  = l;
  high = h;
  
  poly->set_Sturm_seq();
  
  if (n_l < 0)
    num_perm_low = poly->num_Sturm_seq_perm(low);
  else  //  if (n_l >= 0)
    num_perm_low = n_l;
  
  if (n_h < 0)
    num_perm_high = poly->num_Sturm_seq_perm(high);
  else  //  if (n_h >= 0)
    num_perm_high = n_h;
  
  num_roots = num_perm_high - num_perm_low;
  
  if (num_roots == 1)
  {
    int sgn_low, sgn_high;
    
    sgn_low  = poly->sgn_at(low);
    sgn_high = poly->sgn_at(high);
    
    if (sgn_low < 0 && sgn_high > 0)
      sig_low = - 1;
    else if (sgn_low > 0 && sgn_high < 0)
      sig_low = 1;
    else  //  if (sgn_low * sgn_high >= 0)
      sig_low = 0;
    
    unsigned long num_fp_roots;
    double*       fp_roots;
    
    num_fp_roots =
      poly->as_FLOATPOLY().gen_fp_roots(low.as_double(), high.as_double(),
                                        fp_roots);
    
    if (num_fp_roots == 1)
    {
      float_est = fp_roots[0];
      ok_float  = 1;
      delete [] fp_roots;
    }
    else  //  if (num_fp_roots != 1)
    {
      float_est = 0.0;
      ok_float  = 0;
    }
  }
  else  //  if (num_roots != 1)
  {
    sig_low   = 0;
    
    float_est = 0.0;
    ok_float  = 0;
  }
  
  ref_count = 0;
}

//  ROOT1 :: ROOT1(const K_RATPOLY& P,
//                 const bigrational& l, const bigrational& h,
//                 const double d)
//    constructs an instance for the set of the roots of P in [l, h].
//    the floating-point estimate of the root is d if num_roots == 1.

ROOT1 :: ROOT1(const K_RATPOLY& P,
               const bigrational& l, const bigrational& h,
               const double d)
{
  assert(P.num_vars == 1);
  assert(P.deg[0] >= 1);
  assert(l <= h);
  
  poly = new K_RATPOLY(P);
  poly->ref_count++;
  
  low  = l;
  high = h;
  
  poly->set_Sturm_seq();
  
  num_perm_low  = poly->num_Sturm_seq_perm(low);
  num_perm_high = poly->num_Sturm_seq_perm(high);
  
  num_roots = num_perm_high - num_perm_low;
  
  if (num_roots == 1)
  {
    int sgn_l, sgn_h;
    
    sgn_l = poly->sgn_at(low);
    sgn_h = poly->sgn_at(high);
    
    if (sgn_l < 0 && sgn_h > 0)
      sig_low = - 1;
    else if (sgn_l > 0 && sgn_h < 0)
      sig_low = 1;
    else  //  if (sgn_low * sgn_high >= 0)
      sig_low = 0;
    
    float_est = d;
    ok_float  = 1;
  }
  else  //  if (num_roots != 1)
  {
    sig_low   = 0;
    
    float_est = 0.0;
    ok_float  = 0;
  }
  
  ref_count = 0;
}

//  ROOT1 :: ROOT1(const ROOT1& r)
//    the copy constructor

ROOT1 :: ROOT1(const ROOT1& r)
{
  if (poly = r.poly)
    poly->ref_count++;
  
  low  = r.low;
  high = r.high;
  
  num_perm_low  = r.num_perm_low;
  num_perm_high = r.num_perm_high;
  
  num_roots = r.num_roots;
  
  sig_low = r.sig_low;
  
  float_est = r.float_est;
  ok_float  = r.ok_float;
  
  ref_count = 0;
}

//  ROOT1& ROOT1 :: operator =(const ROOT1& r)
//    the assignment operator

ROOT1& ROOT1 :: operator =(const ROOT1& r)
{
  if (this != &r)
  {
    if (poly && !--poly->ref_count)
      delete poly;
    
    if (poly = r.poly)
      poly->ref_count++;
    
    low  = r.low;
    high = r.high;
    
    num_perm_low  = r.num_perm_low;
    num_perm_high = r.num_perm_high;
    
    num_roots = r.num_roots;
    
    sig_low = r.sig_low;
    
    float_est = r.float_est;
    ok_float  = r.ok_float;
  }
  
  return *this;
}

//  ROOT1 :: ~ROOT1()
//    the destructor

ROOT1 :: ~ROOT1()
{
  if (poly && !--poly->ref_count)
    delete poly;
}

ostream& ROOT1 :: output(ostream& o) const
{
  o << "[" << low << ", " << high << "]" << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const ROOT1& r)
{
  return r.output(o);
}

//  int ROOT1 :: cut(const bigrational& b)
//    cuts [low, high] at b s.t. [low, high] becomes
//      either [low, b] or [b, high] that contains some roots of *poly.

int ROOT1 :: cut(const bigrational& b)
{
  if (b > low && b < high)
  {
    if (sig_low > 0)
      if (poly->sgn_at(b) >= 0)
      //  b is of the same sign as low, i.e., b is below the root,
      //  or b is the root.
      //  Assign b to low.
        low = b;
      else  //  if (poly->sgn_at(b) < 0)
      //  b is of the opposite sign to low, i.e., b is above the root.
      //  Assign b to high.
        high = b;
    else if (sig_low < 0)
      if (poly->sgn_at(b) <= 0)
      //  b is of the same sign as low, i.e., b is below the root,
      //  or b is the root.
      //  Assign b to low.
        low = b;
      else  //  if (poly->sgn_at(b) > 0)
      //  b is of the opposite sign to low, i.e., b is above the root.
      //  Assign b to high.
        high = b;
    else  //  if (sig_low == 0)
    {
      unsigned long num_perm_b;
      
      num_perm_b = poly->num_Sturm_seq_perm(b);
      
      if (num_roots == num_perm_b - num_perm_low)
      //  All the roots are between low and b.
      {
        high          = b;
        num_perm_high = num_perm_b;
      }
      else  //  if (num_roots != num_perm_b - num_perm_low)
      //  Some roots are between b and high.
      {
        low          = b;
        num_perm_low = num_perm_b;
      }
    }
    
    if (ok_float)
      if (float_est > low && float_est < high)  //  comparison in bigrational
        ok_float = 1;
      else  //  if (float_est <= low || float_est >= high)
        ok_float = 0;
  }
  
  return 0;
}

//  int ROOT1 :: halve()
//    cuts [low, high] at half = (low + high) / 2 s.t. [low, high] becomes
//      either [low, half] or [half, high] that contains some roots of *poly.

int ROOT1 :: halve()
{
  bigrational half;
  
  half = (low + high) / 2;
  
  if (sig_low > 0)
    if (poly->sgn_at(half) >= 0)
    //  half is of the same sign as low, i.e., half is below the root,
    //  or half is the root.
    //  Assign half to low.
      low = half;
    else  //  if (poly->sgn_at(half) < 0)
    //  half is of the opposite sign to low, i.e., half is above the root.
    //  Assign half to high.
      high = half;
  else if (sig_low < 0)
    if (poly->sgn_at(half) <= 0)
    //  half is of the same sign as low, i.e., half is below the root,
    //  or half is the root.
    //  Assign half to low.
      low = half;
    else  //  if (poly->sgn_at(half) > 0)
    //  half is of the opposite sign to low, i.e., half is above the root.
    //  Assign half to high.
      high = half;
  else  //  if (sig_low == 0)
  {
    unsigned long num_perm_half;
    
    num_perm_half = poly->num_Sturm_seq_perm(half);
    
    if (num_roots == num_perm_half - num_perm_low)
    //  All the roots are between low and half.
    {
      high          = half;
      num_perm_high = num_perm_half;
    }
    else  //  if (num_roots != num_perm_half - num_perm_low)
    //  Some roots are between half and high.
    {
      low          = half;
      num_perm_low = num_perm_half;
    }
  }
  
  if (ok_float)
    if (float_est > low && float_est < high)  //  comparison in bigrational
      ok_float = 1;
    else  //  if (float_est <= low || float_est >= high)
      ok_float = 0;
  
  return 0;
}

//  int ROOT1 :: reduce(const unsigned long num_bits)
//    reduces [low, high] s.t.
//      [low, high] contains all the numbers
//        approximated by float_est to num_bits precision.
//    returns 1 if *this is reduced and
//           0 otherwise.
//    works only when float_est has been computed.
//    returns 0 if the float_est is not good.

int ROOT1 :: reduce(const unsigned long num_bits)
{
  assert(num_roots == 1);
  assert(ok_float);
  
  bigrational l, h;
  int         reduced;
  
  //  Compute l and h s.t. [l, h] contains all the numbers
  //    approximated by float_est to num_bits precision.
  
  double2bigrational_interval(float_est, num_bits, l, h);
  
  //  If [low, high] \sebseteq [l, h] then nothing to do.
  //  Otherwise, see if *this is reducible.
  
  reduced = 1;
  
  if (l >= low || h <= high)
    if (sig_low)
    {
      if (l > low)
      {
        int sgn_l;
        
        sgn_l = poly->sgn_at(l);
        
        if (sig_low < 0 && sgn_l <= 0 || sig_low > 0 && sgn_l >= 0)
        //  l is of the same sign as low, i.e., l is below the root,
        //  or l is the root.
        //  Assign l to low.
          low = l;
        else
          reduced = 0;
      }
      
      //  Do the following even if reduced has been set to be 0
      //  since better high or low might be obtained.
      
      if (h < high)
      {
        int sgn_h;
        
        sgn_h = poly->sgn_at(h);
        
        if (sig_low < 0 && sgn_h > 0 || sig_low > 0 && sgn_h < 0)
        //  h is of the opposite sign to low, i.e., h is above the root.
        //  Assign h to high.
          high = h;
        else if (!sgn_h)
        //  h is the root.
        //  Assign h to low.
          low = h;
        else
          reduced = 0;
      }
    }
    else  //  sig_low == 0
    {
      unsigned long num_perm_l, num_perm_h;
      
      //  Refine [l, h] to be [low, high] \cap [l, h].
      
      if (l < low)
      {
        l          = low;
        num_perm_l = num_perm_low;
      }
      else  //  if (l >= low)
        num_perm_l = poly->num_Sturm_seq_perm(l);
      
      if (h > high)
      {
        h          = high;
        num_perm_h = num_perm_high;
      }
      else  //  if (h <= high)
        num_perm_h = poly->num_Sturm_seq_perm(h);
      
      //  See if all the roots are between (refined) l and h.
      
      if (num_roots == num_perm_h - num_perm_l)
      {
        low           = l;
        high          = h;
        num_perm_low  = num_perm_l;
        num_perm_high = num_perm_h;
      }
      else  //  if (num_roots != num_prem_h - num_perm_l)
        reduced = 0;
    }
  
  return reduced;
}

//  int ROOT1 :: contract(const bigrational& tol)
//    contracts *this s.t. [low, high] is no larger than tol.

int ROOT1 :: contract(const bigrational& tol)
{
  assert(tol > 0);
  
  //  If high - low < tol then nothing to do.
  //  Otherwise, contract *this.
  
  if (high - low > tol)
  {
    if (ok_float)
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
    
    while (high - low > tol)
      halve();
  }
  
  return 0;
}

//  int ROOT1 :: shrink(const bigrational& fac)
//    shrinks *this s.t. [low, high] becomes smaller by fac.
//    e.g. fac == 1/10 => [low, high] becomes 1/10.

int ROOT1 :: shrink(const bigrational& fac)
{
//  assert(sgn(fac) > 0);
  
  return contract(fac * (high - low));
}

int sort_fp_pts(const unsigned long l, double* const D)
{
  unsigned long i, j;
  unsigned long i_l;
  double        d_l, d;
  
  for (i = 0; i < l - 1; i++)
  {
    d_l = D[i_l = i];
    
    for (j = i + 1; j < l; j++)
      if (D[j] < d_l)
        d_l = D[i_l = j];
    
    if (i_l != i)
    {
      d      = D[i];
      D[i]   = D[i_l];
      D[i_l] = d;
    }
  }
  
  return 0;
}

//  unsigned long ROOT1 :: isolate_roots(ROOT1*& R, const bigrational& tol)
//    isolates the roots of *poly in [low, high] and store them in R.
//    each R[i] contains 1 and only 1 root of *poly, and
//    (R[i].low, R[i].high) is no larger than tol.

unsigned long ROOT1 :: isolate_roots(ROOT1*& R, const bigrational& tol)
{
#ifdef _EXPERIMENT
  num_root1_isolate++;
#endif
  
  if (num_roots == 0)
    R = 0;
  else  //  if (num_roots > 0)
  {
    R = new ROOT1 [num_roots];
    
    unsigned long num_fp_roots;
    double*       fp_roots;
    
    //  Compute fp approximations for the roots of *this.
    
    num_fp_roots =
      poly->as_FLOATPOLY().gen_fp_roots(low.as_double(), high.as_double(),
                                        fp_roots);
    
    if (num_roots == 1)
    //  If num_roots == 1 then contract *this s.t.
    //    [low, high] is no longer than tol.
    {
      //  Set float_est since it might be used to contract *this.
      
      if (num_fp_roots == 1)
      {
        float_est = fp_roots[0];
        ok_float  = 1;
      }
      else  //  if (num_fp_roots != 1)
      {
        float_est = 0.0;
        ok_float  = 0;
      }
      
      if (tol > 0)
        contract(tol);
      
      R[0] = ROOT1(*this);
    }
    else  //  if (num_roots > 1)
    //  If num_roots > 1 then isolate the roots of *this.
    {
      unsigned long i, j;
      int           use_exact;
      
      use_exact = 0;
      
      if (num_fp_roots != num_roots)
      //  Fp estimations do not work. Use exact arithmetic.
        use_exact = 1;
      else  //  if (num_fp_roots == num_roots)
      {
        double         d_left;
        double         dist;
        double         d_tol;
        unsigned long* num_bits;
        
        //  Sort fp_roots in increasing order.
        
        sort_fp_pts(num_fp_roots, fp_roots);
        
        //  Set dist to be the distance
        //    from fp_roots[i] to the nearest other roots or endpoints.
        
        i        = 0;
        d_left   = low.as_double();
        dist     = 1.0;
        d_tol    = tol.as_double();
        num_bits = new unsigned long [num_roots];
        
        while (!use_exact && i < num_roots - 1)
        {
          if (fp_roots[i] - d_left < fp_roots[i + 1] - fp_roots[i])
            dist = (fp_roots[i] - d_left) / 2.01;
          else  //  if (fp_roots[i] - d_left >= fp_roots[i + 1] - fp_roots[i])
            dist = (fp_roots[i + 1] - fp_roots[i]) / 2.01;
          
          if (tol > 0 && dist > d_tol)
            dist = d_tol;
          
          if (dist == 0.0)
            use_exact = 1;
          else  //  if (dist > 0.0)
          //  Compute an upper bound on dist in bits.
          {
            num_bits[i] = 0;
            
            while (dist <= 1.0)
            {
              dist *= 2.0;
              num_bits[i]++;
            }
            
            //  Prepare for the next iteration.
            
            d_left = fp_roots[i++];
          }
        }
        
        if (!use_exact)
        {
          if (fp_roots[num_roots - 1] - d_left
              <
              high.as_double() - fp_roots[num_roots - 1])
            dist = (fp_roots[num_roots - 1] - d_left) / 2.01;
          else
            dist = (high.as_double() - fp_roots[num_roots - 1]) / 2.01;
          
          if (tol > 0 && dist > d_tol)
            dist = d_tol;
          
          if (dist == 0.0)
            use_exact = 1;
          else  //  if (dist > 0.0)
          //  Compute an upper bound on dist in bits.
          {
            num_bits[i] = 0;
            
            while (dist <= 1.0)
            {
              dist *= 2.0;
              num_bits[i]++;
            }
            
            //  Use rational arithmetic to double-check whether or not
            //    the intervals containing fp_roots[i] of length num_bits[i]
            //      are disjoint.
            
            for (i = 0; !use_exact && i < num_roots; i++)
            {
              bigrational r_l, r_h;
              
              double2bigrational_interval(fp_roots[i], num_bits[i], r_l, r_h);
              R[i] = ROOT1(*poly, r_l, r_h, fp_roots[i]);
              
              if (r_l < low)
              {
                R[i].cut(low);
                
                if (R[i].low < low)
                  use_exact = 1;
              }
              
              if (r_h > high)
              {
                R[i].cut(high);
                
                if (R[i].high > high)
                  use_exact = 1;
              }
              
              ok_float = 0;
              
              if (R[i].num_roots != 1)
                use_exact = 1;
            }
          }
        }
        
        delete [] num_bits;
      }
      
      if (use_exact)
      //  Fp estimations do not work. Use exact arithmetic.
      {
#ifdef _EXPERIMENT
        num_root1_exact_isolate++;
#endif
        
        bigrational   half;
        unsigned long num_perm_half, diff_num_perm_l, diff_num_perm_h;
        ROOT1         r0;
        unsigned long num_roots_r0;
        ROOT1*        R0;
        
        //  Bisect [low, high] and recurse.
        
        half = (low + high) / 2;
        
        num_perm_half   = poly->num_Sturm_seq_perm(half);
        diff_num_perm_l = num_perm_half - num_perm_low;
        diff_num_perm_h = num_perm_high - num_perm_half;
        
        i = 0;
        
        if (diff_num_perm_l == 1)
        //  if there is 1 and only 1 root of *poly in [low, half]
        {
          R[i] = ROOT1(*poly, low, half, num_perm_low, num_perm_half);
          
          if (tol > 0)
            R[i].contract(tol);
          
          ok_float = 0;
          i++;
        }
        else if (diff_num_perm_l > 1)
        //  If there are more than 1 roots of *poly in [low, half] then
        //    recurse.
        {
          r0           = ROOT1(*poly, low, half, num_perm_low, num_perm_half);
          num_roots_r0 = r0.isolate_roots(R0, tol);
          
          assert(r0.num_roots == num_roots_r0);
          
          for (j = 0; j < num_roots_r0; j++)
            R[i++] = R0[j];
          
          if (num_roots_r0)
            delete [] R0;
        }
        
        if (diff_num_perm_h == 1)
        //  if there is 1 and only 1 root of *poly in [half, high]
        {
          R[i] = ROOT1(*poly, half, high, num_perm_half, num_perm_high);
          
          if (tol > 0)
            R[i].contract(tol);
          
          ok_float = 0;
          i++;
        }
        else if (diff_num_perm_h > 1)
        //  If there are more than 1 roots of *poly in [half, high] then
        //    recurse.
        {
          r0           =
            ROOT1(*poly, half, high, num_perm_half, num_perm_high);
          num_roots_r0 = r0.isolate_roots(R0, tol);
          
          assert(r0.num_roots == num_roots_r0);
          
          for (j = 0; j < num_roots_r0; j++)
            R[i++] = R0[j];
          
          if (num_roots_r0)
            delete [] R0;
        }
        
        assert(i == num_roots);
      }
    }
    
    if (num_fp_roots)
      delete [] fp_roots;
  }
  
  return num_roots;
}

