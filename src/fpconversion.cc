//  file:    fpconversion.cc
//  update:  10/14/02

#include <fpconversion.h>

#if (defined(_Linux_i386_))
#define LITTLE_ENDIAN_DOUBLE
#elif ((defined(_SunOS_)) || (defined(_IRIX_)))
#define BIG_ENDIAN_DOUBLE
#endif

#if (defined(BIG_ENDIAN_DOUBLE))
#define HI_HALF(x) (((unsigned int*)(&(x)))[0])
#define LO_HALF(x) (((unsigned int*)(&(x)))[1])
#elif (defined(LITTLE_ENDIAN_DOUBLE))
#define HI_HALF(x) (((unsigned int*)(&(x)))[1])
#define LO_HALF(x) (((unsigned int*)(&(x)))[0])
#endif

//  get_bits(): see K&R, 2nd ed., p. 49
static unsigned int get_bits(const unsigned int x,
                             const unsigned int p,
                             const unsigned int n)
{
  return (x >> (p + 1 - n)) & ~(~0 << n);
}

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

int double2bigrational_interval(const double x,
                                unsigned int n,
                                bigrational& l, bigrational& r)
{
  bigint num_l, den_l, num_r, den_r;  //  numerators & denominators of l & r
  
  unsigned int num_bits;  //  number of bits of [l, r]
  
  unsigned int x_h, x_l;  //  upper && lower words of x
  
  unsigned int s, e, f_h, f_l;  //  contents of the IEEE double fields
                                //    sign, exponent,
                                //    upper && lower words of mantissa
                                //  note that f takes up 53 bits
                                //  and must be split
  int e_signed;  //  the exponent, interpreted
  
  int          num_use_bits;          //  number of bits of `f' we' ll use
  unsigned int num_use_bits_h;        //  portions of `num_use_bits'
  unsigned int num_use_bits_l;        //    that applies to f_h & f_l
  unsigned int num_use_bits_further;  //  portion of `num_use_bits'
                                      //    beyond the known significand
                                      //  -- >0 is returned if this is positive
  
  //  Ensure that the result will have at least one bit of accuracy.
  
  if ((num_bits = n) < 1)
    num_bits = 1;
  
  x_h = HI_HALF(x);
  x_l = LO_HALF(x);
  
  //  Parse x
  
  s   = get_bits(x_h, 31, 1);
  e   = get_bits(x_h, 30, 11);
  f_h = get_bits(x_h, 19, 20);
  f_l = x_l;
  
  if (e == 2047)  //  x is NaN, Inf, ...
    return - 1;
  else
  {
    if (!e && !f_h && !f_l)  //  x == 0.0
    {
      if (num_bits > 53)
      {
        num_use_bits_further = num_bits - 53;
        num_bits             = 53;
      }
      else
        num_use_bits_further = 0;
      
      den_l = 1;
      den_l <<= num_bits + 1;
      den_r = den_l;
      num_l = - 1;
      num_r = 1;
      
      s = 0;
    }
    else  //  x != 0.0
    {
      //  Interpret e.
      
      if (e)  //  x is normalized
      {
        f_h      |= (unsigned int)1 << 20;  //  Pad the hidden MSB of f.
        e_signed  = e - 1023;
      }
      else  //  x is "subnormal"
        e_signed = - 1022;
      
      //  Compute l & r.
      
      den_l = 1;
      den_l <<= num_bits;
      den_r = 1;
      den_r <<= num_bits;
      
      //  We'll use num_use_bits leading bits from the significand
      //  & pad num_use_bits_further many 0's to right.
      
      num_use_bits = num_bits + e_signed + 1;
      
      if (num_use_bits <= 53)
        num_use_bits_further = 0;
      else
      {
        num_use_bits_further = num_use_bits - 53;
        num_use_bits         = 53;
      }
      
      if (num_use_bits < 1)  //  No significant bits are needed.
      {
        num_l = 0;
        num_r = 1;
      }
      else  //  We'll use some (posiibly all) bits of significand.
      {
        num_use_bits_h = num_use_bits <= 21 ? num_use_bits : 21;
        num_use_bits_l = num_use_bits <= 21 ? 0 : num_use_bits - 21;
        
        num_l = get_bits(f_h, 20, num_use_bits_h);
        num_l <<= num_use_bits_l;
        
        if (num_use_bits_l != 32)
          num_l += get_bits(f_l, 31, num_use_bits_l);
        else  //  f_l must be getting cast to a signed int. Fix it.
          num_l += f_l;
        
        num_r = num_l + 1;
        num_l <<= num_use_bits_further;
        num_r <<= num_use_bits_further;
      }
    }
    
    //  Take care of signs.
    
    bigrational l0(num_l, den_l);
    bigrational r0(num_r, den_r);
    
    if (!s)  //  s is 0 iff x >= 0.0
    {
      l = l0;
      r = r0;
    }
    else  //  s is 1 iff x < 0.0
    {
      l = - r0;
      r = - l0;
    }
    
    return num_use_bits_further ? 1 : 0;
  }
}

bigrational as_bigrational(const float f)
{
  bigrational x, l, h;
  
  if (f == 0.0)
    x = 0;
  else
  {
    double2bigrational_interval(f, 23, l, h);
    
    if (f < 0.0)
      x = h;
    else  //  if (f > 0.0)
      x = l;
  }
  
  return x;
}

