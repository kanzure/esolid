#include <bigint.h>

bigint :: bigint()
{
  mpz_init(rep);
}

bigint :: bigint(const unsigned int u)
{
  mpz_init_set_ui(rep, u);
}

bigint :: bigint(const signed int s)
{
  mpz_init_set_si(rep, s);
}

bigint :: bigint(const unsigned long u)
{
  mpz_init_set_ui(rep, u);
}

bigint :: bigint(const signed long s)
{
  mpz_init_set_si(rep, s);
}

bigint :: bigint(const double d)
{
  mpz_init_set_d(rep, d);
}

bigint :: bigint(const char* str)
{
  mpz_init_set_str(rep, str, 0);
}

bigint :: bigint(const mpz_t z)
{
  mpz_init_set(rep, z);
}

bigint :: bigint(const bigint& I)
{
  mpz_init_set(rep, I.rep);
}

bigint& bigint :: operator =(const bigint& I)
{
  mpz_set(rep, I.rep);
  
  return *this;
}

bigint :: ~bigint()
{
  mpz_clear(rep);
}

//bigint Rand()
//{
//  bigint x;
//  
//  mpz_urandomb(x.rep, rand_state, rand_length);
//  
//  return x;
//}

unsigned long bigint :: len() const
{
  return mpz_sizeinbase(rep, 2);
}

int bigint :: fits_long() const
{
  return mpz_fits_slong_p(rep);
}

long bigint :: as_long() const
{
  return mpz_get_si(rep);
}

int bigint :: fits_unsigned_long() const
{
  return mpz_fits_ulong_p(rep);
}

unsigned long bigint :: as_unsigned_long() const
{
  return mpz_get_ui(rep);
}

double bigint :: frexp(long* x) const
{
  return mpz_get_d_2exp(x, rep);
}

double frexp(const bigint& x, long* y)
{
  return x.frexp(y);
}

long lb4lgAbs(const bigint& x)
{
  long y;
  
  x.abs().frexp(&y);
  
  return y - 1;
}

long ub4lgAbs(const bigint& x)
{
  double y;
  long   z;
  
  y = x.abs().frexp(&z);
  
  if (y == .5)
    return z - 1;
  else
    return z;
}

bigint bigint :: add(const bigint& x) const
{
  bigint y;
  
  mpz_add(y.rep, rep, x.rep);
  
  return y;
}

bigint operator +(const bigint& x, const bigint& y)
{
  return x.add(y);
}

bigint bigint :: sub(const bigint& x) const
{
  bigint y;
  
  mpz_sub(y.rep, rep, x.rep);
  
  return y;
}

bigint operator -(const bigint& x, const bigint& y)
{
  return x.sub(y);
}

bigint bigint :: mul(const bigint& x) const
{
  bigint y;
  
  mpz_mul(y.rep, rep, x.rep);
  
  return y;
}

bigint operator *(const bigint& x, const bigint& y)
{
  return x.mul(y);
}

bigint bigint :: div(const bigint& x) const
{
  bigint y;
  
  mpz_tdiv_q(y.rep, rep, x.rep);
  
  return y;
}

bigint operator /(const bigint& x, const bigint& y)
{
  return x.div(y);
}

void bigint :: inc()
{
  mpz_add_ui(rep, rep, 1);
}

bigint& bigint :: operator ++()  //  prefix
{
  inc();
  
  return *this;
}

bigint bigint :: operator ++(int)  //  postfix
{
  bigint x(*this);
  
  inc();
  
  return x;
}

void bigint :: dec()
{
  mpz_sub_ui(rep, rep, 1);
}

bigint& bigint :: operator --()  //  prefix
{
  dec();
  
  return *this;
}

bigint bigint :: operator --(int)  //  postfix
{
  bigint x(*this);
  
  dec();
  
  return x;
}

bigint bigint :: lshift(const unsigned long x) const
{
  bigint y;
  
  mpz_mul_2exp(y.rep, rep, x);
  
  return y;
}

bigint operator <<(const bigint& x, const unsigned long y)
{
  return x.lshift(y);
}

bigint bigint :: rshift(const unsigned long x) const
{
  bigint y;
  
  mpz_fdiv_q_2exp(y.rep, rep, x);
  
  return y;
}

bigint operator >>(const bigint& x, const unsigned long y)
{
  return x.rshift(y);
}

bigint bigint :: neg() const
{
  bigint y;
  
  mpz_neg(y.rep, rep);
  
  return y;
}

bigint operator -(const bigint& x)
{
  return x.neg();
}

bigint bigint :: abs() const
{
  bigint y;
  
  mpz_abs(y.rep, rep);
  
  return y;
}

bigint abs(const bigint& x)
{
  return x.abs();
}

bigint bigint :: root(const unsigned long x) const
{
  bigint y;
  
  mpz_root(y.rep, rep, x);
  
  return y;
}

bigint root(const bigint& x, const unsigned long y)
{
  return x.root(y);
}

bigint bigint :: sqrt() const
{
  bigint x;
  
  mpz_sqrt(x.rep, rep);
  
  return x;
}

bigint sqrt(const bigint& x)
{
  return x.sqrt();
}

int bigint :: cmp(const bigint& x) const
{
  return mpz_cmp(rep, x.rep);
}

int operator ==(const bigint& x, const bigint& y)
{
  return x.cmp(y) == 0;
}

int operator !=(const bigint& x, const bigint& y)
{
  return x.cmp(y) != 0;
}

int operator <(const bigint& x, const bigint& y)
{
  return x.cmp(y) < 0;
}

int operator <=(const bigint& x, const bigint& y)
{
  return x.cmp(y) <= 0;
}

int operator >(const bigint& x, const bigint& y)
{
  return x.cmp(y) > 0;
}

int operator >=(const bigint& x, const bigint& y)
{
  return x.cmp(y) >= 0;
}

int bigint :: sgn() const
{
  return mpz_sgn(rep);
}

int sgn(const bigint& x)
{
  return x.sgn();
}

int operator &&(const bigint& x, const bigint& y)
{
  return x.sgn() && y.sgn();
}

int operator ||(const bigint& x, const bigint& y)
{
  return x.sgn() || y.sgn();
}

int operator !(const bigint& x)
{
  return x.sgn() == 0;
}

bigint bigint :: AND(const bigint& x) const
{
  bigint y;
  
  mpz_and(y.rep, rep, x.rep);
  
  return y;
}

bigint operator &(const bigint& x, const bigint& y)
{
  return x.AND(y);
}

bigint bigint :: IOR(const bigint& x) const
{
  bigint y;
  
  mpz_ior(y.rep, rep, x.rep);
  
  return y;
}

bigint operator |(const bigint& x, const bigint& y)
{
  return x.IOR(y);
}

bigint bigint :: XOR(const bigint& x) const
{
  bigint y;
  
  mpz_xor(y.rep, rep, x.rep);
  
  return y;
}

bigint operator ^(const bigint& x, const bigint& y)
{
  return x.XOR(y);
}

bigint bigint :: gcd(const bigint& x) const
{
  bigint y;
  
  mpz_gcd(y.rep, rep, x.rep);
  
  return y;
}

bigint gcd(const bigint& x, const bigint& y)
{
  return x.gcd(y);
}

bigint bigint :: lcm(const bigint& x) const
{
  bigint y;
  
  mpz_lcm(y.rep, rep, x.rep);
  
  return y;
}

bigint lcm(const bigint& x, const bigint& y)
{
  return x.lcm(y);
}

bigint& bigint :: operator +=(const bigint& x)
{
  mpz_add(rep, rep, x.rep);
  
  return *this;
}

bigint& bigint :: operator -=(const bigint& x)
{
  mpz_sub(rep, rep, x.rep);
  
  return *this;
}

bigint& bigint :: operator *=(const bigint& x)
{
  mpz_mul(rep, rep, x.rep);
  
  return *this;
}

bigint& bigint :: operator /=(const bigint& x)
{
  mpz_tdiv_q(rep, rep, x.rep);
  
  return *this;
}

bigint& bigint :: operator <<=(const unsigned long x)
{
  mpz_mul_2exp(rep, rep, x);
  
  return *this;
}

bigint& bigint :: operator >>=(const unsigned long x)
{
  mpz_div_2exp(rep, rep, x);
  
  return *this;
}

bigint& bigint :: operator &=(const bigint& x)
{ 
  mpz_and(rep, rep, x.rep);
  
  return *this;
}

bigint& bigint :: operator |=(const bigint& x)
{
  mpz_ior(rep, rep, x.rep);
  
  return *this;
}

bigint& bigint :: operator ^=(const bigint& x)
{
  mpz_xor(rep, rep, x.rep);
  
  return *this;
}

int bigint :: bit(const unsigned long x) const
{
  return mpz_tstbit(rep, x);
}

//ostream& operator <<(ostream& o, const bigint& x)
//{
//  char* s = new char [256];
//  
//  mpz_get_str(s, 10, x.rep);
//  o << s;
//  
//  delete [] s;
//  
//  return o;
//}

ostream& operator <<(ostream& o, const bigint& x)
{
  char* s;
  
  if (s = mpz_get_str(NULL, 10, x.rep))
  {
    o << s;
    
    delete [] s;
  }
  
  return o;
}

istream& operator >>(istream& i, bigint& x)
{
  char c;
  int  s;
  
  s = 1;
  x = 0;
  
  i >> c;
  
  if (c == '-')
  {
    s = - 1;
    i >> c;
  }
  
  while (c >= '0' && c <= '9')
  {
    x = x * 10 + (long)(c - '0');
    i.get(c);
  }
  
  if (s < 0)
    x = - x;
  
  return i;
}

