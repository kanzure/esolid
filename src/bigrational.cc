#include <bigrational.h>

bigrational :: bigrational()
{
  mpq_init(rep);
}

bigrational :: bigrational(const unsigned int u)
{
  mpq_init(rep);
  mpq_set_ui(rep, u, 1);
}

bigrational :: bigrational(const signed int s)
{
  mpq_init(rep);
  mpq_set_si(rep, s, 1);
}

bigrational :: bigrational(const unsigned long u)
{
  mpq_init(rep);
  mpq_set_ui(rep, u, 1);
}

bigrational :: bigrational(const signed long s)
{
  mpq_init(rep);
  mpq_set_si(rep, s, 1);
}

bigrational :: bigrational(const double d)
{
  mpq_init(rep);
  mpq_set_d(rep, d);
  mpq_canonicalize(rep);
}

bigrational :: bigrational(const bigint& I)
{
  mpq_init(rep);
  mpq_set_z(rep, I.rep);
}

bigrational :: bigrational(const bigint& N, const bigint& D)
{
  mpq_init(rep);
  mpz_set(mpq_numref(rep), N.rep);
  mpz_set(mpq_denref(rep), D.rep);
  mpq_canonicalize(rep);
}

bigrational :: bigrational(const bigrational& Q)
{
  mpq_init(rep);
  mpq_set(rep, Q.rep);
}

bigrational& bigrational :: operator =(const bigrational& Q)
{
  mpq_set(rep, Q.rep);
  
  return *this;
}

bigrational :: ~bigrational()
{
  mpq_clear(rep);
}

//static bigint DRand_den = bigint(1) << rand_length;
//
//bigrational DRand()
//{
//  return bigrational(Rand(), DRand_den);
//}

unsigned long bigrational :: len() const
{
  unsigned long l_n, l_d;
  
  l_n = mpz_sizeinbase(mpq_numref(rep), 2);
  l_d = mpz_sizeinbase(mpq_denref(rep), 2);
  
  return l_n + l_d;
}

double bigrational :: as_double() const
{
  return mpq_get_d(rep);
}

int bigrational :: fits_bigint() const
{
  return mpz_cmp_ui(mpq_denref(rep), 1) == 0;
}

bigint bigrational :: as_bigint() const
{
  bigint x;
  
  mpz_divexact(x.rep, mpq_numref(rep), mpq_denref(rep));
  
  return x;
}

bigint bigrational :: num() const
{
  return bigint(mpq_numref(rep));
}

bigint bigrational :: den() const
{
  return bigint(mpq_denref(rep));
}

bigint num(const bigrational& x)
{
  return x.num();
}

bigint den(const bigrational& x)
{
  return x.den();
}

bigint bigrational :: ceil() const
{
  bigint x;
  
  mpz_cdiv_q(x.rep, mpq_numref(rep), mpq_denref(rep));
  
  return x;
}

bigint ceil(const bigrational& x)
{
  return x.ceil();
}

bigint bigrational :: floor() const
{
  bigint x;
  
  mpz_fdiv_q(x.rep, mpq_numref(rep), mpq_denref(rep));
  
  return x;
}

bigint floor(const bigrational& x)
{
  return x.floor();
}

bigint bigrational :: trunc() const
{
  bigint x;
  
  mpz_tdiv_q(x.rep, mpq_numref(rep), mpq_denref(rep));
  
  return x;
}

bigint trunc(const bigrational& x)
{
  return x.trunc();
}

long bigrational :: lgAbs() const
{
  long ld, ln;
  
  num().abs().frexp(&ln);
  den().abs().frexp(&ld);
  
  return ln - ld;
}

long lb4lgAbs(const bigrational& x)
{
  return x.lgAbs() - 1;
}

long ub4lgAbs(const bigrational& x)
{
  return x.lgAbs() + 1;
}

bigrational bigrational :: add(const bigrational& x) const
{
  bigrational y;
  
  mpq_add(y.rep, rep, x.rep);
  mpq_canonicalize(y.rep);
  
  return y;
}

bigrational operator +(const bigrational& x, const bigrational& y)
{
  return x.add(y);
}

bigrational bigrational :: sub(const bigrational& x) const
{
  bigrational y;
  
  mpq_sub(y.rep, rep, x.rep);
  mpq_canonicalize(y.rep);
  
  return y;
}

bigrational operator -(const bigrational& x, const bigrational& y)
{
  return x.sub(y);
}

bigrational bigrational :: mul(const bigrational& x) const
{
  bigrational y;
  
  mpq_mul(y.rep, rep, x.rep);
  mpq_canonicalize(y.rep);
  
  return y;
}

bigrational operator *(const bigrational& x, const bigrational& y)
{
  return x.mul(y);
}

bigrational bigrational :: div(const bigrational& x) const
{
  bigrational y;
  
  mpq_div(y.rep, rep, x.rep);
  mpq_canonicalize(y.rep);
  
  return y;
}

bigrational operator /(const bigrational& x, const bigrational& y)
{
  return x.div(y);
}

bigrational bigrational :: neg() const
{
  bigrational x;
  
  mpq_neg(x.rep, rep);
  
  return x;
}

bigrational operator -(const bigrational& x)
{
  return x.neg();
}

bigrational bigrational :: abs() const
{
  bigrational x;
  
  mpq_abs(x.rep, rep);
  
  return x;
}

bigrational abs(const bigrational& x)
{
  return x.abs();
}

bigrational bigrational :: inv() const
{
  bigrational x;
  
  mpq_inv(x.rep, rep);
  
  return x;
}

bigrational inv(const bigrational& x)
{
  return x.inv();
}

int bigrational :: equal(const bigrational& x) const
{
  return mpq_equal(rep, x.rep);
}

int operator ==(const bigrational& x, const bigrational& y)
{
  return x.equal(y) != 0;
}

int operator !=(const bigrational& x, const bigrational& y)
{
  return x.equal(y) == 0;
}

int bigrational :: cmp(const bigrational& x) const
{
  return mpq_cmp(rep, x.rep);
}

int operator <(const bigrational& x, const bigrational& y)
{
  return x.cmp(y) < 0;
}

int operator <=(const bigrational& x, const bigrational& y)
{
  return x.cmp(y) <= 0;
}

int operator >(const bigrational& x, const bigrational& y)
{
  return x.cmp(y) > 0;
}

int operator >=(const bigrational& x, const bigrational& y)
{
  return x.cmp(y) >= 0;
}

int bigrational :: sgn() const
{
  return mpq_sgn(rep);
}

int sgn(const bigrational& x)
{
  return x.sgn();
}

int operator &&(const bigrational& x, const bigrational& y)
{
  return x.sgn() && y.sgn();
}

int operator ||(const bigrational& x, const bigrational& y)
{
  return x.sgn() || y.sgn();
}

int operator !(const bigrational& x)
{
  return x.sgn() == 0;
}

bigrational& bigrational :: operator +=(const bigrational& x)
{
  mpq_add(rep, rep, x.rep);
  mpq_canonicalize(rep);
  
  return *this;
}

bigrational& bigrational :: operator -=(const bigrational& x)
{
  mpq_sub(rep, rep, x.rep);
  mpq_canonicalize(rep);
  
  return *this;
}

bigrational& bigrational :: operator *=(const bigrational& x)
{
  mpq_mul(rep, rep, x.rep);
  mpq_canonicalize(rep);
  
  return *this;
}

bigrational& bigrational :: operator /=(const bigrational& x)
{
  mpq_div(rep, rep, x.rep);
  mpq_canonicalize(rep);
  
  return *this;
}

//ostream& operator <<(ostream& o, const bigrational& x)
//{
//  char* s = new char [256];
//  mpz_t n, d;
//  
//  mpz_init_set(n, mpq_numref(x.rep));
//  mpz_init_set(d, mpq_denref(x.rep));
//  
//  mpz_get_str(s, 10, n);
//  o << s;
//  
//  if (mpz_cmp_ui(d, 1))
//  {
//    mpz_get_str(s, 10, d);
//    o << "/" << s;
//  }
//  
//  mpz_clear(n);
//  mpz_clear(d);
//  delete [] s;
//  
//  return o;
//}

ostream& operator <<(ostream& o, const bigrational& x)
{
  char* s;
  
  if (s = mpq_get_str(NULL, 10, x.rep))
  {
    o << s;
    delete [] s;
  }
  
  return o;
}

istream& operator >>(istream& i, bigrational& x)
{
  char   c;
  int    s;
  bigint n, d;
  
  s = 1;
  n = d = 0;
  
  i >> c;
  
  if (c == '-')
  {
    s = - 1;
    i >> c;
  }
  
  while (c >= '0' && c <= '9')
  {
    n = n * 10 + (long)(c - '0');
    i.get(c);
  }
  
  while (c == ' ' || c == '\t')
    i.get(c);
  
  assert(c == '/');
  
  i >> d;
  
  x = s * bigrational(n, d);
  
  return i;
}

