//  file:    kratpoly.cc
//  update:  10/04/02, 04/15/04

#ifdef _EXPERIMENT
#include <counter.h>
#endif

#include <kratpoly.h>

#include <fpfilter.h>
#include <VDM.h>

//  K_RATPOLY :: K_RATPOLY()
//    constructs a zero polynomial.

K_RATPOLY :: K_RATPOLY()
{
  num_vars = 0;
  deg      = 0;
  
  coeffs    = new bigrational [num_coeffs = 1];
  coeffs[0] = 0;
  
  Sturm_seq = 0;
  
  ref_count = 0;
}

//  K_RATPOLY :: K_RATPOLY(const unsigned long nv, const long* d)
//    constructs a polynomial
//      of nv variables of degree d[0], d[1], ..., d[nv - 1]
//      with anonymous coefficients.

K_RATPOLY :: K_RATPOLY(const unsigned long nv, const long* const d)
{
  unsigned long i;
  
  if ((num_vars = nv) > 0)
    deg = new long [num_vars];
  else  //  if (num_vars == 0)
    deg = 0;
  
  num_coeffs = 1;
  
  for (i = 0; i < num_vars; i++)
  {
    deg[i]      = d[i];
    num_coeffs *= d[i] + 1;
  }
  
  coeffs = new bigrational [num_coeffs];
  
  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = 0;
  
  Sturm_seq = 0;
  
  ref_count = 0;
}

//  K_RATPOLY :: K_RATPOLY(const unsigned long nv, const long* const d,
//                         const unsigned long nc, const bigrational* const c)
//    constructs a polynomial
//      of nv variables of degree d[0], d[1], ..., d[nv - 1]
//      with coefficients c[0], c[1], ..., c[nc - 1].

K_RATPOLY :: K_RATPOLY(const unsigned long nv, const long* const d,
                       const unsigned long nc, const bigrational* const c)
{
  unsigned long i;
  
  if ((num_vars = nv) > 0)
    deg = new long [num_vars];
  else  //  if (num_vars == 0)
    deg = 0;
  
  num_coeffs = 1;
  
  for (i = 0; i < num_vars; i++)
  {
    deg[i]      = d[i];
    num_coeffs *= d[i] + 1;
  }
  
  assert(nc == num_coeffs);
  
  coeffs = new bigrational [num_coeffs];
  
  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = c[i];
  
  Sturm_seq = 0;
  
  ref_count = 0;
  
  reduce_deg();
}

//  K_RATPOLY :: K_RATPOLY(const unsigned long nv,
//                         const unsigned long i, const bigrational b)
//    constructs a polynomial X_i - b of nv variables.

K_RATPOLY :: K_RATPOLY(const unsigned long nv,
                       const unsigned long i, const bigrational& b)
{
  assert(i < nv);  //  i < nv => nv > 0
  
  unsigned long j;
  
  deg = new long [num_vars = nv];
  
  for (j = 0; j < i; j++)
    deg[j] = 0;
  
  deg[i] = 1;
  
  for (j = i + 1; j < num_vars; j++)
    deg[j] = 0;
  
  coeffs    = new bigrational [num_coeffs = 2];
  coeffs[0] = 1;
  coeffs[1] = - b;
  
  Sturm_seq = 0;
  
  ref_count = 0;
}

//  K_RATPOLY :: K_RATPOLY(const K_RATPOLY& P)
//    the copy constructor

K_RATPOLY :: K_RATPOLY(const K_RATPOLY& P)
{
  unsigned long i;
  
  if ((num_vars = P.num_vars) > 0)
  {
    deg = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      deg[i] = P.deg[i];
  }
  else  //  if (num_vars == 0)
    deg = 0;
  
  coeffs = new bigrational [num_coeffs = P.num_coeffs];
  
  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = P.coeffs[i];
  
  if (Sturm_seq = P.Sturm_seq)
    Sturm_seq->ref_count++;
  
  ref_count = 0;
}

//  K_RATPOLY& K_RATPOLY :: operator =(const K_RATPOLY& P)
//    the assignment operator

K_RATPOLY& K_RATPOLY :: operator =(const K_RATPOLY& P)
{
  if (this != &P)
  {
    unsigned long i;
    
    if (deg)
      delete [] deg;
    
    if (coeffs)
      delete [] coeffs;
    
    if (Sturm_seq && !--Sturm_seq->ref_count)
      delete Sturm_seq;
    
    if ((num_vars = P.num_vars) > 0)
    {
      deg = new long [num_vars];
      
      for (i = 0; i < num_vars; i++)
        deg[i] = P.deg[i];
    }
    else  //  if (num_vars == 0)
      deg = 0;
    
    coeffs = new bigrational [num_coeffs = P.num_coeffs];
    
    for (i = 0; i < num_coeffs; i++)
      coeffs[i] = P.coeffs[i];
    
    if (Sturm_seq = P.Sturm_seq)
      Sturm_seq->ref_count++;
  }
  
  return *this;
}

//  K_RATPOLY :: ~K_RATPOLY()
//    the destructor

K_RATPOLY :: ~K_RATPOLY()
{
  if (deg)
    delete [] deg;
  
  if (coeffs)
    delete [] coeffs;
  
  if (Sturm_seq && !--Sturm_seq->ref_count)
    delete Sturm_seq;
}

ostream& K_RATPOLY :: output(ostream& o) const
{
  unsigned long i, j, k;
  long*         p;
  
  for (i = 0; i < num_coeffs && !sgn(coeffs[i]); i++)
    ;
  
  if (i < num_coeffs)
  {
    //  *this is NOT a zero polynomial.
    
    for (j = 0; j < num_coeffs; j++)
      if (sgn(coeffs[j]))
      {
        assert(p = index_to_powers(j));
        
        o << coeffs[j] << " X^(";
        
        for (k = 0; k < num_vars - 1; k++)
          o << p[k] << ',';
        
        o << p[num_vars - 1] << ')' << endl;
        delete [] p;
      }
    
    o << flush;
  }
  else  //  if (i == num_coeffs)
  {
    //  *this is a zero polynomial.
    
    o << "0 (";
    
    for (k = 0; k < num_vars - 1; k++)
      o << "0,";
    
    o << "0)" << endl;
  }
  
  return o;
}

ostream& operator <<(ostream& o, const K_RATPOLY& P)
{
  return P.output(o);
}

//  unsigned long K_RATPOLY :: get_num_vars() const
//    returns the number of variables of *this.

unsigned long K_RATPOLY :: get_num_vars() const
{
  return num_vars;
}

//  long* K_RATPOLY :: index_to_powers(const unsigned long i) const
//    returns p = (p[0], p1[1], ..., p[num_vars - 1]) s.t.
//      coeffs[i] is the coefficient of the monomial X^p of *this.

long* K_RATPOLY :: index_to_powers(const unsigned long i) const
{
  assert(i < num_coeffs);
  
  unsigned long j;
  long          k;
  long*         p;
  
  if (num_vars > 0)
  {
    j = i;
    p = new long [num_vars];
    
    for (k = num_vars - 1; k >= 0; k--)
    {
      p[k]  = deg[k] - j % (deg[k] + 1);
      j    /= deg[k] + 1;
    }
  }
  else  //  if (num_vars == 0)
    p = 0;
  
  return p;
}

//  unsigned long K_RATPOLY :: index_to_total_deg(const unsigned long i) const
//    returns t s.t.
//      coeffs[i] is the coefficient of the monomial of *this
//                                        of total degree t.

unsigned long K_RATPOLY :: index_to_total_deg(const unsigned long i) const
{
  assert(i < num_coeffs);
  
  unsigned long j;
  long          k;
  unsigned long t;
  
  j = i;
  t = 0;
  
  for (k = num_vars - 1; k >= 0; k--)
  {
    t += deg[k] - j % (deg[k] + 1);
    j /= deg[k] + 1;
  }
  
  return t;
}

//  unsigned long K_RATPOLY :: get_total_deg() const
//    returns the total degree of *this.

unsigned long K_RATPOLY :: get_total_deg() const
{
  unsigned long i, ti;
  unsigned long t;
  
  t = 0;
  
  for (i = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]) && (ti = index_to_total_deg(i)) > t)
      t = ti;
  
  return t;
}

//  unsigned long K_RATPOLY :: powers_to_index(const long* const p) const
//    returns i s.t.
//      coeffs[i] is the coefficient of the monomial X^p of *this.

unsigned long K_RATPOLY :: powers_to_index(const long* const p) const
{
  unsigned long j, k;
  unsigned long i;
  
  i = 0;
  
  if (p)
    for (j = 0, k = num_coeffs; j < num_vars; j++)
    {
      assert(p[j] <= deg[j]);
      
      k /= deg[j] + 1;
      i += k * (deg[j] - p[j]);
    }
  
  assert(i < num_coeffs);
  
  return i;
}

//  bigrational& K_RATPOLY :: get_coeff(const long* const p) const
//    returns the coefficient of the monomial X^p of *this.

bigrational& K_RATPOLY :: get_coeff(const long* const p) const
{
  return coeffs[powers_to_index(p)];
}

//  int K_RATPOLY :: reduce_deg()
//    reduces deg[0], deg[1], ..., deg[num_vars - 1]
//      to max. degrees in X_0, X_1, ..., X_{num_vars - 1}, resp.
//    returns 1 if some deg[i] is reduced, and
//            0 otherwise.

int K_RATPOLY :: reduce_deg()
{
  unsigned long i, j, k, l;
  long*         d;
  long*         p;
  unsigned long nc;
  bigrational*  c;
  int           reduced;
  
  reduced = 0;
  
  if (num_vars > 0)
  {
    d = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      d[i] = 0;
    
    for (i = 0; i < num_coeffs; i++)
      if (sgn(coeffs[i]) && (p = index_to_powers(i)))
      {
        for (j = 0; j < num_vars; j++)
          if (d[j] < p[j])
            d[j] = p[j];
        
        delete [] p;
      }
    
    i = 0;
    
    while (!reduced && i < num_vars)
      if (d[i] < deg[i])
        reduced = 1;
      else  //  if (d[i] >= deg[i])
        i++;
    
    if (reduced)
    {
      for (i = 0, nc = 1; i < num_vars; i++)
        nc *= d[i] + 1;
      
      c = new bigrational [nc];
      
      for (i = 0; i < num_coeffs; i++)
        if (sgn(coeffs[i]))
        {
          p = index_to_powers(i);
          
          for (j = k = 0, l = nc; k < num_vars; k++)
          {
            assert(p[k] <= d[k]);
            
            l /= d[k] + 1;
            j += l * (d[k] - p[k]);
          }
          
          c[j] = coeffs[i];
          delete [] p;  //  num_vars > 0 => p != 0
        }
      
      delete [] coeffs;
      num_coeffs = nc;
      coeffs     = c;
    }
    
    delete [] deg;  //  num_vars > 0 => deg != 0
    deg = d;
  }
  
  return reduced;
}

//  K_RATPOLY K_RATPOLY :: add_var(const unsigned long i) const
//    returns a polynomial of num_vars + 1 variables of degree
//      deg[0], deg[1], ..., deg[i - 1], 0, deg[i], ..., deg[num_vars - 1].

K_RATPOLY K_RATPOLY :: add_var(const unsigned long i) const
{
  assert(i < num_vars + 1);
  
  unsigned long j;
  long*         dp;
//  K_RATPOLY     P;
  
  dp = new long [num_vars + 1];
  
  for (j = 0; j < i; j++)
    dp[j] = deg[j];
  
  dp[i] = 0;
  
  for (j = i; j < num_vars; j++)
    dp[j + 1] = deg[j];
  
//  P = K_RATPOLY(num_vars + 1, dp);
  K_RATPOLY P(num_vars + 1, dp);
  
  if (dp)
    delete [] dp;
  
  for (j = 0; j < num_coeffs; j++)
    P.coeffs[j] = coeffs[j];
  
  return P;
}

//  K_RATPOLY K_RATPOLY :: remove_var(const unsigned long i) const
//    PROVIDED deg[i] == 0,
//    returns a polynomial of num_vars - 1 variables of degree
//      deg[0], deg[1], ..., deg[i - 1], deg[i + 1], ..., deg[num_vars - 1].

K_RATPOLY K_RATPOLY :: remove_var(const unsigned long i) const
{
  assert(i < num_vars);  //  i < num_vars => num_vars > 0
  assert(deg[i] == 0);
  
  unsigned long  j;
  long*          dp;
//  K_RATPOLY      P;
  
  if (num_vars > 1)
  {
    dp = new long [num_vars - 1];
    
    for (j = 0; j < i; j++)
      dp[j] = deg[j];
    
    for (j = i + 1; j < num_vars; j++)
      dp[j - 1] = deg[j];
  }
  else  //  if (num_vars == 1)
    dp = 0;
  
//  P = K_RATPOLY(num_vars - 1, dp);
  K_RATPOLY P(num_vars - 1, dp);
  
  if (dp)
    delete [] dp;
  
  for (j = 0; j < num_coeffs; j++)
    P.coeffs[j] = coeffs[j];
  
  return P;
}

//  K_RATPOLY K_RATPOLY :: add(const K_RATPOLY& P) const
//    returns *this + P.

K_RATPOLY K_RATPOLY :: add(const K_RATPOLY& P) const
{
  assert(num_vars == P.num_vars);
  
  unsigned long i, j;
  long*         dq;
  long*         pq;
//  K_RATPOLY     Q;
  
  if (num_vars > 0)
  {
    dq = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dq[i] = deg[i] > P.deg[i] ? deg[i] : P.deg[i];
  }
  else  //  if (num_vars == 0)
    dq = 0;
  
//  Q = K_RATPOLY(num_vars, dq);
  K_RATPOLY Q(num_vars, dq);
  
  if (dq)
    delete [] dq;
  
  for (i = 0; i < num_coeffs; i++)
  {
    pq          = index_to_powers(i);
    j           = Q.powers_to_index(pq);
    Q.coeffs[j] = coeffs[i];
    
    if (pq)
      delete [] pq;
  }
  
  for (i = 0; i < P.num_coeffs; i++)
  {
    pq           = P.index_to_powers(i);
    j            = Q.powers_to_index(pq);
    Q.coeffs[j] += P.coeffs[i];
    
    if (pq)
      delete [] pq;
  }
  
  Q.reduce_deg();
  
  return Q;
}

K_RATPOLY operator +(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  return P1.add(P2);
}

//  K_RATPOLY K_RATPOLY :: sub(const K_RATPOLY& P) const
//    returns *this - P.

K_RATPOLY K_RATPOLY :: sub(const K_RATPOLY& P) const
{
  assert(num_vars == P.num_vars);
  
  unsigned long i, j;
  long*         dq;
  long*         pq;
//  K_RATPOLY     Q;
  
  if (num_vars > 0)
  {
    dq = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dq[i] = deg[i] > P.deg[i] ? deg[i] : P.deg[i];
  }
  else  //  if (num_vars == 0)
    dq = 0;
  
//  Q = K_RATPOLY(num_vars, dq);
  K_RATPOLY Q(num_vars, dq);
  
  if (dq)
    delete [] dq;
  
  for (i = 0; i < num_coeffs; i++)
  {
    pq          = index_to_powers(i);
    j           = Q.powers_to_index(pq);
    Q.coeffs[j] = coeffs[i];
    
    if (pq)
      delete [] pq;
  }
  
  for (i = 0; i < P.num_coeffs; i++)
  {
    pq           = P.index_to_powers(i);
    j            = Q.powers_to_index(pq);
    Q.coeffs[j] -= P.coeffs[i];
    
    if (pq)
      delete [] pq;
  }
  
  Q.reduce_deg();
  
  return Q;
}

K_RATPOLY operator -(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  return P1.sub(P2);
}

//  K_RATPOLY K_RATPOLY :: mul(const K_RATPOLY& P) const
//    returns *this * P.

K_RATPOLY K_RATPOLY :: mul(const K_RATPOLY& P) const
{
  assert(num_vars == P.num_vars);
  
  unsigned long i, j, k;
  long*         dq;
  long*         p;
  long*         pp;
  long*         pq;
//  K_RATPOLY     Q;
  
  if (num_vars > 0)
  {
    dq = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dq[i] = deg[i] + P.deg[i];
  }
  else  //  if (num_vars == 0)
    dq = 0;
  
//  Q = K_RATPOLY(num_vars, dq);
  K_RATPOLY Q(num_vars, dq);
  
  for (i = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]))
    {
      p = index_to_powers(i);
      
      for (j = 0; j < P.num_coeffs; j++)
        if (sgn(P.coeffs[j]))
        {
          pp = P.index_to_powers(j);
          
          if (num_vars > 0)
          {
            pq = new long [num_vars];
            
            for (k = 0; k < num_vars; k++)
              pq[k] = p[k] + pp[k];
          }
          else  //  if (num_vars == 0)
            pq = 0;
          
          k            = Q.powers_to_index(pq);
          Q.coeffs[k] += coeffs[i] * P.coeffs[j];
          
          if (pq)
            delete [] pq;
          
          if (pp)
            delete [] pp;
        }
      
      if (p)
        delete [] p;
    }
  
  Q.reduce_deg();
  
  return Q;
}

K_RATPOLY operator *(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  return P1.mul(P2);
}

//  K_RATPOLY K_RATPOLY :: mul(const bigrational& b) const
//    returns *this * b.

K_RATPOLY K_RATPOLY :: mul(const bigrational& b) const
{
  unsigned long i;
  long*         dp;
//  K_RATPOLY     P;
  
  if (num_vars > 0)
  {
    dp = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dp[i] = sgn(b) ? deg[i] : 0;
  }
  else  //  if (num_vars == 0)
    dp = 0;
  
//  P = K_RATPOLY(num_vars, dp);
  K_RATPOLY P(num_vars, dp);
  
  if (dp)
    delete [] dp;
  
  if (sgn(b))
    for (i = 0; i < num_coeffs; i++)
      P.coeffs[i] = coeffs[i] * b;
//  Assume the constructor for P initializes all its coefficients to be 0.
//  else  //  if (!sgn(b))
//    for (i = 0; i < P.num_coeffs; i++)
//      P.coeffs[i] = 0;
  
  return P;
}

K_RATPOLY operator *(const K_RATPOLY& P, const bigrational& b)
{
  return P.mul(b);
}

K_RATPOLY operator *(const bigrational& b, const K_RATPOLY& P)
{
  return P.mul(b);
}

//  K_RATPOLY K_RATPOLY :: neg()
//    returns - *this.

K_RATPOLY K_RATPOLY :: neg() const
{
  unsigned long i;
//  K_RATPOLY     P;
  
//  P = K_RATPOLY(num_vars, deg);
  K_RATPOLY P(num_vars, deg);
  
  for (i = 0; i < num_coeffs; i++)
    P.coeffs[i] = - coeffs[i];
  
  return P;
}

K_RATPOLY operator -(const K_RATPOLY& P)
{
  return P.neg();
}

//  K_RATPOLY K_RATPOLY :: exact_div(const K_RATPOLY& P) const
//    returns *this / P PROVIDED *this is divisible by P.

K_RATPOLY K_RATPOLY :: exact_div(const K_RATPOLY& P) const
{
  assert(num_vars == P.num_vars);
  
  unsigned long i, j, k, l, m;
  long*         dq;
  long*         pp0;
  long*         pt0;
  long*         pq;
  long*         pp1;
  long*         pt1;
  int           zero_divisor, zero_dividend;
//  K_RATPOLY     T;
//  K_RATPOLY     Q;
  
  if (num_vars > 0)
    dq = new long [num_vars];
  else  //  if (num_vars == 0)
    dq = 0;
  
  //  Make sure that
  //    in each variable, dividend is of a higher degree than divisor.
  
  for (i = 0; i < num_vars; i++)
    assert((dq[i] = deg[i] - P.deg[i]) >= 0);
  
//  Q = K_RATPOLY(num_vars, dq);
  K_RATPOLY Q(num_vars, dq);
  
  if (dq)
    delete [] dq;
  
  //  Find the leading term of divisor.
  
  zero_divisor = 1;
  i            = 0;
  
  while (zero_divisor && i < P.num_coeffs)
    if (sgn(P.coeffs[i]))
      zero_divisor = 0;
    else  //  if (!sgn(P.coeffs[i]))
      i++;
  
  //  Make sure divisor is non-zero.
  
  assert(!zero_divisor);  //  !zero_divisor => i < P.num_coeffs
  
  pp0 = P.index_to_powers(i);
  
  //  Find the leading term of dividend.
  
//  T             = *this;
  K_RATPOLY T(*this);
  zero_dividend = 1;
  j             = 0;
  
  while (zero_dividend && j < T.num_coeffs)
    if (sgn(T.coeffs[j]))
      zero_dividend = 0;
    else  //  if (!sgn(T.coeffs[j]))
      j++;
  
  while (!zero_dividend && j < T.num_coeffs)
  {
    pt0 = T.index_to_powers(j);
    
    if (num_vars > 0)
    {
      pq = new long [num_vars];
      
      for (k = 0; k < num_vars; k++)
        pq[k] = pt0[k] - pp0[k];
    }
    else  //  if (num_vars == 0)
      pq = 0;
    
    //  Compute another non-vanishing coefficient of quotient.
    
    k           = Q.powers_to_index(pq);
    Q.coeffs[k] = T.coeffs[j] / P.coeffs[i];
    
    //  dividend -= divisor * the non-vanishing coefficient of quotient.
    
    for (l = i; l < P.num_coeffs; l++)
      if (sgn(P.coeffs[l]))
      {
        pp1 = P.index_to_powers(l);
        
        if (num_vars > 0)
        {
          pt1 = new long [num_vars];
          
          for (m = 0; m < num_vars; m++)
            pt1[m] = pp1[m] + pq[m];
        }
        else  //  if (num_vars == 0)
          pt1 = 0;
        
        m            = T.powers_to_index(pt1);
        T.coeffs[m] -= P.coeffs[l] * Q.coeffs[k];
        
        if (pt1)
          delete [] pt1;
        
        if (pp1)
          delete [] pp1;
      }
    
    if (pq)
      delete [] pq;
    
    if (pt0)
      delete [] pt0;
    
    //  Find the leading term of "updated" dividend.
    
    zero_dividend = 1;
    
    while (zero_dividend && j < T.num_coeffs)
      if (sgn(T.coeffs[j]))
        zero_dividend = 0;
      else  //  if (!sgn(T.coeffs[j]))
        j++;
  }
  
  //  Make sure division is exact.
  
  assert(zero_dividend);
  
  if (pp0)
    delete [] pp0;
  
  Q.reduce_deg();
  
  return Q;
}

//  K_RATPOLY exact_div(const K_RATPOLY& P1, const K_RATPOLY& P2)
//    returns P1 / P2 PROVIDED P1 is divisible by P2.

K_RATPOLY exact_div(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  return P1.exact_div(P2);
}

//  K_RATPOLY K_RATPOLY :: derivative(const unsigned long i) const
//     returns the partial derivative of *this w.r.t. X_i.

K_RATPOLY K_RATPOLY :: derivative(const unsigned long i) const
{
  assert(num_vars == 0 || i < num_vars);
  
  unsigned long j, k;
  unsigned long c;
  long*         dp;
  long*         p;
//  K_RATPOLY     P;
  
  if (num_vars > 0)
  {
    dp = new long [num_vars];
    
    if (deg[i] > 0)
    {
      for (j = 0; j < i; j++)
        dp[j] = deg[j];
      
      dp[i] = deg[i] - 1;
      
      for (j = i + 1; j < num_vars; j++)
        dp[j] = deg[j];
    }
    else  //  if (deg[i] == 0)
      for (j = 0; j < num_vars; j++)
        dp[j] = 0;
  }
  else  //  if (num_vars == 0)
    dp = 0;
  
//  P = K_RATPOLY(num_vars, dp);
  K_RATPOLY P(num_vars, dp);
  
  if (dp)
    delete [] dp;
  
  for (j = 0; j < num_coeffs; j++)
    if (sgn(coeffs[j]) && (p = index_to_powers(j)))
    {
      if ((c = p[i]) > 0)
      {
        p[i]--;
        k           = P.powers_to_index(p);
        P.coeffs[k] = coeffs[j] * c;
      }
      
      delete [] p;
    }
  
  P.reduce_deg();
  
  return P;
}

//  int K_RATPOLY :: cmp(const K_RATPOLY& P) const
//    returns 0     if *this and P are the same. i.e.,
//                       deg and coeffs of *this and P are the same, and
//    returns non-0 otherwise.

int K_RATPOLY :: cmp(const K_RATPOLY& P) const
{
  unsigned long i;
  int           c;
  
  c = 0;
  
  if (num_vars != P.num_vars)
    c = 1;
  else  //  if (num_vars == P.num_vars)
  {
    i = 0;
    
    while (!c && i < num_vars)
      if (deg[i] != P.deg[i])
        c = 1;
      else  // if (deg[i] == P.deg[i])
        i++;
    
    if (!c)
    {
      if (num_coeffs != P.num_coeffs)
        c = 1;
      else  //  if (num_coeffs == P.num_coeffs)
      {
        i = 0;
        
        while (!c && i < num_coeffs)
          if (coeffs[i] != P.coeffs[i])
            c = 1;
          else  //  if (coeffs[i] == P.coeffs[i])
            i++;
      }
      
//      if (!c && Sturm_seq != P.Sturm_seq)
//        c = 1;
    }
  }
  
  return c;
}

int operator ==(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  return P1.cmp(P2) == 0;
}

//  int K_RATPOLY :: is_zero() const
//    returns 1 if *this is a zero polynomial, and
//            0 otherwise.

int K_RATPOLY :: is_zero() const
{
  unsigned long i;
  int           z;
  
  for (i = 0; i < num_coeffs && !sgn(coeffs[i]); i++)
    ;
  
  if (i < num_coeffs)
    z = 0;
  else  //  if (i == num_coeffs)
    z = 1;
  
  return z;
}

//  int K_RATPOLY :: eq_upto_const(const K_RATPOLY& P) const
//    returns 1 if *this == P * b for some bigrational b, and
//            0  otherwise.

int K_RATPOLY :: eq_upto_const(const K_RATPOLY& P) const
{
  long i;
  int  e;
  
  e = 1;
  
  if (num_vars != P.num_vars)
    e = 0;
  else  //  if (num_vars == P.num_vars)
  {
    i = 0;
    
    while (e && i < num_vars)
      if (deg[i] != P.deg[i])
        e = 0;
      else  //  if (deg[i] == P.deg[i])
        i++;
    
    if (e)
      if (num_coeffs != P.num_coeffs)
        e = 0;
      else  //  if (num_coeffs == P.num_coeffs)
      {
        bigrational lc, lcp;
        
        for (i = 0; i < num_coeffs && !sgn(coeffs[i]); i++)
          ;
        
        if (i < num_coeffs)
        {
          //  *this is NOT a zero polynomial.
          
          lc = coeffs[i];
          
          if (!sgn(lcp = P.coeffs[i]))
          //  *this is NOT a zero polynomial and P is.
            e = 0;
          else
          {
            i++;
            
            while (e && i < num_coeffs)
              if (coeffs[i] * lcp != lc * P.coeffs[i])
                e = 0;
              else  //  if (coeffs[i] * lcp == lc * P.coeffs[i])
                i++;
          }
        }
        else  //  if (i == num_coeffs)
        {
          //  *this is a zero polynomial.
          
          i--;
          
          while (e && i >= 0)
            if (sgn(P.coeffs[i]))
            //  * this is a zero polynomial and P is not.
              e = 0;
            else  //  if (!sgn(P.coeffs[i]))
              i--;
        }
      }
  }
  
  return e;
}

//  K_RATPOLY K_RATPOLY :: subst_val_first_var_proto(const bigrational& b,
//                                                   const int reduced) const
//    returns a polynomial of num_vars - 1 variables obtained by
//    computing *this(b, X_1, ..., X_{num_vars - 1}), and
//    renaming X_1, ..., X_{num_vars - 1} to X_0, X_1, ..., X_{num_vars - 2}.
//    reduce_deg() is applied if reduced == 1.

K_RATPOLY K_RATPOLY :: subst_val_first_var_proto(const bigrational& b,
                                                 const int reduced) const
{
  unsigned long i, j;
  long*         dp;
  K_RATPOLY     P;
  
  if (num_vars > 0)
  {
    if (num_vars > 1)
    {
      dp = new long [num_vars - 1];
      
      for (i = 1; i < num_vars; i++)
        dp[i - 1] = deg[i];
    }
    else  //  if (num_vars == 1)
      dp = 0;
    
    P = K_RATPOLY(num_vars - 1, dp);
    
    if (dp)
      delete [] dp;
    
    for (i = 0; i < P.num_coeffs; i++)
      P.coeffs[i] = coeffs[i];
    
    for (i = 1; i <= deg[0]; i++)
      for (j = 0; j < P.num_coeffs; j++)
      {
        P.coeffs[j] *= b;
        P.coeffs[j] += coeffs[j + i * P.num_coeffs];
      }
    
    if (reduced)
      P.reduce_deg();
  }
  else  //  if (num_vars == 0)
    P = *this;
  
  return P;
}

//  K_RATPOLY K_RATPOLY :: subst_val_first_var(const bigrational& b) const
//    returns a polynomial of num_vars - 1 variables obtained by
//    computing *this(b, X_1, ..., X_{num_vars - 1}), and
//    renaming X_1, ..., X_{num_vars - 1} to X_0, X_1, ..., X_{num_vars - 2}.

K_RATPOLY K_RATPOLY :: subst_val_first_var(const bigrational& b) const
{
  return subst_val_first_var_proto(b, 1);
}

//  bigrational K_RATPOLY :: evaluate(const bigrational& b) const
//    returns *this(b) PROVIDED *this is a univariate polynomial.

bigrational K_RATPOLY :: evaluate(const bigrational& b) const
{
  assert(num_vars == 1);
  
  return subst_val_first_var(b).coeffs[0];
}

//  bigrational K_RATPOLY :: evaluate(const bigrational* const B) const
//    returns *this(B[0], B[1], ..., B[num_vars - 1]).

bigrational K_RATPOLY :: evaluate(const bigrational* const B) const
{
  unsigned long i;
//  K_RATPOLY     P;
  
//  P = *this;
  K_RATPOLY P(*this);
  
  for (i = 0; i < num_vars; i++)
    P = P.subst_val_first_var(B[i]);
  
  return P.coeffs[0];
}

//  bigrational K_RATPOLY :: evaluate(const bigrational_vector& B) const
//    returns *this(B[0], B[1], ..., B[num_vars - 1]).

bigrational K_RATPOLY :: evaluate(const bigrational_vector& B) const
{
  assert(num_vars == B.get_dim());
  
  unsigned long i;
//  K_RATPOLY     P;
  
//  P = *this;
  K_RATPOLY P(*this);
  
  for (i = 0; i < num_vars; i++)
    P = P.subst_val_first_var(B[i]);
  
  return P.coeffs[0];
}

//  int K_RATPOLY :: fp_sgn_at(const bigrational& b) const
//    PROVIDED *this is a univariate polynomial,
//    returns   1 if *this(b) is certainly positive,
//            - 1 if *this(b) is certainly negative, and
//              0 otherwise.

int K_RATPOLY :: fp_sgn_at(const bigrational& b) const
{
  assert(num_vars == 1);
  
  unsigned long i;
  double*       c;
  double        val, err;
  int           s;
  
  c = new double [num_coeffs];
  
  for (i = 0; i < num_coeffs && finite(c[i] = coeffs[i].as_double()); i++)
    ;
  
  if (i < num_vars)
  //  if some c[i] is not a finite floating number
    s = 0;
  else  //  if (i == num_vars)
  //  if all c[i]'s are finite floating numbers
    if (horner_with_err(c, deg[0], RAT_TO_DBL_ERR,
                        b.as_double(), RAT_TO_DBL_ERR,
                        val, err))
      s = sign_given_err(val, err);
    else  //  if an overflow occurs at horner_with_err
      s = 0;
  
  delete [] c;
  
  return s;
}

//  int K_RATPOLY :: sgn_at(const bigrational& b) const
//    PROVIDED *this is a univariate polynomial,
//    returns the sign of *this(b).

int K_RATPOLY :: sgn_at(const bigrational& b) const
{
  assert(num_vars == 1);
  
#ifdef _EXPERIMENT
  num_kratpoly_sgn_at++;
#endif
  
  int fs, es, s;
  
  if (fs = fp_sgn_at(b))
    s = fs;
  else  //  if (!fs)
  {
#ifdef _EXPERIMENT
    num_kratpoly_exact_sgn_at++;
#endif
    
    s = sgn(evaluate(b));
  }
  
  return s;
}

//  int K_RATPOLY :: sgn_at(const bigrational* const B)
//    return the sign of *this(B[0], B[1], ..., B[num_vars - 1]).

int K_RATPOLY :: sgn_at(const bigrational* const B) const
{
  if (num_vars == 1)
    return sgn_at(B[0]);
  else  //  if (num_vars > 1)
    return sgn(evaluate(B));
}

//  int K_RATPOLY :: sgn_at(const bigrational_vector& B)
//    return the sign of *this(B[0], B[1], ..., B[num_vars - 1]).

int K_RATPOLY :: sgn_at(const bigrational_vector& B) const
{
  if (num_vars == 1)
    return sgn_at(B[0]);
  else  //  if (num_vars > 1)
    return sgn(evaluate(B));
}

//  K_RATPOLY K_RATPOLY :: subst_val_proto(const unsigned long i,
//                                         const bigrational& b,
//                                         const int reduced) const
//    returns a polynomial of num_vars - 1 variables obtained by
//    computing
//      *this(X_0, X_1, ..., X_{i - 1}, b, X_{i + 1}, ..., X_{num_vars - 1})
//    and, renaming
//           X_{i + 1}, ..., X_{num_vars - 1} to X_i, ..., X_{num_vars - 2}.
//    reduce_deg() is applied if reduced == 1.

K_RATPOLY K_RATPOLY :: subst_val_proto(const unsigned long i,
                                       const bigrational& b,
                                       const int reduced) const
{
  assert(num_vars == 0 || i < num_vars);
  
  unsigned long j, k;
  long*         dp;
  bigrational*  bp;
  long*         p;
  long*         pp;
  K_RATPOLY     P;
  
  if (num_vars > 0)
  {
    if (i == 0)
      P = subst_val_first_var_proto(b, reduced);
    else  //  if (i > 0)
    {
      if (num_vars > 1)
      {
        dp = new long [num_vars - 1];
        
        for (j = 0; j < i; j++)
          dp[j] = deg[j];
        
        for (j = i + 1; j < num_vars; j++)
          dp[j - 1] = deg[j];
      }
      else  //  if (num_vars == 1)
        dp = 0;
      
      P = K_RATPOLY(num_vars - 1, dp);
      
      if (dp)
        delete [] dp;
      
      bp    = new bigrational [deg[i] + 1];
      bp[0] = 1;
      
      if (deg[i] > 0)
      {
        bp[1] = b;
        
        for (j = 1; j <= deg[i]; j++)
          bp[j] = b * bp[j - 1];
      }
      
      for (j = 0; j < num_coeffs; j++)
        if (sgn(coeffs[j]))
        {
          p = index_to_powers(j);
          
          if (num_vars > 1)
          {
            pp = new long [num_vars - 1];
            
            for (k = 0; k < i; k++)
              pp[k] = p[k];
            
            for (k = i + 1; k < num_vars; k++)
              pp[k - 1] = p[k];
          }
          else  //  if (num_vars == 1)
            pp = 0;
          
          k            = P.powers_to_index(pp);
          P.coeffs[k] += coeffs[j] * bp[p[i]];
          
          if (pp)
            delete [] pp;
          
          delete [] p;  //  num_vars > 0 => p != 0
        }
      
      delete [] bp;
      
      if (reduced)
        P.reduce_deg();
    }
  }
  else  //  if (num_vars == 0)
    P = *this;
  
  return P;
}

//  K_RATPOLY K_RATPOLY :: subst_val(const unsigned long i,
//                                   const bigrational& b) const
//    returns a polynomial of num_vars - 1 variables obtained by
//    computing
//      *this(X_0, X_1, ..., X_{i - 1}, b, X_{i + 1}, ..., X_{num_vars - 1})
//    and, renaming
//           X_{i + 1}, ..., X_{num_vars - 1} to X_i, ..., X_{num_vars - 2}.

K_RATPOLY K_RATPOLY :: subst_val(const unsigned long i,
                                 const bigrational& b) const
{
  return subst_val_proto(i, b, 1);
}

//  K_RATPOLY K_RATPOLY :: subst_expr(const unsigned long i,
//                                    const K_RATPOLY& P) const
//    returns a polynomial of num_vars - 1 variables obtained by
//    computing
//      *this(X_0, X_1, ..., X_{i - 1},
//            P (X_0, X_1, ..., X_{i - 1}, X_{i + 1}, X_{num_vars - 1}),
//            X_{i + 1}, ..., X_{num_vars - 1})
//    and, renaming
//           X_{i + 1}, ..., X_{num_vars - 1} to X_i, ..., X_{num_vars - 2}.

K_RATPOLY K_RATPOLY :: subst_expr(const unsigned long i,
                                  const K_RATPOLY& P) const
{
  assert(i < num_vars);                //  i < num_vars => num_vars > 0
  assert(num_vars == P.num_vars + 1);
  
  unsigned long j, k, l;
  long*         dq;
  K_RATPOLY*    PP;
  long*         dpp0;
  long*         p;
  long*         p1;
  long*         ppp;
  long*         pq;
//  K_RATPOLY     Q;
  
  if (num_vars > 1)
  {
    dq = new long [num_vars - 1];
    
    for (j = 0; j < i; j++)
      dq[j] = deg[j] + deg[i] * P.deg[j];
    
    for (j = i + 1; j < num_vars; j++)
      dq[j - 1] = deg[j] + deg[i] * P.deg[j - 1];
  }
  else  //  if (num_vars == 1)
    dq = 0;
  
//  Q = K_RATPOLY(num_vars - 1, dq);
  K_RATPOLY Q(num_vars - 1, dq);
  
  if (dq)
    delete [] dq;
  
  PP = new K_RATPOLY [deg[i] + 1];
  
  if (P.num_vars > 0)
  {
    dpp0 = new long [P.num_vars];
    
    for (j = 0; j < P.num_vars; j++)
      dpp0[j] = 0;
  }
  else  //  if (P.num_vars == 0)
    dpp0 = 0;
  
  PP[0] = K_RATPOLY(P.num_vars, dpp0);
  
  if (dpp0)
    delete [] dpp0;
  
  PP[0].coeffs[0] = 1;
  
  if (deg[i] > 0)
  {
    PP[1] = P;
    
    for (j = 2; j <= deg[i]; j++)
      PP[j] = P * PP[j - 1];
  }
  
  for (j = 0; j < num_coeffs; j++)
    if (sgn(coeffs[j]))
    {
      p = index_to_powers(j);  //  num_vars > 0 => p != 0
      
      if (num_vars > 1)
      {
        p1 = new long [num_vars - 1];
        
        for (k = 0; k < i; k++)
          p1[k] = p[k];
        
        for (k = i + 1; k < num_vars; k++)
          p1[k - 1] = p[k];
      }
      else  //  if (num_vars == 1)
        p1 = 0;
      
      for (k = 0; k < PP[p[i]].num_coeffs; k++)
        if (sgn(PP[p[i]].coeffs[k]))
        {
          ppp = PP[p[i]].index_to_powers(k);
          
          if (num_vars > 1)
          {
            pq = new long [num_vars - 1];
            
            for (l = 0; l < num_vars - 1; l++)
              pq[l] = p1[l] + ppp[l];
          }
          else  //  if (num_vars == 1)
            pq = 0;
          
          l            = Q.powers_to_index(pq);
          Q.coeffs[l] += coeffs[j] * PP[p[i]].coeffs[k];
          
          if (pq)
            delete [] pq;
          
          if (ppp)
            delete [] ppp;
        }
      
      if (p1)
        delete [] p1;
      
      delete [] p;  //  num_vars > 0 => p != 0
    }
  
  delete [] PP;
  Q.reduce_deg();
  
  return Q;
}

//  K_RATPOLY K_RATPOLY :: subst_expr(const unsigned long i,
//                                    const K_RATPOLY& N,
//                                    const K_RATPOLY& D) const
//    returns a polynomial of num_vars - 1 variables obtained by
//    computing
//      *this(X_0, X_1, ..., X_{i - 1},
//            N/D (X_0, X_1, ..., X_{i - 1}, X_{i + 1}, X_{num_vars - 1}),
//            X_{i + 1}, ..., X_{num_vars - 1}),
//    and, renaming
//           X_{i + 1}, ..., X_{num_vars - 1} to X_i, ..., X_{num_vars - 2}.

K_RATPOLY K_RATPOLY :: subst_expr(const unsigned long i,
                                  const K_RATPOLY& N,
                                  const K_RATPOLY& D) const
{
  assert(i < num_vars);                //  i < num_vars => num_vars > 0
  assert(num_vars == N.num_vars + 1);
  assert(num_vars == D.num_vars + 1);
  
  unsigned long j, k, l, m;
  long*         dq;
  K_RATPOLY*    NP;
  K_RATPOLY*    DP;
  long*         dnp0;
  long*         ddp0;
  long*         p;
  long*         p1;
  long*         pnp;
  long*         pdp;
  long*         pq;
//  K_RATPOLY     Q;
  
  if (num_vars > 1)
  {
    dq = new long [num_vars - 1];
    
    for (j = 0; j < i; j++)
      if (N.deg[j] > D.deg[j])
        dq[j] = deg[j] + deg[i] * N.deg[j];
      else  //  if (N.deg[j] <= D.deg[j])
        dq[j] = deg[j] + deg[i] * D.deg[j];
    
    for (j = i + 1; j < num_vars; j++)
      if (N.deg[j - 1] > D.deg[j - 1])
        dq[j - 1] = deg[j] + deg[i] * N.deg[j - 1];
      else  //  if (N.deg[j - 1] <= D.deg[j - 1])
        dq[j - 1] = deg[j] + deg[i] * D.deg[j - 1];
  }
  else  //  if (num_vars == 1)
    dq = 0;
  
//  Q = K_RATPOLY(num_vars - 1, dq);
  K_RATPOLY Q(num_vars - 1, dq);
  
  if (dq)
    delete [] dq;
  
  NP = new K_RATPOLY [deg[i] + 1];
  
  if (N.num_vars > 0)
  {
    dnp0 = new long [N.num_vars];
    
    for (j = 0; j < N.num_vars; j++)
      dnp0[j] = 0;
  }
  else  //  if (N.num_vars == 0)
    dnp0 = 0;
  
  NP[0] = K_RATPOLY(N.num_vars, dnp0);
  
  if (dnp0)
    delete [] dnp0;
  
  NP[0].coeffs[0] = 1;
  
  if (deg[i] > 0)
  {
    NP[1] = N;
    
    for (j = 2; j <= deg[i]; j++)
      NP[j] = N * NP[j - 1];
  }
  
  DP = new K_RATPOLY [deg[i] + 1];
  
  if (D.num_vars > 0)
  {
    ddp0 = new long [D.num_vars];
    
    for (j = 0; j < D.num_vars; j++)
      ddp0[j] = 0;
  }
  else  //  if (D.num_vars == 0)
    ddp0 = 0;
  
  DP[0] = K_RATPOLY(D.num_vars, ddp0);
  
  if (ddp0)
    delete [] ddp0;
  
  DP[0].coeffs[0] = 1;
  
  if (deg[i] > 0)
  {
    DP[1] = D;
    
    for (j = 2; j <= deg[i]; j++)
      DP[j] = D * DP[j - 1];
  }
  
  for (j = 0; j < num_coeffs; j++)
    if (sgn(coeffs[j]))
    {
      p = index_to_powers(j);  //  i < num_vars => p != 0
      
      if (num_vars > 1)
      {
        p1 = new long [num_vars - 1];
        
        for (k = 0; k < i; k++)
          p1[k] = p[k];
        
        for (k = i + 1; k < num_vars; k++)
          p1[k - 1] = p[k];
      }
      else
        p1 = 0;
      
      for (k = 0; k < NP[p[i]].num_coeffs; k++)
        if (sgn(NP[p[i]].coeffs[k]))
        {
          pnp = NP[p[i]].index_to_powers(k);
          
          for (l = 0; l < DP[deg[i] - p[i]].num_coeffs; l++)
            if (sgn(DP[deg[i] - p[i]].coeffs[l]))
            {
              pdp = DP[deg[i] - p[i]].index_to_powers(l);
              
              if (num_vars > 1)
              {
                pq = new long [num_vars - 1];
                
                for (m = 0; m < num_vars - 1; m++)
                  pq[m] = p1[m] + pnp[m] + pdp[m];
              }
              else  //  if (num_vars == 1)
                pq = 0;
              
              m            = Q.powers_to_index(pq);
              Q.coeffs[m] += coeffs[j] * NP[p[i]].coeffs[k]
                                       * DP[deg[i] - p[i]].coeffs[l];
              
              if (pq)
                delete [] pq;
              
              if (pdp)
                delete [] pdp;
            }
          
          if (pnp)
            delete [] pnp;
        }
      
      if (p1)
        delete [] p1;
      
      delete [] p;  //  num_vars > 0 => p != 0
    }
  
  delete [] NP;
  delete [] DP;
  Q.reduce_deg();
  
  return Q;
}

//  K_RATPOLY K_RATPOLY :: subst_param_expr(const K_RATPOLY& X,
//                                          const K_RATPOLY& Y,
//                                          const K_RATPOLY& Z,
//                                          const K_RATPOLY& W) const
//    PROVIDED *this is a 3-variate polynomial and
//             X, Y, Z, W are bivariate polynomials,
//    returns *this(X/W, Y/W, Z/W).

K_RATPOLY K_RATPOLY :: subst_param_expr(const K_RATPOLY& X,
                                        const K_RATPOLY& Y,
                                        const K_RATPOLY& Z,
                                        const K_RATPOLY& W) const
{
  assert(num_vars == 3);
  assert(X.num_vars == 2);
  assert(Y.num_vars == 2);
  assert(Z.num_vars == 2);
  assert(W.num_vars == 2);
  
  unsigned long i, j;
  long          d[2];
  long          td, t, dd;
  long*         p;
  K_RATPOLY     S;
//  K_RATPOLY     I;
  
  d[0] = 0;
  d[1] = 0;
  
//  I = K_RATPOLY(2, d);
  K_RATPOLY I(2, d);
  
  for (i = 0, td = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]))
    {
      p = index_to_powers(i);
      
      if ((t = p[0] + p[1] + p[2]) > td)
        td = t;
      
      delete [] p;  //  num_vars == 3 => p != 0
    }
  
  for (i = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]))
    {
      S = K_RATPOLY(2, d);
      
      S.coeffs[0] = coeffs[i];
      p           = index_to_powers(i);
      dd          = td - p[0] - p[1] - p[2];
      
      for (j = 0; j < p[0]; j++)
        S = S * X;
      
      for (j = 0; j < p[1]; j++)
        S = S * Y;
      
      for (j = 0; j < p[2]; j++)
        S = S * Z;
      
      for (j = 0; j < dd; j++)
        S = S * W;
      
      I = I + S;
      
      delete [] p;  //  num_vars == 3 => p != 0
    }
  
  I.reduce_deg();
  I.reduce_num_coeffs();
  
  return I;
}

//  K_FLOATPOLY K_RATPOLY :: as_FLOATPOLY() const
//    returns a polynomial whose coefficients are
//                           floating point approximations of those of *this.

K_FLOATPOLY K_RATPOLY :: as_FLOATPOLY() const
{
  unsigned long i;
  bigrational   b;
  double*       c;
//  K_FLOATPOLY   F;
  
  b = 0;
  
  for (i = 0; i < num_coeffs; i++)
    if (abs(coeffs[i]) > b)
      b = abs(coeffs[i]);
  
  c = new double [num_coeffs];
  
  for (i = 0; i < num_coeffs; i++)
    if (sgn(b))
      c[i] = (coeffs[i] / b).as_double();
    else  //  if (!sgn(b))
      c[i] = 0.0;
  
//  F = K_FLOATPOLY(num_vars, deg, c);
  K_FLOATPOLY F(num_vars, deg, c);
  
  delete [] c;
  
  return F;
}

//  K_RATPOLY K_RATPOLY :: div(const K_RATPLOY& P, K_RATPOLY& R) const
//    PROVIDED *this and P are univariate polynomials,
//    computes univariate polynomials Q and R s.t.
//      *this = P * Q + R with deg R < deg P
//    returns Q.

K_RATPOLY K_RATPOLY :: div(const K_RATPOLY& P, K_RATPOLY& R) const
{
  assert(num_vars == 1);
  assert(P.num_vars == 1);
  assert(sgn(P.coeffs[0]));
  
  unsigned long i, j;
  long          dq[1];
  bigrational   q;
  K_RATPOLY     Q;
  
  if (deg[0] >= P.deg[0])
  {
    dq[0] = deg[0] - P.deg[0];
    Q     = K_RATPOLY(1, dq);
    
    R = *this;
    
    for (i = 0; i <= Q.deg[0]; i++)
    {
      Q.coeffs[i] = q = R.coeffs[i] / P.coeffs[0];
      R.coeffs[i] = 0;
      
      for (j = 1; j <= P.deg[0]; j++)
        R.coeffs[i + j] -= P.coeffs[j] * q;
    }
    
    Q.reduce_deg();
    R.reduce_deg();
  }
  else  //  if (deg[0] < P.deg[0])
  {
    dq[0]       = 0;
    Q           = K_RATPOLY(1, dq);
    Q.coeffs[0] = 0;
    
    R = *this;
  }
  
  return Q;
}

//  K_RATPOLY div(const K_RATPOLY& P1, const K_RATPLOY& P2, K_RATPOLY& R) const
//    PROVIDED P1 and P2 are univariate polynomials,
//    computes univariate polynomials Q and R s.t.
//      P1 = P2 * Q + R with deg R < deg P2
//    returns Q.

K_RATPOLY div(const K_RATPOLY& P1, const K_RATPOLY& P2, K_RATPOLY& R)
{
  return P1.div(P2, R);
}

K_RATPOLY operator /(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  K_RATPOLY R;
  
  return P1.div(P2, R);
}

K_RATPOLY rem(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  K_RATPOLY R;
  
  P1.div(P2, R);
  
  return R;
}

//  K_RATPOLY* K_RATPOLY :: set_Sturm_seq()
//    PROVIDED *this is a univariate polynomial,
//    generates Sturm sequence of *this.

K_RATPOLY* K_RATPOLY :: set_Sturm_seq()
{
  assert(num_vars == 1);
  
  if (deg[0] && !Sturm_seq)
  {
    K_RATPOLY* P;
    K_RATPOLY* Q;
    K_RATPOLY* R;
    
    Q = new K_RATPOLY(derivative(0));
    
    if (Q->deg[0] || sgn(Q->coeffs[0]))
    {
      P            = this;
      P->Sturm_seq = Q;
      P->Sturm_seq->ref_count++;
      
      while ((R = new K_RATPOLY(rem(*P, *Q)))->deg[0])
      {
        P = Q;
        Q = new K_RATPOLY(R->neg());
        delete R;
        
        P->Sturm_seq = Q;
        P->Sturm_seq->ref_count++;
      }
      
      if (sgn(R->coeffs[0]))
      {
        Q->Sturm_seq = new K_RATPOLY(R->neg());
        Q->Sturm_seq->ref_count++;
        Q->Sturm_seq->Sturm_seq = 0;
      }
      else  //  if *R is a zero polynomial
        Q->Sturm_seq = 0;
      
      delete R;
    }
    else  //  if this->derivative(0) is a zero polynomial
      delete Q;
  }
  
//  cerr << endl << " kratpoly: set_Sturm_seq: -------------------- " << endl << flush;
//  
//  unsigned long i;
//  K_RATPOLY*    S;
//  
//  i = 0;
//  S = this;
//  
//  while (S)
//  {
//    cerr << i++ << "th = " << endl << flush;
//    cerr << *S << flush;
//    cerr << " ref_count = " << S->ref_count << endl << endl << flush;
//    S = S->Sturm_seq;
//  }
//  
//  cerr << endl << " kratpoly: set_Sturm_seq: ==================== " << endl << flush;
  
  return Sturm_seq;
}

//  long K_RATPOLY :: num_Sturm_seq_perm(const bigrational& b) const
//    PROVIDED *this is a univariate polynomial,
//    returns
//      the number of sign permanencies in the Sturm sequence of *this at b.

long K_RATPOLY :: num_Sturm_seq_perm(const bigrational& b) const
{
  assert(num_vars == 1);
  
  long count;
  
  if (Sturm_seq)
  {
    int        prev_sgn, curr_sgn;
    K_RATPOLY* curr;
    K_RATPOLY* T;
    
    count    = 0;
    prev_sgn = sgn_at(b);
    curr     = Sturm_seq;  //  Sturm_seq != 0 => curr != 0
    
    if (!prev_sgn)  //  if b is a root of *this
    {
      curr_sgn = curr->sgn_at(b);
      
      if (!curr_sgn)
      //  If b is a common root of *this and its first derivative,
      //  i.e., b is a root of *this of multiplicity > 1
      //  then remove multiplicity.
      {
        K_RATPOLY B;
        
        T = new K_RATPOLY(*this);
        B = K_RATPOLY(1, 0, b);
        
        while (!curr_sgn)
        {
          *T       = T->exact_div(B);
          curr_sgn = T->derivative(0).sgn_at(b);
        }
        
        curr = T->set_Sturm_seq();
      }
      else  //  if b is a root of *this and not a root of its first derivative.
        T = 0;
      
      prev_sgn = curr_sgn;
      curr     = curr->Sturm_seq;
    }
    else  //  if b is not a root of *this
      T = 0;
    
    assert(prev_sgn);
    
    while (curr)
    {
      curr_sgn = curr->sgn_at(b);
      
      if (curr_sgn > 0  && prev_sgn > 0
          ||
          curr_sgn == 0 && prev_sgn > 0
          ||
          curr_sgn < 0  && prev_sgn < 0)
        count++;
      else  //  if (curr_sgn * prev_sgn < 0)
        prev_sgn = curr_sgn ? curr_sgn : 1;  //  prev_sgn != 0
      
      curr = curr->Sturm_seq;
    }
    
    if (T)
      delete T;
  }
  else  //  if (!Sturm_seq)
    count = - 1;
  
  return count;
}

//  K_RATPOLY K_RATPOLY :: gcd1(const K_RATPOLY& P) const
//    PROVIDED *this and P are univariate polynomials,
//    return their gcd computed by using Euclidean algorithm.

K_RATPOLY K_RATPOLY :: gcd1(const K_RATPOLY& P) const
{
  assert(num_vars == 1);
  assert(P.num_vars == 1);
  
  unsigned long i, j;
  K_RATPOLY     G;
  
  if (!sgn(coeffs[0]))            //  if *this is 0
    G = P;
  else if (!sgn(P.coeffs[0]))     //  if P is 0
    G = *this;
  else if (!deg[0] || !P.deg[0])  //  if *this or P is constant
  {
    long dg[1];
    
    dg[0]       = 0;
    G           = K_RATPOLY(1, dg);
    G.coeffs[0] = 1;
  }
  else
  //  compute gcd of *this and P by using Euclidean algorithm
  {
    K_RATPOLY Q0, Q1, Q2;
    
    if (deg[0] >= P.deg[0])
    {
      Q0 = *this;
      Q1 = P;
    }
    else  //  if (deg[0] < P.deg[0])
    {
      Q0 = P;
      Q1 = *this;
    }
    
    while ((Q2 = rem(Q0, Q1)).deg[0])
    {
      Q0 = Q1;
      Q1 = Q2;
    }
    
    if (sgn(Q2.coeffs[0]))
      G = Q2;
    else  //  if Q2 is a zero polynomial
      G = Q1;
  }
  
  return G;
}

//  K_RATPOLY gcd(const K_RATPOLY& P1, const K_RATPOLY& P2)
//    PROVIDED P1 and P2 are uni or bivariate polynomials,
//    returns their gcd.

K_RATPOLY gcd(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  assert(P1.num_vars == 1 && P2.num_vars == 1
         ||
         P1.num_vars == 2 && P2.num_vars == 2);
  
  K_RATPOLY G;
  
  if (P1.num_vars == 1 && P2.num_vars == 1)
    G = P1.gcd1(P2);
  else  //  if P1 and P2 are bivariate polynomials
    G = P1.gcd2(P2);
  
  return G;
}

//  bigrational K_RATPOLY :: Sylvester1(const K_RATPOLY& P) const
//    PROVIDED *this and P are univariate polynomials,
//    returns Sylvester resultant for *this and P.

bigrational K_RATPOLY :: Sylvester1(const K_RATPOLY& P) const
{
  assert(num_vars == 1);
  assert(P.num_vars == 1);
  
  unsigned long i, j;
  unsigned long n;
  bigrational   r;
  
  if (n = deg[0] + P.deg[0])
  {
    bigrational_matrix S(n, n);
    
    for (i = 0; i < P.deg[0]; i++)
      for (j = 0; j <= deg[0]; j++)
        S(i, i + j) = coeffs[j];
    
    for (i = 0; i < deg[0]; i++)
      for (j = 0; j <= P.deg[0]; j++)
        S(P.deg[0] + i, i + j) = P.coeffs[j];
    
    r = det(S);
  }
  else if (sgn(coeffs[0]) || sgn(P.coeffs[0]))
  //  if *this or P is a non-zero constant polynomial
    r = 1;
  else
  //  if both *this and P are zero polynomials
    r = 0;
  
  return r;
}

bigrational Sylvester1(const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  return P1.Sylvester1(P2);
}

//  bigrational_vector gen_VDM_row1(n)
//    returns the first row of Vandermonde matrix, i.e.,
//            a vector of n distinct rational numbers.

bigrational_vector gen_VDM_row1(const unsigned long n)
{
  long               i;
  bigrational_vector R(n);
  
  for (i = 0; i < n / 2; i++)
  {
    R[2 * i]     = i + 1;
    R[2 * i + 1] = - i - 1;
  }
  
  if (n % 2)
    R[n - 1] = n / 2 + 1;
  
  return R;
}

//  K_RATPOLY K_RATPOLY :: Sylvester2(const K_RATPOLY& P,
//                                    const unsigned long i) const
//    PROVIDED *this and P are bivariate polynomials,
//    returns Sylvester resultant for *this and P w.r.t. X_i.

K_RATPOLY K_RATPOLY :: Sylvester2(const K_RATPOLY& P,
                                  const unsigned long i) const
{
  assert(num_vars == 2);
  assert(P.num_vars == 2);
  assert(i == 0 || i == 1);
  
  static unsigned long      num_vdm_rows = 128;
  static bigrational_vector U0           = gen_VDM_row1(num_vdm_rows);
  
  unsigned long j, k;
  unsigned long n;
  int           singular;
  K_RATPOLY     T1, P1;
  long          dt1, dp1;
  long          dr[1];
//  K_RATPOLY     R;
  
  if ((n = deg[0] * P.deg[1] + deg[1] * P.deg[0] + 1) > num_vdm_rows)
  {
    U0 = gen_VDM_row1(num_vdm_rows = n);
//    cerr << " kratpoly: Sylvester2: num_vdm_rows = " << num_vdm_rows << endl << flush;
  }
  
  bigrational_vector RHS(n);
  bigrational_vector U(n);
  bigrational_vector Sol;
  
  for (j = 0; j < n; j++)
  {
    singular = 1;
    
    while (singular)
    {
      T1  = subst_val_proto(1 - i, U0[j], 0);
      P1  = P.subst_val_proto(1 - i, U0[j], 0);
      dt1 = T1.deg[0];
      dp1 = P1.deg[0];
      T1.reduce_deg();
      P1.reduce_deg();
      
      if (dt1 == T1.deg[0] && dp1 == P1.deg[0])
        singular = 0;
      else  //  if (dt1 > T1.deg[0] || dp1 > P1.deg[0])
        for (k = j + 1; k < num_vdm_rows; k++)
          U0[k - 1] = U0[k];
    }
    
    RHS[j] = T1.Sylvester1(P1);
  }
  
  for (j = 0; j < n; j++)
    U[j] = U0[j];
  
  Sol = Solve_VDM_GE(U, RHS);
  
  dr[0] = n - 1;
  
//  R = K_RATPOLY(1, dr);
  K_RATPOLY R(1, dr);
  
  for (j = 0; j < n; j++)
    R.coeffs[j] = Sol[n - j - 1];
  
  R.reduce_deg();
  R.reduce_num_coeffs();
  
  return R;
}

//  K_RATPOLY K_RATPOLY :: GoodSylvester(const K_RATPOLY& P,
//                                       const unsigned long i) const
//    PROVIDED *this and P are bivariate polynomials,
//    returns Sylvetser resultant for *this and P w.r.t. X_i.

K_RATPOLY K_RATPOLY :: GoodSylvester(const K_RATPOLY& P,
                                     const unsigned long i) const
{
  assert(num_vars == 2);
  assert(P.num_vars == 2);
  assert(i == 0 || i == 1);
  
  K_RATPOLY R;
  
  if (!deg[i] && !P.deg[i])
  //  if *this and P are actually univariate polynomials only in X_{1 - i}
    R = remove_var(i).gcd1(P.remove_var(i));
  else if (!deg[1 - i] && !P.deg[1 - i])
  //  if *this and P are actually univariate polynomials only in X_i.
  {
    long dr[1];
    
    dr[0]       = 0;
    R           = K_RATPOLY(1, dr);
    R.coeffs[0] = remove_var(1 - i).Sylvester1(P.remove_var(1 - i));
  }
  else  //  if *this and P are non-trivial bivariate polynomials
  {
    R = Sylvester2(P, i);
    
    if (!R.deg[0] && !sgn(R.coeffs[0]))
    //  If a univariate polynomial R is 0, i.e.,
    //  if *this and P have a non-trivial common divisor.
//    {
//      bigrational y(3, 7);
//      
//      K_RATPOLY T1 = subst_val(1 - i, y);
//      K_RATPOLY P1 = P.subst_val(1 - i, y);
//      
//      K_RATPOLY G1 = T1.gcd1(P1);
//      K_RATPOLY G2 = G1.add_var(1 - i);
//      
//      K_RATPOLY T2 = exact_div(G2);
//      K_RATPOLY P2 = P.exact_div(G2);
//      
//      T2.reduce_num_coeffs();
//      P2.reduce_num_coeffs();
//      
//      if (!T2.deg[i] && !P2.deg[i])
//      //  If univariate polynomials T2 and P2 are both 0.
//        R = T2.remove_var(i).gcd1(P2.remove_var(i));
//      else
//        R = T2.Sylvester2(P2, i);
//      cerr << " kratpoly: GoodSylvester: refined Res = " << endl << R << endl << flush;
//    }
    {
      K_RATPOLY G, T1, P1;
      
      G = gcd2(P);
//      cerr << " kratpoly: GoodSylvester: G = " << endl << G << endl << flush;
      
      T1 = exact_div(G);
      P1 = P.exact_div(G);
      T1.reduce_num_coeffs();
      P1.reduce_num_coeffs();
//      cerr << " kratpoly: GoodSylvester: T1 = " << endl << T1 << endl << flush;
//      cerr << " kratpoly: GoodSylvester: P1 = " << endl << P1 << endl << flush;
      
      R = T1.Sylvester2(P1, i);
//      cerr << " kratpoly: GoodSylvester: refined Res = " << endl << R << endl << flush;
    }
  }
  
  return R;
}

K_RATPOLY GoodSylvester(const K_RATPOLY& P1, const K_RATPOLY& P2,
                        const unsigned long i)
{
  return P1.GoodSylvester(P2, i);
}

//  K_RATPOLY K_RATPOLY :: monic() const
//    PROVIDED *this is a univariate polynomial,
//    returns a monic polynomial obtained by
//      dividing all the coeffs of *this by the leading coeff *this.

K_RATPOLY K_RATPOLY :: monic() const
{
  assert(num_vars == 1);
  
  unsigned long i;
  bigrational   lc;
//  K_RATPOLY     P;
  
//  P = K_RATPOLY(1, deg);
  K_RATPOLY P(1, deg);
  
  if (sgn(lc = coeffs[0]))
  {
    for(i = 0; i < num_coeffs; i++)
      if (sgn(coeffs[i]))
        P.coeffs[i] = coeffs[i] / lc;
  }
  else  //  if *this is a zero polnomial
    P.coeffs[0] = 0;
  
  return P;
}

//  K_RATPOLY K_RATPOLY :: monic_gcd1(const K_RATPOLY& P) const
//    PROVIDED *this and P are univariate polynomials,
//    returns monic gcd of *this and P.

K_RATPOLY K_RATPOLY :: monic_gcd1(const K_RATPOLY& P) const
{
  assert(num_vars == 1);
  assert(P.num_vars == 1);
  
  K_RATPOLY G;
  
  if (!sgn(coeffs[0]))            //  if *this is s zero polynomial
    G =  P.monic();
  else if (!sgn(P.coeffs[0]))     //  if P is a zero polynomial
    G = monic();
  else if (!deg[0] || !P.deg[0])  //  if *this or P is constant
  {
    long dg[1];
    
    dg[0]       = 0;
    G           = K_RATPOLY(1, dg);
    G.coeffs[0] = 1;
  }
  else  //  Non-trivial cases: perform Euclidean algorithm
  {
    K_RATPOLY Q0, Q1, Q2;
    
    if (deg[0] >= P.deg[0])
    {
      Q0 = monic();
      Q1 = P.monic();
    }
    else  //  if (deg[0] < P.deg[0])
    {
      Q0 = P.monic();
      Q1= monic();
    }
    
    while ((Q2 = rem(Q0, Q1)).deg[0])
    {
      Q0 = Q1;
      Q1 = Q2.monic();
    }
    
    if (sgn(Q2.coeffs[0]))
      G = Q2.monic();
    else  //  if Q2 is a zero polynomial
      G = Q1;
  }
  
  return G;
}

//  K_RATPOLY K_RATPOLY :: gcd2_pp(const K_RATPOLY& P) const
//     PROVIDED *this and P are bivariate primitive polynomials,
//     returns gcd(*this, P).

K_RATPOLY K_RATPOLY :: gcd2_pp(const K_RATPOLY& P) const
{
  assert(num_vars == 2);
  assert(P.num_vars == 2);
  
  unsigned long i, j, k;
  long          d[2];
  long          p[2];
  int           c0, c1;
  long          dlc[1];
  long          plc[1];
  K_RATPOLY     LCT, LCP, B;
  unsigned long dx, dy, l;
  unsigned long s0, s1, s2;
  bigrational*  S;
  bigrational   u0;
  bigrational   u[1];
  K_RATPOLY*    V_S;
  unsigned long e, dT, dP;
  int*          good_S;
  K_RATPOLY     W, CW, PPW;
  
//  cerr << " kratpoly: gcd2_pp: *this = " << endl << *this << endl << flush;
//  cerr << " kratpoly: gcd2_pp: P = " << endl << P << endl << flush;
  
  //  1. Compute B = gcd(lc_X(*this), lc_X(P)) and
  //             l = max{ deg_Y *this, deg_Y P } + 1 + deg_Y B.
  
  dlc[0] = deg[1];
  LCT    = K_RATPOLY(1, dlc);
  
  for (i = 0; i <= deg[1]; i++)
  {
    p[0]               = deg[0];
    p[1]               = i;
    plc[0]             = i;
    LCT.get_coeff(plc) = get_coeff(p);
  }
  
  LCT.reduce_deg();
  
  dlc[0] = P.deg[1];
  LCP    = K_RATPOLY(1, dlc);
  
  for (i = 0; i <= P.deg[1]; i++)
  {
    p[0]               = P.deg[0];
    p[1]               = i;
    plc[0]             = i;
    LCP.get_coeff(plc) = P.get_coeff(p);
  }
  
  LCP.reduce_deg();
  
  B = LCT.monic_gcd1(LCP);
  
  dx = deg[0] > P.deg[0] ? deg[0] : P.deg[0];
  dy = deg[1] > P.deg[1] ? deg[1] : P.deg[1];
  l  = B.deg[0] + dy + 1;                      //  l > 0
  
//  cerr << " kratpoly: gcd2_pp: LCT = " << endl << LCT << endl << flush;
//  cerr << " kratpoly: gcd2_pp: LCP = " << endl << LCP << endl << flush;
//  cerr << " kratpoly: gcd2_pp: B = " << endl << B << endl << flush;
//  cerr << " kratpoly: gcd2_pp: dx = " << dx << ", dy = " << dy << ", l = " << l << endl << flush;
  
  S  = new bigrational [s0 = 2 * l];  //  s0 > l > 0
  u0 = 0;
  c0 = 1;
  
  while (c0)
  {
    //  4-1.  Compute S = { u | B(u) != 0 } and s1 = #S.
    
    for (i = s1 = 0; i < s0; i++)
    {
      u[0] = u0;
      
      if (sgn(B.evaluate(u)))
        S[s1++] = u0;
      
      u0 = u0 + 1;
    }
    
//    cerr << " kratpoly: gcd2_pp: s1 = " << s1 << ", S = { ";
//    if (s1)
//    {
//      for (i = 0; i < s1 - 1; i++)
//        cerr << S[i]  << ", ";  
//      cerr << S[s1 - 1];
//    }
//    cerr << " }" << endl << flush;
    
    if (s1 >= l)
    {
      //  4-2.  Compute V_S[i] = gcd(*this(X, S[i]), P(X, S[i])).
      
      V_S = new K_RATPOLY [s1];  //  s1 >= l > 0
      e   = dx;
      
      for (i = 0; i < s1; i++)
      {
        V_S[i] = subst_val(1, S[i]).monic_gcd1(P.subst_val(1, S[i]));
        
        if (e > V_S[i].deg[0])
          e = V_S[i].deg[0];
      }
      
//      for (i = 0; i < s1; i++)
//        cerr << " kratpoly: gcd2_pp: V_S[" << i << "] = " << endl << V_S[i] << endl << flush;
//      cerr << " kratpoly: gcd2_pp: e = " << e << endl << flush;
      
      //  5.  Compute e = min{ deg V_S[i] } and
      //      set good_S[i] = 1 iff deg V_S[i] == e.
      
      good_S = new int [s1];  //  s1 >= l > 0
      
      for (i = s2 = 0 ; i < s1; i++)
        if (e == V_S[i].deg[0])
        {
          good_S[i] = 1;
          s2++;
        }
        else
          good_S[i] = 0;
      
//      cerr << " kratpoly: gcd2_pp: s2 = " << s2 << ", good_S = { ";
//      for (i = 0; i < s1 - 1; i++)
//        cerr << good_S[i]  << ", ";
//      cerr << good_S[s1 - 1] << " }" << endl << flush;
      
      if (s2 >= l)
      {
        //  6.  Use interpolation to compute the polynomials
        //        W         in variable X and Y, and
        //        T2 and P2 in variable X
        //      s.t.
        //        W      = B(S[i]) * V_S[i],
        //        T2 * W = B(S[i]) * *this(X, S[i]) and
        //        P2 * W = B(S[i]) * P(X, S[i]).
        
        dT = deg[0] - e;
        dP = P.deg[0] - e;
        
        bigrational_vector U(l);
        bigrational_matrix WU(l, e + 1);
        bigrational_matrix TU(l, dT + 1);
        bigrational_matrix PU(l, dP + 1);
        
        c1 = 1;
        
        for (i = j = 0; c1 && i < s1 && j < l; i++)
          if (good_S[i])
          {
            u[0] = U[j] = S[i];
            
            K_RATPOLY T1 = subst_val(1, S[i]).exact_div(V_S[i]);
            K_RATPOLY P1 = P.subst_val(1, S[i]).exact_div(V_S[i]);
            
//            cerr << " kratpoly: gcd2_pp: u = " << u[0] << endl << flush;
//            cerr << " kratpoly: gcd2_pp: T1 = " << endl << T1 << endl << flush;
//            cerr << " kratpoly: gcd2_pp: T1.deg[0] = " << T1.deg[0] << ", dT = " << dT << endl << flush;
//            cerr << " kratpoly: gcd2_pp: P1 = " << endl << P1 << endl << flush;
//            cerr << " kratpoly: gcd2_pp: P1.deg[0] = " << P1.deg[0] << ", dP = " << dP << endl << flush;
            
            if (T1.deg[0] == dT && P1.deg[0] == dP)
            {
              for (k = 0; k <= e; k++)
              {
                p[0]     = k;
                WU(j, k) = B.evaluate(u) * V_S[i].get_coeff(p);
              }
              
              for (k = 0; k <= dT; k++)
              {
                p[0]     = k;
                TU(j, k) = T1.get_coeff(p);
              }
              
              for (k = 0; k <= dP; k++)
              {
                p[0]     = k;
                PU(j, k) = P1.get_coeff(p);
              }
              
              j++;
            }
            else
              c1 = 0;
          }
        
        if (c1)
        {
          d[0] = e;
          d[1] = l - 1;
          
          W = K_RATPOLY(2, d);
          
          for (i = 0; i <= e; i++)
          {
            bigrational_vector WL(l);
            bigrational_vector WR(l);
            
            for (j = 0; j < l; j++)
              WR[j] = WU(j, i);
            
            WL = Solve_VDM_GE(U, WR);
            
            for (j = 0; j < l; j++)
            {
              p[0]           = i;
              p[1]           = j;
              W.get_coeff(p) = WL[j];
            }
          }
          
          W.reduce_deg();
//          cerr << " kratpoly: gcd2_pp: W = " << endl << W << endl << flush;
          
          d[0] = dT;
          d[1] = l - 1;
          
          K_RATPOLY T2(2, d);
          
          for (i = 0; i <= dT; i++)
          {
            bigrational_vector TL(l);
            bigrational_vector TR(l);
            
            for (j = 0; j < l; j++)
              TR[j] = TU(j, i);
            
            TL = Solve_VDM_GE(U, TR);
            
            for (j = 0; j < l; j++)
            {
              p[0]            = i;
              p[1]            = j;
              T2.get_coeff(p) = TL[j];
            }
          }
          
          T2.reduce_deg();
          
          K_RATPOLY T2_W = T2 * W;
          K_RATPOLY B_T  = B.add_var(0) * *this;
          
//          cerr << " kratpoly: gcd2_pp: T2 = " << endl << T2 << endl << flush;
//          cerr << " kratpoly: gcd2_pp: T2_W = " << endl << T2_W << endl << flush;
//          cerr << " kratpoly: gcd2_pp: B_T = " << endl << B_T << endl << flush;
//          cerr << " kratpoly: gcd2_pp: deg_Y T2_W = " << T2_W.deg[1] << ", deg_Y B_T = " << B_T.deg[1] << endl << flush;
          
          d[0] = dP;
          d[1] = l - 1;
          
          K_RATPOLY P2(2, d);
          
          for (i = 0; i <= dP; i++)
          {
            bigrational_vector PL(l);
            bigrational_vector PR(l);
            
            for (j = 0; j < l; j++)
              PR[j] = PU(j, i);
            
            PL = Solve_VDM_GE(U, PR);
            
            for (j = 0; j < l; j++)
            {
              p[0]            = i;
              p[1]            = j;
              P2.get_coeff(p) = PL[j];
            }
          }
          
          P2.reduce_deg();
          
          K_RATPOLY P2_W = P2 * W;
          K_RATPOLY B_P  = B.add_var(0) * P;
          
//          cerr << " kratpoly: gcd2_pp: P2 = " << endl << P2 << endl << flush;
//          cerr << " kratpoly: gcd2_pp: P2_W = " << endl << P2_W << endl << flush;
//          cerr << " kratpoly: gcd2_pp: B_P = " << endl << B_P << endl << flush;
//          cerr << " kratpoly: gcd2_pp: deg_Y P2_W = " << P2_W.deg[1] << ", deg_Y B_P = " << B_P.deg[1] << endl << flush;
          
          //  7.  If deg_Y (T2 * W) == deg_Y (B * *this) and
          //         deg_Y (P2 * w) == deg_Y (B * P) then exit the loop.
          
          if (T2_W.deg[1] == B_T.deg[1] && P2_W.deg[1] == B_P.deg[1])
            c0 = 0;
        }
      }
      
      delete [] good_S;  //  s1 >= l > 0 => good_S != 0
      delete [] V_S;     //  s1 >= l > 0 => V_S != 0
    }
  }
  
  delete [] S;  //  s0 > 0 => S != 0
  
  //  8.  Compute pp_X (W).
  
  K_RATPOLY CW0[W.deg[0] + 1];
  
  dlc[0] = W.deg[1];
  
  for (i = 0; i <= W.deg[0]; i++)
  {
    CW0[i] = K_RATPOLY(1, dlc);
    
    for (j = 0; j <= W.deg[1]; j++)
    {
      p[0]                  = i;
      p[1]                  = j;
      plc[0]                = j;
      CW0[i].get_coeff(plc) = W.get_coeff(p);
    }
    
    CW0[i].reduce_deg();
  }
  
//  for (i = 0; i <= W.deg[0]; i++)
//    cerr << " kratpoly: gcd2_pp: CW0[" << i << "] = " << endl << CW0[i] << endl << flush;
  
  CW = CW0[0];
  
  for (i= 1; i <= W.deg[0]; i++)
    CW = CW.monic_gcd1(CW0[i]);
  
  PPW = W.exact_div(CW.add_var(0));
  
//  cerr << " kratpoly: gcd2_pp: CW = " << endl << CW << endl << flush;
//  cerr << " kratpoly: gcd2_pp: PPW = " << endl << PPW << endl << flush;
  
  return PPW;
}

//  K_RATPOLY K_RATPOLY :: gcd2(const K_RATPOLY& P) const
//    PROVIDED *this and P are bivariate polynomials,
//    returns gcd(*this, P).

K_RATPOLY K_RATPOLY :: gcd2(const K_RATPOLY& P) const
{
  assert(num_vars == 2);
  assert(P.num_vars == 2);
  
  unsigned long i, j;
  long          dc[1];
  long          p[2];
  long          pc[1];
  
//  cerr << " kratpoly: gcd2: *this = " << endl << *this << endl << flush;
//  cerr << " kratpoly: gcd2: P = " << endl << P << endl << flush;
  
  //  1-1. Compute the content C and primitive part PP of *this.
  
  K_RATPOLY C0[deg[0] + 1];
  
  dc[0] = deg[1];
  
  for (i = 0; i <= deg[0]; i++)
  {
    C0[i] = K_RATPOLY(1, dc);
    
    for (j = 0; j <= deg[1]; j++)
    {
      p[0]                = i;
      p[1]                = j;
      pc[0]               = j;
      C0[i].get_coeff(pc) = get_coeff(p);
    }
    
    C0[i].reduce_deg();
//    cerr << " kratpoly: gcd2: C0[" << i << "] = " << endl << C0[i] << endl << flush;
  }
  
  K_RATPOLY C(C0[0]);
  
  for (i= 1; i <= deg[0]; i++)
    C = C.monic_gcd1(C0[i]);
  
//  cerr << " kratpoly: gcd2: C = " << endl << C << endl << flush;
  
  K_RATPOLY PP = exact_div(C.add_var(0));
  
//  cerr << " kratpoly: gcd2: PP = " << endl << PP << endl << flush;
  
  //  1-2. Compute the content CP and primitive part PPP of P.
  
  K_RATPOLY CP0[P.deg[0] + 1];
  
  dc[0] = P.deg[1];
  
  for (i = 0; i <= P.deg[0]; i++)
  {
    CP0[i] = K_RATPOLY(1, dc);
    
    for (j = 0; j <= P.deg[1]; j++)
    {
      p[0]                 = i;
      p[1]                 = j;
      pc[0]                = j;
      CP0[i].get_coeff(pc) = P.get_coeff(p);
    }
    
    CP0[i].reduce_deg();
//    cerr << " kratpoly: gcd2: CP0[" << i << "] = " << endl << CP0[i] << endl << flush;
  }
  
  K_RATPOLY CP(CP0[0]);
  
  for (i= 1; i <= P.deg[0]; i++)
    CP = CP.monic_gcd1(CP0[i]);
  
//  cerr << " kratpoly: gcd2: CP = " << endl << CP << endl << flush;
  
  K_RATPOLY PPP = P.exact_div(CP.add_var(0));
  
//  cerr << " kratpoly: gcd2: PPP = " << endl << PPP << endl << flush;
  
  //  2-1. Compute monic_gcd CG of C and CP.
  
  K_RATPOLY CG = C.monic_gcd1(CP);
  
//  cerr << " kratpolyY: gcd2: CG = " << endl << CG << endl << flush;
  
  //  2-2. Compute gcd PPG of PP and PPP.
  
  K_RATPOLY PPG = PP.gcd2_pp(PPP);
  
//  cerr << " kratpoly: gcd2: PPG = " << endl << PPG << endl << flush;
  
  //  2-3. Compute gcd g of *this and P.
  
  K_RATPOLY G = CG.add_var(0) * PPG;
  
//  cerr << " kratpoly: gcd2: G = " << endl << G << endl << flush;
  
  return G;
}

//  int K_RATPOLY :: eval_range(const bigrational* const low_in,
//                              const bigrational* const high_in,
//                              bigrational& low_out,
//                              bigrational& high_out)
//    computes the range [low_out, high_out]
//    which *this is evaluated at some value in [low_in, high_in]
//    using affine arithmetic.

int K_RATPOLY :: eval_range(const bigrational* const low_in,
                            const bigrational* const high_in,
                            bigrational& low_out,
                            bigrational& high_out)
{
  unsigned long i, j;
  long          d[] = { 1 };
  K_RATPOLY     N;
  K_RATPOLY     T(*this);
  bigrational   err;
  
  for (i = 0; i < num_vars; i++)
  {
    N           = K_RATPOLY(1, d);
    N.coeffs[0] = (high_in[i] - low_in[i]) / 2;
    N.coeffs[1] = (low_in[i] + high_in[i]) / 2;
    T           = T.add_var(i);
    
    for (j = 0; j < i; j++)
      N = N.add_var(j);
    
    for (j = i + 1; j < num_vars; j++)
      N = N.add_var(j);
    
    T = T.subst_expr(i + 1, N);
  }
  
  err = 0;
  
  for (i = 0; i < T.num_coeffs - 1; i++)
    err += abs(T.coeffs[i]);
  
  low_out  = T.coeffs[T.num_coeffs - 1] - err;
  high_out = T.coeffs[T.num_coeffs - 1] + err;
  
  return 0;
}

//  int K_RATPOLY :: eval_range(const bigrational_vector& low_in,
//                              const bigrational_vector& high_in,
//                              bigrational& low_out,
//                              bigrational& high_out)
//    computes the range [low_out, high_out]
//    which *this is evaluated at some value in [low_in, high_in]
//    using affine arithmetic.

int K_RATPOLY :: eval_range(const bigrational_vector& low_in,
                            const bigrational_vector& high_in,
                            bigrational& low_out,
                            bigrational& high_out)
{
  assert(num_vars == low_in.get_dim());
  assert(num_vars == high_in.get_dim());
  
  unsigned long i, j;
  long          d[] = { 1 };
  K_RATPOLY     N;
  K_RATPOLY     T(*this);
  bigrational   err;
  
  for (i = 0; i < num_vars; i++)
  {
    N           = K_RATPOLY(1, d);
    N.coeffs[0] = (high_in[i] - low_in[i]) / 2;
    N.coeffs[1] = (low_in[i] + high_in[i]) / 2;
    T           = T.add_var(i);
    
    for (j = 0; j < i; j++)
      N = N.add_var(j);
    
    for (j = i + 1; j < num_vars; j++)
      N = N.add_var(j);
    
    T = T.subst_expr(i + 1, N);
  }
  
  err = 0;
  
  for (i = 0; i < T.num_coeffs - 1; i++)
    err += abs(T.coeffs[i]);
  
  low_out  = T.coeffs[T.num_coeffs - 1] - err;
  high_out = T.coeffs[T.num_coeffs - 1] + err;
  
  return 0;
}

//  bigrational K_RATPOLY :: get_Mignotte_bd() const
//    PROVIDED *this is a univariate polynomial,
//    returns an upper bound on the size of the largest root of *this.

bigrational K_RATPOLY :: get_Mignotte_bd() const
{
  assert(num_vars == 1);
  
  unsigned long i;
  bigrational   lc, mac, ac;
  bigrational   M;
  
  lc = M = 0;
  
  for (i = 0; i < num_coeffs && !sgn(lc); i++)
    if (sgn(coeffs[i]))
      lc = coeffs[i];
  
  if (i < num_coeffs)
  {
    mac = 0;
    
    for (; i < num_coeffs; i++)
      if (sgn(coeffs[i]) && (ac = abs(coeffs[i])) > mac)
        mac = ac;
    
    M = 1 + mac / abs(lc);
  }
  
  return M;
}

int K_RATPOLY :: reduce_num_coeffs()
{
  unsigned long i, j;
  bigint        d, g;
  
  d = 1;
  
  for (i = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]))
      d *= den(coeffs[i]);
  
//  if (d != 1)
  {
    for (i = 0; i < num_coeffs; i++)
      if (sgn(coeffs[i]))
        coeffs[i] *= d;
    
    for (i = 0; i < num_coeffs && !sgn(coeffs[i]); i++)
      ;
    
    g = 1;
    
    if (i < num_coeffs)
    {
      g = num(coeffs[i]);
      
      for (j = i + 1; j < num_coeffs && g != 1; j++)
        if (sgn(coeffs[j]))
          g = gcd(g, num(coeffs[j]));
      
      if (g != 1)
        for (j = 0; j < num_coeffs; j++)
          if (sgn(coeffs[j]))
            coeffs[j] /= g;
    }
  }
  
  return d != 1;
}

int K_RATPOLY :: reduce_coeffs()
{
  unsigned long i, j;
  bigint        g, l, gl;
  
  g = l = 1;
  
  for (i = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]))
    {
      g = gcd(g, num(coeffs[i]));
      l = lcm(l, den(coeffs[i]));
    }
  
  if ((gl = g * l) != 1)
    for (i = 0; i < num_coeffs; i++)
      if (sgn(coeffs[i]))
        coeffs[i] = num(coeffs[i] / gl);
  
  return gl != 1;
}

K_RATPOLY K_RATPOLY :: conv_to_Bernstein(const long deg_s,
                                         const long deg_t) const
{
  assert(num_vars == 2);
  
  long          i, j, k;
  long          m, n;
  bigint        num, den;
  long*         p;
  long          d[2];
//  K_RATPOLY     B;
  
  m = deg[0];
  n = deg[1];
  
  bigrational_matrix M(m + 1, m + 1);
  bigrational_matrix N(n + 1, n + 1);
  bigrational_matrix C(m + 1, n + 1);
  
  for (i = 0; i <= m; i++)
    for (j = i; j <= m; j++)
    {
      num = 1;
      
      for (k = j; k > i; k--)
        num *= k;
      
      for (k = m - i; k > 1; k--)
        num *= k;
      
      den = 1;
      
      for (k = m; k > i; k--)
        den *= k;
      
      for (k = j - i; k > 1; k--)
        den *= k;
      
      M(j, i) = bigrational(num, den);
    }
  
  for (i = 0; i <= n; i++)
    for (j = i; j <= n; j++)
    {
      num = 1;
      
      for (k = j; k > i; k--)
        num *= k;
      
      for (k = n - i; k > 1; k--)
        num *= k;
      
      den = 1;
      
      for (k = n; k > i; k--)
        den *= k;
      
      for (k = j - i; k > 1; k--)
        den *= k;
      
      N(i, j) = bigrational(num, den);
    }
  
  for (k = 0; k < num_coeffs; k++)
  {
    p             = index_to_powers(k);
    C(p[0], p[1]) = coeffs[k];
    delete [] p;
  }
  
  bigrational_matrix R, Raise;
  
  R = M * C * N;
  
  for (; m < deg_s; m++)
  {
    Raise = bigrational_matrix(m + 2, m + 1);
    
    for (i = 0; i <= m + 1; i++)
    {
      if (i > 0)
        Raise(i, i - 1) = bigrational(i, m + 1);
      
      if (i < m + 1)
        Raise(i, i) = bigrational(m + 1 - i, m + 1);
    }
    
    R = Raise * R;
  }
  
  for (; n < deg_t; n++)
  {
    Raise = bigrational_matrix(n + 1, n + 2);
    
    for (i = 0; i <= n + 1; i++)
    {
      if (i > 0)
        Raise(i - 1, i) = bigrational(i, n + 1);
      
      if (i < n + 1)
        Raise(i, i) = bigrational(n + 1 - i, n + 1);
    }
    
    R = R * Raise;
  }
  
  d[0] = deg_s;
  d[1] = deg_t;
  
//  B = K_RATPOLY(num_vars, d);
  K_RATPOLY B(num_vars, d);
  
  p = new long [2];
  
  for (i = 0; i <= deg_s; i++)
    for (j = 0; j <= deg_t; j++)
    {
      p[0]                         = i;
      p[1]                         = j;
      B.coeffs[B.powers_to_index(p)] = R(i, j);
    }
  
  delete [] p;  //  p != 0
  
  return B;
}

K_RATPOLY K_RATPOLY :: transform_Impl(const bigrational_matrix& T) const
{
  assert(num_vars == 3);
  
  unsigned long i;
  long          d[3];
  long          p[3];
  long          t;
  K_RATPOLY*    Xp;
  K_RATPOLY*    Yp;
  K_RATPOLY*    Zp;
  K_RATPOLY*    Wp;
  long*         q;
//  K_RATPOLY     I;
  
  d[0] = d[1] = d[2] = 1;
  
  K_RATPOLY base_X(3, d);
  K_RATPOLY base_Y(3, d);
  K_RATPOLY base_Z(3, d);
  K_RATPOLY base_W(3, d);
  
  p[0]                                      = 1;
  p[1]                                      = 0;
  p[2]                                      = 0;
  base_X.coeffs[base_X.powers_to_index(p)] = T(0, 0);
  base_Y.coeffs[base_Y.powers_to_index(p)] = T(1, 0);
  base_Z.coeffs[base_Z.powers_to_index(p)] = T(2, 0);
  base_W.coeffs[base_W.powers_to_index(p)] = T(3, 0);
  
  p[0]                                      = 0;
  p[1]                                      = 1;
  p[2]                                      = 0;
  base_X.coeffs[base_X.powers_to_index(p)] = T(0, 1);
  base_Y.coeffs[base_Y.powers_to_index(p)] = T(1, 1);
  base_Z.coeffs[base_Z.powers_to_index(p)] = T(2, 1);
  base_W.coeffs[base_W.powers_to_index(p)] = T(3, 1);
  
  p[0]                                      = 0;
  p[1]                                      = 0;
  p[2]                                      = 1;
  base_X.coeffs[base_X.powers_to_index(p)] = T(0, 2);
  base_Y.coeffs[base_Y.powers_to_index(p)] = T(1, 2);
  base_Z.coeffs[base_Z.powers_to_index(p)] = T(2, 2);
  base_W.coeffs[base_W.powers_to_index(p)] = T(3, 2);
  
  p[0]                                      = 0;
  p[1]                                      = 0;
  p[2]                                      = 0;
  base_X.coeffs[base_X.powers_to_index(p)] = T(0, 3);
  base_Y.coeffs[base_Y.powers_to_index(p)] = T(1, 3);
  base_Z.coeffs[base_Z.powers_to_index(p)] = T(2, 3);
  base_W.coeffs[base_W.powers_to_index(p)] = T(3, 3);
  
  Xp = new K_RATPOLY [deg[0] + 1];
  Yp = new K_RATPOLY [deg[1] + 1];
  Zp = new K_RATPOLY [deg[2] + 1];
  Wp = new K_RATPOLY [(t = get_total_deg()) + 1];
  
  Xp[0] = K_RATPOLY(3, d);
  Yp[0] = K_RATPOLY(3, d);
  Zp[0] = K_RATPOLY(3, d);
  Wp[0] = K_RATPOLY(3, d);
  
  p[0]                                   = 1;
  p[1]                                   = 1;
  p[2]                                   = 1;
  Xp[0].coeffs[Xp[0].powers_to_index(p)] = 1;
  Yp[0].coeffs[Yp[0].powers_to_index(p)] = 1;
  Zp[0].coeffs[Zp[0].powers_to_index(p)] = 1;
  Wp[0].coeffs[Wp[0].powers_to_index(p)] = 1;
  
  for (i = 1; i <= deg[0]; i++)
    Xp[i] = base_X * Xp[i - 1];
  
  for (i = 1; i <= deg[1]; i++)
    Yp[i] = base_Y * Yp[i - 1];
  
  for (i = 1; i <= deg[2]; i++)
    Zp[i] = base_Z * Zp[i - 1];
  
  for (i = 1; i <= t; i++)
    Wp[i] = base_W * Wp[i - 1];
  
//  I = K_RATPOLY(3, d);
  K_RATPOLY I(3, d);
  
  for (i = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]))
    {
      q = index_to_powers(i);
      I = I +
          coeffs[i] * Xp[q[0]] * Yp[q[1]] * Zp[q[2]] *
          Wp[t - q[0] - q[1] - q[2]];
    }
  
  delete [] Xp;
  delete [] Yp;
  delete [] Zp;
  delete [] Wp;
  delete [] q;
  
  return I;
}

K_RATPOLY read_poly(istream& in_fs)
{
  unsigned long i;
  char          s[1024];
  unsigned long nv;
  long*         d;
  long*         p;
//  K_RATPOLY     P;
  
  in_fs >> s;
  assert(!strcmp(s, "kratpoly"));
  in_fs >> s;
  assert(!strcmp(s, "{"));
  in_fs >> s;
  assert(!strcmp(s, "nvars"));
  in_fs >> nv;
  
  if (nv > 0)
  {
    d = new long [nv];
    
    in_fs >> s;
    assert(!strcmp(s, "degree"));
    
    for (i = 0; i < nv; i++)
      in_fs >> d[i];
  }
  else  //  if (nv == 0)
    d = 0;
  
//  P = K_RATPOLY(nv, d);
  K_RATPOLY P(nv, d);
  
  if (d)
    delete [] d;
  
  if (nv > 0)
    p = new long [nv];
  
  in_fs >> s;
  
  while (!strcmp(s, "coeff"))
  {
    for (i = 0; i < nv; i++)
      in_fs >> p[i];
    
    in_fs >> P.get_coeff(p);
    in_fs >> s;
  }
  
  if (p)
    delete [] p;
  
  P.reduce_deg();
  
  return P;
}

