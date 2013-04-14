#include <kfloatpoly.h>

K_FLOATPOLY :: K_FLOATPOLY()
{
  num_vars = 0;
  deg      = 0;
  
  coeffs    = new double [num_coeffs = 1];
  coeffs[0] = 0.0;
}

K_FLOATPOLY :: K_FLOATPOLY(const unsigned long n, const long* const d)
{
  unsigned long i;
  
  if ((num_vars = n) > 0)
    deg = new long [num_vars];
  else
    deg = 0;
  
  num_coeffs = 1;
  
  for (i = 0; i < num_vars; i++)
  {
    deg[i]      = d[i];
    num_coeffs *= d[i] + 1;
  }
  
  coeffs = new double [num_coeffs];
  
  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = 0.0;
}

K_FLOATPOLY :: K_FLOATPOLY(const unsigned long n, const long* const d,
                           const double* const x)
{
  unsigned long i;
  
  if ((num_vars = n) > 0)
    deg = new long [num_vars];
  else
    deg = 0;
  
  num_coeffs = 1;
  
  for (i = 0; i < num_vars; i++)
  {
    deg[i]      = d[i];
    num_coeffs *= d[i] + 1;
  }
  
  coeffs = new double [num_coeffs];
  
  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = x[i];
}

K_FLOATPOLY :: K_FLOATPOLY(const K_FLOATPOLY& X)
{
  unsigned long i;
  
  if ((num_vars = X.num_vars) > 0)
  {
    deg = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      deg[i] = X.deg[i];
  }
  else
    deg = 0;
  
  coeffs = new double [num_coeffs = X.num_coeffs];
  
  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = X.coeffs[i];
}

K_FLOATPOLY& K_FLOATPOLY :: operator =(const K_FLOATPOLY& X)
{
  if (this != &X)
  {
    unsigned long i;
    
    if (deg)
      delete [] deg;
    
    if (coeffs)
      delete [] coeffs;
    
    if ((num_vars = X.num_vars) > 0)
    {
      deg = new long [num_vars];
      
      for (i = 0; i < num_vars; i++)
        deg[i] = X.deg[i];
    }
    else
      deg = 0;
    
    coeffs = new double [num_coeffs = X.num_coeffs];
    
    for (i = 0; i < num_coeffs; i++)
      coeffs[i] = X.coeffs[i];
  }
  
  return *this;
}

K_FLOATPOLY :: ~K_FLOATPOLY()
{
  if (deg)
    delete [] deg;
  
  if (coeffs)
    delete [] coeffs;
}

ostream& K_FLOATPOLY :: output(ostream& o) const
{
  unsigned long i, j;
  long*         p;
  
  for (i = 0; i < num_coeffs; i++)
//    if (coeffs[i] != 0.0)
      if (p = index_to_powers(i))
      {
        o << coeffs[i] << " X^(";
        
        for (j = 0; j < num_vars - 1; j++)
          o << p[j] << ',';
        
        o << p[num_vars - 1] << ')' << endl;
        delete [] p;
      }
      else
        o << coeffs[i] << endl;
  
  o  << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const K_FLOATPOLY& X)
{
  return X.output(o);
}

long* K_FLOATPOLY :: index_to_powers(const unsigned long i) const
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
  else
    p = 0;
  
  return p;
}

unsigned long K_FLOATPOLY :: powers_to_index(const long* const p) const
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

int K_FLOATPOLY :: reduce_deg()
{
  unsigned long i, j, k, l;
  long*         d;
  long*         p;
  unsigned long m;
  double*       c;
  int           reduced;
  
  reduced = 0;
  
  if (num_vars > 0)
  {
    d = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      d[i] = 0;
    
    for (i = 0; i < num_coeffs; i++)
      if (coeffs[i] != 0.0 && (p = index_to_powers(i)))
      {
        for (j = 0; j < num_vars; j++)
          if (d[j] < p[j])
            d[j] = p[j];
        
        delete [] p;
      }
    
    i = 0;
    
    while (!reduced && i < num_vars)
    {
      if (d[i] < deg[i])
        reduced = 1;
      
      i++;
    }
    
    if (reduced)
    {
      for (i = 0; i < num_vars; i++)
        m *= d[i] + 1;
      
      c = new double [m];
      
      for (i = 0; i < num_coeffs; i++)
        if (coeffs[i] != 0.0)
        {
          p = index_to_powers(i);
          
          for (j = k = 0, l = m; k < num_vars; k++)
          {
            assert(p[k] <= d[k]);
            
            l /= d[k] + 1;
            j += l * (d[k] - p[k]);
          }
          
          c[j] = coeffs[i];
          delete [] p;  //  num_vars > 0 => p != 0
        }
      
      delete [] coeffs;
      num_coeffs = m;
      coeffs     = c;
    }
    
    delete [] deg;  //  num_vars > 0 => deg != 0
    deg = d;
  }
  
  return reduced;
}

K_FLOATPOLY K_FLOATPOLY :: add(const K_FLOATPOLY& X) const
{
  assert(num_vars == X.num_vars);
  
  unsigned long i, j;
  long*         dy;
  long*         py;
  
  if (num_vars > 0)
  {
    dy = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dy[i] = deg[i] > X.deg[i] ? deg[i] : X.deg[i];
  }
  else
    dy = 0;
  
  K_FLOATPOLY Y(num_vars, dy);
  
  if (dy)
    delete [] dy;
  
  for (i = 0; i < num_coeffs; i++)
  {
    py          = index_to_powers(i);
    j           = Y.powers_to_index(py);
    Y.coeffs[j] = coeffs[i];
    
    if (py)
      delete [] py;
  }
  
  for (i = 0; i < X.num_coeffs; i++)
  {
    py           = X.index_to_powers(i);
    j            = Y.powers_to_index(py);
    Y.coeffs[j] += X.coeffs[i];
    
    if (py)
      delete [] py;
  }
  
  Y.reduce_deg();
  
  return Y;
}

K_FLOATPOLY operator +(const K_FLOATPOLY& X, const K_FLOATPOLY& Y)
{
  return X.add(Y);
}

K_FLOATPOLY K_FLOATPOLY :: sub(const K_FLOATPOLY& X) const
{
  assert(num_vars == X.num_vars);
  
  unsigned long i, j;
  long*         dy;
  long*         py;
  
  if (num_vars > 0)
  {
    dy = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dy[i] = deg[i] > X.deg[i] ? deg[i] : X.deg[i];
  }
  else
    dy = 0;
  
  K_FLOATPOLY Y(num_vars, dy);
  
  if (dy)
    delete [] dy;
  
  for (i = 0; i < num_coeffs; i++)
  {
    py          = index_to_powers(i);
    j           = Y.powers_to_index(py);
    Y.coeffs[j] = coeffs[i];
    
    if (py)
      delete [] py;
  }
  
  for (i = 0; i < X.num_coeffs; i++)
  {
    py           = X.index_to_powers(i);
    j            = Y.powers_to_index(py);
    Y.coeffs[j] -= X.coeffs[i];
    
    if (py)
      delete [] py;
  }
  
  Y.reduce_deg();
  
  return Y;
}

K_FLOATPOLY operator -(const K_FLOATPOLY& X, const K_FLOATPOLY& Y)
{
  return X.sub(Y);
}

K_FLOATPOLY K_FLOATPOLY :: mul(const K_FLOATPOLY& X) const
{
  assert(num_vars == X.num_vars);
  
  unsigned long i, j, k;
  long*         dy;
  long*         p;
  long*         px;
  long*         py;
  
  if (num_vars > 0)
  {
    dy = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dy[i] = deg[i] + X.deg[i];
  }
  else
    dy = 0;
  
  K_FLOATPOLY Y(num_vars, dy);
  
  for (i = 0; i < num_coeffs; i++)
    if (coeffs[i] != 0.0)
    {
      p = index_to_powers(i);
      
      for (j = 0; j < X.num_coeffs; j++)
        if (X.coeffs[j] != 0.0)
        {
          px = X.index_to_powers(j);
          
          if (num_vars > 0)
          {
            py = new long [num_vars];
            
            for (k = 0; k < num_vars; k++)
              py[k] = p[k] + px[k];
          }
          else
            py = 0;
          
          k            = Y.powers_to_index(py);
          Y.coeffs[k] += coeffs[i] * X.coeffs[j];
          
          if (py)
            delete [] py;
          
          if (px)
            delete [] px;
        }
      
      if (p)
        delete [] p;
    }
  
  Y.reduce_deg();
  
  return Y;
}

K_FLOATPOLY operator *(const K_FLOATPOLY& X, const K_FLOATPOLY& Y)
{
  return X.mul(Y);
}

K_FLOATPOLY K_FLOATPOLY :: mul(const double x) const
{
  unsigned long i;
  long*         dy;
  
  if (num_vars > 0)
  {
    dy = new long [num_vars];
    
    for (i = 0; i < num_vars; i++)
      dy[i] = x != 0 ? deg[i] : 0;
  }
  else
    dy = 0;
  
  K_FLOATPOLY Y(num_vars, dy);
  
  if (dy)
    delete [] dy;
  
  if (x != 0)
    for (i = 0; i < num_coeffs; i++)
      Y.coeffs[i] = coeffs[i] * x;
//  Assume the constructor for Y initializes all its coefficients to be 0.0.
//  else
//    for (i = 0; i < Y.num_coeffs; i++)
//      Y.coeffs[i] = 0.0;
  
  return Y;
}

K_FLOATPOLY operator *(const K_FLOATPOLY& X, const double y)
{
  return X.mul(y);
}

K_FLOATPOLY operator *(const double x, const K_FLOATPOLY& Y)
{
  return Y.mul(x);
}

K_FLOATPOLY K_FLOATPOLY :: neg() const
{
  unsigned long i;
  
  K_FLOATPOLY X(num_vars, deg);
  
  for (i = 0; i < num_coeffs; i++)
    X.coeffs[i] = - coeffs[i];
  
  return X;
}

K_FLOATPOLY operator -(const K_FLOATPOLY& X)
{
  return X.neg();
}

K_FLOATPOLY K_FLOATPOLY :: derivative(const unsigned long i) const
{
  assert(!num_vars || i < num_vars);
  
  unsigned long j, k;
  unsigned long c;
  long*         dx;
  long*         p;
  
  if (num_vars > 0)
  {
    dx = new long [num_vars];
    
    if (deg[i] > 0)
    {
      for (j = 0; j < i; j++)
        dx[j] = deg[j];
      
      dx[i] = deg[i] - 1;
      
      for (j = i + 1; j < num_vars; j++)
        dx[j] = deg[j];
    }
    else
      for (j = 0; j < num_vars; j++)
        dx[j] = 0;
  }
  else
    dx = 0;
  
  K_FLOATPOLY X(num_vars, dx);
  
  if (dx)
    delete [] dx;
  
  for (j = 0; j < num_coeffs; j++)
    if (coeffs[j] != 0.0 && (p = index_to_powers(j)))
    {
      if ((c = p[i]) > 0)
      {
        p[i]--;
        k           = X.powers_to_index(p);
        X.coeffs[k] = coeffs[j] * c;
      }
      
      delete [] p;
    }
  
  return X;
}

K_FLOATPOLY K_FLOATPOLY :: subst_first_var(const double x) const
{
  unsigned long i, j;
  long*         dy;
  
  if (num_vars > 0)
  {
    if (num_vars > 1)
    {
      dy = new long [num_vars - 1];
      
      for (i = 1; i < num_vars; i++)
        dy[i - 1] = deg[i];
    }
    else
      dy = 0;
    
    K_FLOATPOLY Y(num_vars - 1, dy);
    
    if (dy)
      delete [] dy;
    
    for (i = 0; i < Y.num_coeffs; i++)
      Y.coeffs[i] = coeffs[i];
    
    for (i = 1; i <= deg[0]; i++)
      for (j = 0; j < Y.num_coeffs; j++)
      {
        Y.coeffs[j] *= x;
        Y.coeffs[j] += coeffs[j + i * Y.num_coeffs];
      }
    
    Y.reduce_deg();
    
    return Y;
  }
  else
    return K_FLOATPOLY(*this);
}

double K_FLOATPOLY :: evaluate(const double* const x) const
{
  unsigned long i;
  
  K_FLOATPOLY Y(*this);
  
  for (i = 0; i < num_vars; i++)
    Y = Y.subst_first_var(x[i]);
  
  return Y.coeffs[0];
}

K_FLOATPOLY K_FLOATPOLY :: subst_val(const unsigned long i,
                                     const double x) const
{
  assert(!num_vars || i < num_vars);
  
  unsigned long j, k;
  long*         dy;
  double*       xp;
  long*         p;
  long*         py;
  
  if (num_vars > 0)
  {
    if (!i)
      return subst_first_var(x);
    else
    {
      if (num_vars > 1)
      {
        dy = new long [num_vars - 1];
        
        for (j = 0; j < i; j++)
          dy[j] = deg[j];
        
        for (j = i + 1; j < num_vars; j++)
          dy[j - 1] = deg[j];
      }
      else  //  num_vars == 1
        dy = 0;
      
      K_FLOATPOLY Y(num_vars - 1, dy);
      
      if (dy)
        delete [] dy;
      
      xp    = new double [deg[i] + 1];
      xp[0] = 1;
      
      if (deg[i] > 0)
      {
        xp[1] = x;
        
        for (j = 1; j <= deg[i]; j++)
          xp[j] = x * xp[j - 1];
      }
      
      for (j = 0; j < num_coeffs; j++)
        if (coeffs[j] != 0.0)
        {
          p = index_to_powers(j);
          
          if (num_vars > 1)
          {
            py = new long [num_vars - 1];
            
            for (k = 0; k < i; k++)
              py[k] = p[k];
            
            for (k = i + 1; k < num_vars; k++)
              py[k - 1] = p[k];
          }
          else
            py = 0;
          
          k            = Y.powers_to_index(py);
          Y.coeffs[k] += coeffs[j] * xp[p[i]];
          
          if (py)
            delete [] py;
          
          delete [] p;  //  num_vars > 0 => p != 0
        }
      
      delete [] xp;
      
      Y.reduce_deg();
      
      return Y;
    }
  }
  else
    return K_FLOATPOLY(*this);
}

K_FLOATPOLY K_FLOATPOLY :: subst_expr(const unsigned long i,
                                      const K_FLOATPOLY& X) const
{
  assert(i < num_vars);                //  i < num_vars => num_vars > 0
  assert(num_vars == X.num_vars + 1);
  
  unsigned long j, k, l;
  long*         dy;
  K_FLOATPOLY*  XP;
  long*         dxp0;
  long*         p;
  long*         p1;
  long*         pxpe;
  long*         py;
  
  if (num_vars > 1)
  {
    dy = new long [num_vars - 1];
    
    for (j = 0; j < i; j++)
      dy[j] = deg[j] + deg[i] * X.deg[j];
    
    for (j = i + 1; j < num_vars; j++)
      dy[j - 1] = deg[j] + deg[i] * X.deg[j - 1];
  }
  else  //  num_vars == 1
    dy = 0;
  
  K_FLOATPOLY Y(num_vars - 1, dy);
  
  if (dy)
    delete [] dy;
  
  XP = new K_FLOATPOLY [deg[i] + 1];
  
  if (X.num_vars > 0)
  {
    dxp0 = new long [X.num_vars];
    
    for (j = 0; j < X.num_vars; j++)
      dxp0[j] = 0;
  }
  else
    dxp0 = 0;
  
  XP[0] = K_FLOATPOLY(X.num_vars, dxp0);
  
  if (dxp0)
    delete [] dxp0;
  
  XP[0].coeffs[0] = 1.0;
  
  if (deg[i] > 0)
  {
    XP[1] = X;
    
    for (j = 2; j <= deg[i]; j++)
      XP[j] = X * XP[j - 1];
  }
  
  for (j = 0; j < num_coeffs; j++)
    if (coeffs[j] != 0.0)
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
      else  //  num_vars == 1
        p1 = 0;
      
      for (k = 0; k < XP[p[i]].num_coeffs; k++)
        if (XP[p[i]].coeffs[k] != 0.0)
        {
          pxpe = XP[p[i]].index_to_powers(k);
          
          if (num_vars > 1)
          {
            py = new long [num_vars - 1];
            
            for (l = 0; l < num_vars - 1; l++)
              py[l] = p1[l] + pxpe[l];
          }
          else  //  num_vars == 1
            py = 0;
          
          l            = Y.powers_to_index(py);
          Y.coeffs[l] += coeffs[j] * XP[p[i]].coeffs[k];
          
          if (py)
            delete [] py;
          
          if (pxpe)
            delete [] pxpe;
        }
      
      if (p1)
        delete [] p1;
      
      delete [] p;  //  num_vars > 0 => p != 0
    }
  
  delete [] XP;
  Y.reduce_deg();
  
  return Y;
}

K_FLOATPOLY K_FLOATPOLY :: subst_expr(const unsigned long i,
                                      const K_FLOATPOLY& N,
                                      const K_FLOATPOLY& D) const
{
  assert(i < num_vars);                //  i < num_vars => num_vars > 0
  assert(num_vars == N.num_vars + 1);
  assert(num_vars == D.num_vars + 1);
  
  unsigned long j, k, l, m;
  long*         dy;
  K_FLOATPOLY*  NP;
  K_FLOATPOLY*  DP;
  long*         dnp0;
  long*         ddp0;
  long*         p;
  long*         p1;
  long*         pnpe;
  long*         pdpf;
  long*         py;
  
  if (num_vars > 1)
  {
    dy = new long [num_vars - 1];
    
    for (j = 0; j < i; j++)
      if (N.deg[j] > D.deg[j])
        dy[j] = deg[j] + deg[i] * N.deg[j];
      else
        dy[j] = deg[j] + deg[i] * D.deg[j];
    
    for (j = i + 1; j < num_vars; j++)
      if (N.deg[j - 1] > D.deg[j - 1])
        dy[j - 1] = deg[j] + deg[i] * N.deg[j - 1];
      else
        dy[j - 1] = deg[j] + deg[i] * D.deg[j - 1];
  }
  else
    dy = 0;
  
  K_FLOATPOLY Y(num_vars - 1, dy);
  
  if (dy)
    delete [] dy;
  
  NP = new K_FLOATPOLY [deg[i] + 1];
  
  if (N.num_vars > 0)
  {
    dnp0 = new long [N.num_vars];
    
    for (j = 0; j < N.num_vars; j++)
      dnp0[j] = 0;
  }
  else
    dnp0 = 0;
  
  NP[0] = K_FLOATPOLY(N.num_vars, dnp0);
  
  if (dnp0)
    delete [] dnp0;
  
  NP[0].coeffs[0] = 1.0;
  
  if (deg[i] > 0)
  {
    NP[1] = N;
    
    for (j = 2; j <= deg[i]; j++)
      NP[j] = N * NP[j - 1];
  }
  
  DP = new K_FLOATPOLY [deg[i] + 1];
  
  if (D.num_vars > 0)
  {
    ddp0 = new long [D.num_vars];
    
    for (j = 0; j < D.num_vars; j++)
      ddp0[j] = 0;
  }
  else
    ddp0 = 0;
  
  DP[0] = K_FLOATPOLY(D.num_vars, ddp0);
  
  if (ddp0)
    delete [] ddp0;
  
  DP[0].coeffs[0] = 1.0;
  
  if (deg[i] > 0)
  {
    DP[1] = D;
    
    for (j = 2; j <= deg[i]; j++)
      DP[j] = D * DP[j - 1];
  }
  
  for (j = 0; j < num_coeffs; j++)
    if (coeffs[j] != 0.0)
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
        if (NP[p[i]].coeffs[k] != 0.0)
        {
          pnpe = NP[p[i]].index_to_powers(k);
          
          for (l = 0; l < DP[deg[i] - p[i]].num_coeffs; l++)
            if (DP[deg[i] - p[i]].coeffs[l] != 0.0)
            {
              pdpf = DP[deg[i] - p[i]].index_to_powers(l);
              
              if (num_vars > 1)
              {
                py = new long [num_vars - 1];
                
                for (m = 0; m < num_vars - 1; m++)
                  py[m] = p1[m] + pnpe[m] + pdpf[m];
              }
              else  //  num_vars == 1
                py = 0;
              
              m            = Y.powers_to_index(py);
              Y.coeffs[m] += coeffs[j] * NP[p[i]].coeffs[k]
                                       * DP[deg[i] - p[i]].coeffs[l];
              
              if (py)
                delete [] py;
              
              if (pdpf)
                delete [] pdpf;
            }
          
          if (pnpe)
            delete [] pnpe;
        }
      
      if (p1)
        delete [] p1;
      
      delete [] p;  //  num_vars > 0 => p!= 0
    }
  
  delete [] NP;
  delete [] DP;
  Y.reduce_deg();
  
  return Y;
}

extern "C" {
int dgeev_(char*, char*, long*, double*,
           long*,	double*, double*, double*,
           long*, double*, long*, double*,
           long*, long*);
}

//  unsigned long K_FLOATPOLY :: gen_fp_roots(const double l,
//                                            const doubld h,
//                                            double*& X) const
//    PROVIDED *this is a univariate polynomial,
//    compute roots X of *this in [l, h] by using only FP arithmetic.

unsigned long K_FLOATPOLY :: gen_fp_roots(const double l,
                                          const double h,
                                          double*& X) const
{
  assert(num_vars == 1);
  assert(l <= h);
  
  unsigned long i, j;
  long          d;
  unsigned long n;
  double*       Y;
  
  n = 0;
  
  if ((d = deg[0]) > 0)
  {
    double* C;
    
    Y = new double [d];
    C = new double [d * d];
    
    for (i = 0; i < d - 1; i++)
    {
      for (j = 0; j < i + 1; j++)
        C[j * d + i] = 0.0;
      
      C[(i + 1) * d + i] = 1.0;
      
      for (j = i + 2; j < d; j++)
        C[j * d + i] = 0.0;
    }
    
    j = 0;
    
    while (j < d && finite(C[j * d + d - 1] = - coeffs[d - j] / coeffs[0]))
      j++;
    
    if (j == d)
    {
      double* eigenval_re;
      double* eigenval_im;
      double* dummy_l;
      long    l_dummy_l;
      double* dummy_r;
      long    l_dummy_r;
      double* work;
      long    l_work;
      long    OK;
      
      eigenval_re = new double [d];
      eigenval_im = new double [d];
      
      dummy_l = new double [l_dummy_l = 1];
      dummy_r = new double [l_dummy_r = 1];
      
      l_work = 100 * 4 * d * d + 1000;  // We allocate a space for work
      work   = new double [l_work];     // 100 times larger than necessary.
      
      dgeev_("N", "N", &d, C, &d,
             eigenval_re, eigenval_im,
             dummy_l, &l_dummy_l, dummy_r, &l_dummy_r,
             work, &l_work,
             &OK);
      
      delete [] dummy_l;
      delete [] dummy_r;
      
      delete [] work;
      
      if (!OK)  //  OK == 0 iff dgeev_() was sucessful.
        for (i = 0; i < d; i++)
          if (eigenval_im[i] == 0.0)
            if (eigenval_re[i] >= l && eigenval_re[i] <= h)
              Y[n++] = eigenval_re[i];
      
      delete [] eigenval_re;
      delete [] eigenval_im;
    }
    
    delete [] C;
  }
  else  //  if (d == 0)
    Y = 0;
  
  if (n > 0)
  {
    X = new double [n];
    
    for (i = 0; i < n; i++)
      X[i] = Y[i];
  }
  else  //  if (n == 0)
    X = 0;
  
  if (Y)
    delete [] Y;
  
  return n;
}

