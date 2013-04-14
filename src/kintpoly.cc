#include <kintpoly.h>

K_INTPOLY :: K_INTPOLY()
{
  num_vars = 0;
  deg      = 0;

  coeffs    = new bigint [num_coeffs = 1];
  coeffs[0] = 0;
}

K_INTPOLY :: K_INTPOLY(const unsigned long n, const long* const d)
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

  coeffs = new bigint [num_coeffs];

  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = 0;
}

K_INTPOLY :: K_INTPOLY(const unsigned long n, const long* const d,
                       const bigint* const x)
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

  coeffs = new bigint [num_coeffs];

  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = x[i];
}

K_INTPOLY :: K_INTPOLY(const K_INTPOLY& X)
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

  coeffs = new bigint [num_coeffs = X.num_coeffs];

  for (i = 0; i < num_coeffs; i++)
    coeffs[i] = X.coeffs[i];
}

K_INTPOLY& K_INTPOLY :: operator =(const K_INTPOLY& X)
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

    coeffs = new bigint [num_coeffs = X.num_coeffs];

    for (i = 0; i < num_coeffs; i++)
      coeffs[i] = X.coeffs[i];
  }

  return *this;
}

K_INTPOLY :: ~K_INTPOLY()
{
  if (deg)
    delete [] deg;

  if (coeffs)
    delete [] coeffs;
}

ostream& K_INTPOLY :: output(ostream& o) const
{
  unsigned long i, j;
  long*         p;

  for (i = 0; i < num_coeffs; i++)
//    if (sgn(coeffs[i]))
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

ostream& operator <<(ostream& o, const K_INTPOLY& X)
{
  return X.output(o);
}

long* K_INTPOLY :: index_to_powers(const unsigned long i) const
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

unsigned long K_INTPOLY :: powers_to_index(const long* const p) const
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

int K_INTPOLY :: reduce_deg()
{
  unsigned long i, j, k, l;
  long*         d;
  long*         p;
  unsigned long m;
  bigint*       c;
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
    {
      if (d[i] < deg[i])
        reduced = 1;

      i++;
    }

    if (reduced)
    {
      for (i = 0; i < num_vars; i++)
        m *= d[i] + 1;

      c = new bigint [m];

      for (i = 0; i < num_coeffs; i++)
        if (sgn(coeffs[i]))
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

K_INTPOLY K_INTPOLY :: add(const K_INTPOLY& X) const
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

  K_INTPOLY Y(num_vars, dy);

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

K_INTPOLY operator +(const K_INTPOLY& X, const K_INTPOLY& Y)
{
  return X.add(Y);
}

K_INTPOLY K_INTPOLY :: sub(const K_INTPOLY& X) const
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

  K_INTPOLY Y(num_vars, dy);

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

K_INTPOLY operator -(const K_INTPOLY& X, const K_INTPOLY& Y)
{
  return X.sub(Y);
}

K_INTPOLY K_INTPOLY :: mul(const K_INTPOLY& X) const
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

  K_INTPOLY Y(num_vars, dy);

  for (i = 0; i < num_coeffs; i++)
    if (sgn(coeffs[i]))
    {
      p = index_to_powers(i);

      for (j = 0; j < X.num_coeffs; j++)
        if (sgn(X.coeffs[j]))
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

K_INTPOLY operator *(const K_INTPOLY& X, const K_INTPOLY& Y)
{
  return X.mul(Y);
}

K_INTPOLY K_INTPOLY :: mul(const bigint& x) const
{
  unsigned long i;
  long*         dy;

  if (num_vars > 0)
  {
    dy = new long [num_vars];

    for (i = 0; i < num_vars; i++)
      dy[i] = sgn(x) ? deg[i] : 0;
  }
  else
    dy = 0;

  K_INTPOLY Y(num_vars, dy);

  if (dy)
    delete [] dy;

  if (sgn(x))
    for (i = 0; i < num_coeffs; i++)
      Y.coeffs[i] = coeffs[i] * x;
//  Assume the constructor for Y initializes all its coefficients to be 0.
//  else
//    for (i = 0; i < Y.num_coeffs; i++)
//      Y.coeffs[i] = 0;

  return Y;
}

K_INTPOLY operator *(const K_INTPOLY& X, const bigint& y)
{
  return X.mul(y);
}

K_INTPOLY operator *(const bigint& x, const K_INTPOLY& Y)
{
  return Y.mul(x);
}

K_INTPOLY K_INTPOLY :: neg() const
{
  unsigned long i;

  K_INTPOLY X(num_vars, deg);

  for (i = 0; i < num_coeffs; i++)
    X.coeffs[i] = - coeffs[i];

  return X;
}

K_INTPOLY operator -(const K_INTPOLY& X)
{
  return X.neg();
}

K_INTPOLY K_INTPOLY :: derivative(const unsigned long i) const
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

  K_INTPOLY X(num_vars, dx);

  if (dx)
    delete [] dx;

  for (j = 0; j < num_coeffs; j++)
    if (sgn(coeffs[j]) && (p = index_to_powers(j)))
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

K_INTPOLY K_INTPOLY :: subst_first_var(const bigint& x) const
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

    K_INTPOLY Y(num_vars - 1, dy);

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
    return K_INTPOLY(*this);
}

bigint K_INTPOLY :: evaluate(const bigint* const x) const
{
  unsigned long i;

  K_INTPOLY Y(*this);

  for (i = 0; i < num_vars; i++)
    Y = Y.subst_first_var(x[i]);

  return Y.coeffs[0];
}

int K_INTPOLY :: sgn_at(const bigint* const x) const
{
  if (num_vars == 1)
    return sgn_at(x[0]);
  else
    return sgn(evaluate(x));
}

K_INTPOLY K_INTPOLY :: subst_val(const unsigned long i, const bigint& x) const
{
  assert(!num_vars || i < num_vars);

  unsigned long j, k;
  long*         dy;
  bigint*       xp;
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

      K_INTPOLY Y(num_vars - 1, dy);

      if (dy)
        delete [] dy;

      xp    = new bigint [deg[i] + 1];
      xp[0] = 1;

      if (deg[i] > 0)
      {
        xp[1] = x;

        for (j = 1; j <= deg[i]; j++)
          xp[j] = x * xp[j - 1];
      }

      for (j = 0; j < num_coeffs; j++)
        if (sgn(coeffs[j]))
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
    return K_INTPOLY(*this);
}

K_INTPOLY K_INTPOLY :: subst_expr(const unsigned long i,
                                  const K_INTPOLY& X) const
{
  assert(i < num_vars);                //  i < num_vars => num_vars > 0
  assert(num_vars == X.num_vars + 1);

  unsigned long j, k, l;
  long*         dy;
  K_INTPOLY*    XP;
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

  K_INTPOLY Y(num_vars - 1, dy);

  if (dy)
    delete [] dy;

  XP = new K_INTPOLY [deg[i] + 1];

  if (X.num_vars > 0)
  {
    dxp0 = new long [X.num_vars];

    for (j = 0; j < X.num_vars; j++)
      dxp0[j] = 0;
  }
  else
    dxp0 = 0;

  XP[0] = K_INTPOLY(X.num_vars, dxp0);

  if (dxp0)
    delete [] dxp0;

  XP[0].coeffs[0] = 1;

  if (deg[i] > 0)
  {
    XP[1] = X;

    for (j = 2; j <= deg[i]; j++)
      XP[j] = X * XP[j - 1];
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
      else  //  num_vars == 1
        p1 = 0;

      for (k = 0; k < XP[p[i]].num_coeffs; k++)
        if (sgn(XP[p[i]].coeffs[k]))
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

K_INTPOLY K_INTPOLY :: subst_expr(const unsigned long i,
                                  const K_INTPOLY& N,
                                  const K_INTPOLY& D) const
{
  assert(i < num_vars);                //  i < num_vars => num_vars > 0
  assert(num_vars == N.num_vars + 1);
  assert(num_vars == D.num_vars + 1);

  unsigned long j, k, l, m;
  long*         dy;
  K_INTPOLY*    NP;
  K_INTPOLY*    DP;
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

  K_INTPOLY Y(num_vars - 1, dy);

  if (dy)
    delete [] dy;

  NP = new K_INTPOLY [deg[i] + 1];

  if (N.num_vars > 0)
  {
    dnp0 = new long [N.num_vars];

    for (j = 0; j < N.num_vars; j++)
      dnp0[j] = 0;
  }
  else
    dnp0 = 0;

  NP[0] = K_INTPOLY(N.num_vars, dnp0);

  if (dnp0)
    delete [] dnp0;

  NP[0].coeffs[0] = 1;

  if (deg[i] > 0)
  {
    NP[1] = N;

    for (j = 2; j <= deg[i]; j++)
      NP[j] = N * NP[j - 1];
  }

  DP = new K_INTPOLY [deg[i] + 1];

  if (D.num_vars > 0)
  {
    ddp0 = new long [D.num_vars];

    for (j = 0; j < D.num_vars; j++)
      ddp0[j] = 0;
  }
  else
    ddp0 = 0;

  DP[0] = K_INTPOLY(D.num_vars, ddp0);

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
          pnpe = NP[p[i]].index_to_powers(k);

          for (l = 0; l < DP[deg[i] - p[i]].num_coeffs; l++)
            if (sgn(DP[deg[i] - p[i]].coeffs[l]))
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

      delete [] p;  //  num_vars > 0 => p != 0
    }

  delete [] NP;
  delete [] DP;
  Y.reduce_deg();

  return Y;
}

int K_INTPOLY :: sgn_at(const bigint& x) const
{
  assert(num_vars == 1);

  return sgn(subst_first_var(x).coeffs[0]);
}

