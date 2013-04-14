#include <bigrational_vector.h>

bigrational_vector :: bigrational_vector(const unsigned long n)
  : dim(n)
{
  if (dim > 0)
    rep = new bigrational [dim];
  else
    rep = 0;
}

bigrational_vector :: bigrational_vector(const unsigned long n,
                                         const bigrational* X)
  : dim(n)
{
  unsigned long i;

  if (dim > 0)
  {
    rep = new bigrational [dim];

    for (i = 0; i < dim; i++)
      rep[i] = X[i];
  }
  else
    rep = 0;
}

//bigrational_vector :: bigrational_vector(const bigint_vector& X)
//  : dim(X.dim)
//{
//  unsigned long i;
//
//  if (dim > 0)
//  {
//    rep = new bigrational [dim];
//
//    for (i = 0; i < dim; i++)
//      rep[i] = X[i];
//  }
//  else
//    rep = 0;
//}

bigrational_vector :: bigrational_vector(const bigrational_vector& X)
  : dim(X.dim)
{
  unsigned long i;

  if (dim > 0)
  {
    rep = new bigrational [dim];

    for (i = 0; i < dim; i++)
      rep[i] = X[i];
  }
  else
    rep = 0;
}

bigrational_vector&
bigrational_vector :: operator =(const bigrational_vector& X)
{
  unsigned long i;

  if (this != &X)
  {
    if (dim != X.dim)
    {
      if (rep)
        delete [] rep;

      if (dim = X.dim)
        rep = new bigrational [dim];
      else
        rep = 0;
    }

    for (i = 0; i < dim; i++)
      rep[i] = X[i];
  }

  return *this;
}

bigrational_vector :: ~bigrational_vector()
{
  if (rep)
    delete [] rep;
}

bigrational_vector bigrational_vector :: add(const bigrational_vector& X) const
{
  assert(dim == X.dim);

  unsigned long      i;
  bigrational_vector Y(dim);

  for (i = 0; i < dim; i++)
    Y[i] = this->operator [](i) + X[i];

  return Y;
}

bigrational_vector operator +(const bigrational_vector& X,
                              const bigrational_vector& Y)
{
  return X.add(Y);
}

bigrational_vector bigrational_vector :: sub(const bigrational_vector& X) const
{
  assert(dim == X.dim);

  unsigned long      i;
  bigrational_vector Y(dim);

  for (i = 0; i < dim; i++)
    Y[i] = this->operator [](i) - X[i];

  return Y;
}

bigrational_vector operator -(const bigrational_vector& X,
                              const bigrational_vector& Y)
{
  return X.sub(Y);
}

bigrational_vector bigrational_vector :: neg() const
{
  unsigned long      i;
  bigrational_vector X(dim);

  for (i = 0; i < dim; i++)
    X[i] = - this->operator [](i);

  return X;
}

bigrational_vector operator -(const bigrational_vector& X)
{
  return X.neg();
}

bigrational bigrational_vector :: in_prod(const bigrational_vector& X) const
{
  assert(dim == X.dim);

  unsigned long i;
  bigrational   y;

  y = 0;

  for (i = 0; i < dim; i++)
    y += this->operator [](i) * X[i];

  return y;
}

bigrational in_prod(const bigrational_vector& X, const bigrational_vector& Y)
{
  return X.in_prod(Y);
}

bigrational_vector bigrational_vector :: scalar_mul(const bigrational& x) const
{
  unsigned long      i;
  bigrational_vector X(dim);

  for (i = 0; i < dim; i++)
    X[i] = this->operator [](i) * x;

  return X;
}

bigrational_vector scalar_mul(const bigrational& x,
                              const bigrational_vector& Y)
{
  return Y.scalar_mul(x);
}

int bigrational_vector :: cmp(const bigrational_vector& X) const
{
  assert(dim == X.dim);

  unsigned long i;
  int           c;

  i = 0;

  while (i < dim && this->operator[](i) ==  X[i])
    i++;

  if (i == dim)
    c = 0;
  else  //  if (i < dim)
    c = 1;

  return c;
}

int operator ==(const bigrational_vector& X, const bigrational_vector& Y)
{
  return X.cmp(Y) == 0;
}

int operator !=(const bigrational_vector& X, const bigrational_vector& Y)
{
  return X.cmp(Y) != 0;
}

int bigrational_vector :: is_zero() const
{
  unsigned long i;

  for (i = 0; i < dim && !sgn(this->operator[](i)); i++)
    ;

  return i == dim;
}

int is_zero(const bigrational_vector& X)
{
  return X.is_zero();
}

bigrational_vector&
bigrational_vector :: operator +=(const bigrational_vector& X)
{
  assert(dim == X.dim);

  unsigned long i;

  for (i = 0; i < dim; i++)
    this->operator[](i) += X[i];

  return *this;
}

bigrational_vector&
bigrational_vector :: operator -=(const bigrational_vector& X)
{
  assert(dim == X.dim);

  unsigned long i;

  for (i = 0; i < dim; i++)
    this->operator[](i) -= X[i];

  return *this;
}

ostream& operator <<(ostream& o, const bigrational_vector& X)
{
  unsigned long i;

  if (X.dim > 0)
  {
    o << " ( ";

    for (i = 0; i < X.dim - 1; i++)
      o << X[i] << ", ";

    o << X[X.dim - 1] << " ) " << flush;
  }
  else  //  if (!X.dim)
    o << " NULL " << flush;

  return o;
}

