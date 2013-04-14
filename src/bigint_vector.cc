//  file:    bigint_vector.cc
//  update:  09/25/02

#include <bigint_vector.h>

bigint_vector :: bigint_vector(const unsigned long n)
  : dim(n)
{
  if (dim > 0)
    rep = new bigint [dim];
  else
    rep = 0;
}

bigint_vector :: bigint_vector(const unsigned long n, const bigint* X)
  : dim(n)
{
  unsigned long i;
  
  if (dim > 0)
  {
    rep = new bigint [dim];
    
    for (i = 0; i < dim; i++)
      rep[i] = X[i];
  }
  else
    rep = 0;
}

bigint_vector :: bigint_vector(const bigint_vector& X)
  : dim(X.dim)
{
  unsigned long i;
  
  if (dim > 0)
  {
    rep = new bigint [dim];
    
    for (i = 0; i < dim; i++)
      rep[i] = X[i];
  }
  else
    rep = 0;
}

bigint_vector& bigint_vector :: operator =(const bigint_vector& X)
{
  unsigned long i;
  
  if (this != &X)
  {
    if (dim != X.dim)
    {
      if (rep)
        delete [] rep;
      
      if (dim = X.dim)
        rep = new bigint [dim];
      else
        rep = 0;
    }
    
    for (i = 0; i < dim; i++)
      rep[i] = X[i];
  }
  
  return *this;
}

bigint_vector :: ~bigint_vector()
{
  if (rep)
    delete [] rep;
}

bigint_vector bigint_vector :: add(const bigint_vector& X) const
{
  assert(dim == X.dim);
  
  unsigned long i;
  bigint_vector Y(dim);
  
  for (i = 0; i < dim; i++)
    Y[i] = this->operator [](i) + X[i];
  
  return Y;
}

bigint_vector operator +(const bigint_vector& X, const bigint_vector& Y)
{
  return X.add(Y);
}

bigint_vector bigint_vector :: sub(const bigint_vector& X) const
{
  assert(dim == X.dim);
  
  unsigned long i;
  bigint_vector Y(dim);
  
  for (i = 0; i < dim; i++)
    Y[i] = this->operator [](i) - X[i];
  
  return Y;
}

bigint_vector operator -(const bigint_vector& X, const bigint_vector& Y)
{
  return X.sub(Y);
}

bigint_vector bigint_vector :: neg() const
{
  unsigned long i;
  bigint_vector X(dim);
  
  for (i = 0; i < dim; i++)
    X[i] = - this->operator [](i);
  
  return X;
}

bigint_vector operator -(const bigint_vector& X)
{
  return X.neg();
}

bigint bigint_vector :: in_prod(const bigint_vector& X) const
{
  assert(dim == X.dim);
  
  unsigned long i;
  bigint        y;
  
  y = 0;
  
  for (i = 0; i < dim; i++)
    y += this->operator [](i) * X[i];
  
  return y;
}

bigint in_prod(const bigint_vector& X, const bigint_vector& Y)
{
  return X.in_prod(Y);
}

bigint_vector bigint_vector :: scalar_mul(const bigint& x) const
{
  unsigned long i;
  bigint_vector X(dim);
  
  for (i = 0; i < dim; i++)
    X[i] = this->operator [](i) * x;
  
  return X;
}

bigint_vector scalar_mul(const bigint& x, const bigint_vector& Y)
{
  return Y.scalar_mul(x);
}

int bigint_vector :: cmp(const bigint_vector& X) const
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

int operator ==(const bigint_vector& X, const bigint_vector& Y)
{
  return X.cmp(Y) == 0;
}

int operator !=(const bigint_vector& X, const bigint_vector& Y)
{
  return X.cmp(Y) != 0;
}

int bigint_vector :: is_zero() const
{
  unsigned long i;
  
  for (i = 0; i < dim && !sgn(this->operator[](i)); i++)
    ;
  
  return i == dim;
}

int is_zero(const bigint_vector& X)
{
  return X.is_zero();
}

bigint_vector& bigint_vector :: operator +=(const bigint_vector& X)
{
  assert(dim == X.dim);
  
  unsigned long i;
  
  for (i = 0; i < dim; i++)
    this->operator[](i) += X[i];
  
  return *this;
}

bigint_vector& bigint_vector :: operator -=(const bigint_vector& X)
{
  assert(dim == X.dim);
  
  unsigned long i;
  
  for (i = 0; i < dim; i++)
    this->operator[](i) -= X[i];
  
  return *this;
}

ostream& operator <<(ostream& o, const bigint_vector& X)
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

