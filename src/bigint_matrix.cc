#include <bigint_matrix.h>

bigint_matrix :: bigint_matrix(const unsigned long r, const unsigned long c)
  : num_row(r), num_col(c)
{
  unsigned long n;
  
  if (n = num_row * num_col)
    rep = new bigint [n];
  else
    rep = 0;
}

bigint_matrix :: bigint_matrix(const unsigned long r, const unsigned long c,
                               const bigint* X)
  : num_row(r), num_col(c)
{
  unsigned long i, j;
  unsigned long n;
  
  if (n = num_row * num_col)
  {
    rep = new bigint [n];
    
    for (i = 0; i < num_row; i++)
      for (j = 0; j < num_col; j++)
        rep[i * num_col + j] = X[i * num_col + j];
  }
  else
    rep = 0;
}

bigint_matrix :: bigint_matrix(const bigint_matrix& X)
  : num_row(X.num_row), num_col(X.num_col)
{
  unsigned long i, j;
  unsigned long n;
  
  if (n = num_row * num_col)
  {
    rep = new bigint [n];
    
    for (i = 0; i < num_row; i++)
      for (j = 0; j < num_col; j++)
        rep[i * num_col + j] = X(i, j);
  }
  else
    rep = 0;
}

bigint_matrix& bigint_matrix :: operator =(const bigint_matrix& X)
{
  unsigned long i, j;
  unsigned long n;
  
  if (this != &X)
  {
    if (num_row != X.num_row || num_col != X.num_col)
    {
      if (rep)
        delete [] rep;
      
      num_row = X.num_row;
      num_col = X.num_col;
      
      if (n = num_row * num_col)
        rep = new bigint [n];
      else
        rep = 0;
    }
    
    for (i = 0; i < num_row; i++)
      for (j = 0; j < num_col; j++)
        rep[i * num_col + j] = X(i, j);
  }
  
  return *this;
}

bigint_matrix :: ~bigint_matrix()
{
  if (rep)
    delete [] rep;
}

bigint_matrix bigint_matrix :: add(const bigint_matrix& X) const
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long i, j;
  bigint_matrix Y(num_row, num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      Y(i, j) = this->operator ()(i, j) + X(i, j);
  
  return Y;
}

bigint_matrix operator +(const bigint_matrix& X, const bigint_matrix& Y)
{
  return X.add(Y);
}

bigint_matrix bigint_matrix :: sub(const bigint_matrix& X) const
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long i, j;
  bigint_matrix Y(num_row, num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      Y(i, j) = this->operator ()(i, j) - X(i, j);
  
  return Y;
}

bigint_matrix operator -(const bigint_matrix& X, const bigint_matrix& Y)
{
  return X.sub(Y);
}

bigint_matrix bigint_matrix :: mul(const bigint_matrix& X) const
{
  assert(num_col == X.num_row);
  
  unsigned long i, j, k;
  bigint_matrix Y(num_row, X.num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < X.num_col; j++)
    {
      Y(i, j) = 0;
      
      for (k = 0; k < num_col; k++)
        Y(i, j) += this->operator ()(i, k) * X(k, j);
    }
  
  return Y;
}

bigint_matrix operator *(const bigint_matrix& X, const bigint_matrix& Y)
{
  return X.mul(Y);
}

bigint_matrix bigint_matrix :: neg() const
{
  unsigned long i, j;
  bigint_matrix X(num_row, num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      X(i, j) = - this->operator ()(i, j);
  
  return X;
}

bigint_matrix operator -(const bigint_matrix& X)
{
  return X.neg();
}

//bigint_matrix bigint_matrix :: scalar_mul(const bigint& x) const
//{
//  unsigned long i, j;
//  bigint_matrix X(num_row, num_col);
//  
//  for (i = 0; i < num_row; i++)
//    for (j = 0; j < num_col; j++)
//      X(i, j) = this->operator ()(i, j) * x;
//  
//  return X;
//}
//
//bigint_matrix scalar_mul(const bigint& x, const bigint_matrix& Y)
//{
//  return Y.scalr_mul(x);
//}

int bigint_matrix :: cmp(const bigint_matrix& X) const
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long i;
  unsigned long n;
  int           c;
  
  n = num_row * num_col;
  i = 0;
  
  while (i < n && rep[i] ==  X.rep[i])
    i++;
  
  if (i == n)
    c = 0;
  else  //  if (i < n)
    c = 1;
  
  return c;
}

int operator ==(const bigint_matrix& X, const bigint_matrix& Y)
{
  return X.cmp(Y) == 0;
}

int operator !=(const bigint_matrix& X, const bigint_matrix& Y)
{
  return X.cmp(Y) != 0;
}

bigint_matrix& bigint_matrix :: operator +=(const bigint_matrix& X)
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long i, j;
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      this->operator ()(i, j) += X(i, j);
  
  return *this;
}

bigint_matrix& bigint_matrix :: operator -=(const bigint_matrix& X)
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long i, j;
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      this->operator ()(i, j) -= X(i, j);
  
  return *this;
}

//bigint_matrix& bigint_matrix :: scalar_mul_assign(const bigint& x)
//{
//  unsigned long i, j;
//  
//  for (i = 0; i < num_row; i++)
//    for (j = 0; j < num_col; j++)
//      this->operator ()(i, j) *= x;
//  
//  return *this;
//}

int bigint_matrix :: is_identity() const
{
  assert(num_row == num_col);
  
  unsigned long i, j;
  int           e;
  
  e = 1;
  i = j = 0;
  
  while (e && i < num_row)
  {
    j = 0;
    
    while (e && j < num_col)
    {
      if (i == j)
        e = this->operator ()(i, j) == 1;
      else
        e = this->operator ()(i, j) == 0;
      
      j++;
    }
    
    i++;
  }
  
  return e;
}

bigint bigint_matrix :: det() const
{
  assert(num_row == num_col);
  
  long          i;
  unsigned long j, k;
  bigint_matrix X = *this;
  bigint        p, x;
  int           s;
  bigint        d;
  
  for (i = 0, p = 1, s = 1; i < num_col - 1 && s; i++)
  {
    j = i;
    
//    while (j < num_row && sign(X(j, i)) == 0)
    while (j < num_row && sgn(X(j, i)) == 0)
      j++;
    
    if (j < num_row)
    {
      if (j != i)
      {
        s = - s;
        
        for (k = i; k < num_col; k++)
        {
          x       = X(i, k);
          X(i, k) = X(j, k);
          X(j, k) = x;
        }
      }
      
      for (j = i + 1; j < num_row; j++)
      {
        for (k = i + 1; k < num_col; k++)
        {
          X(j, k) = X(i, i) * X(j, k) - X(j, i) * X(i, k);
          X(j, k) /= p;
        }
        
        X(j, i) = 0;
      }
      
      p = X(i, i);
    }
    else
      s = 0;
  }
  
  if (s)
    d = bigint(s) * X(num_row - 1, num_col - 1);
  else
    d = 0;
  
  return d;
}

bigint det(const bigint_matrix& X)
{
  return X.det();
}

ostream& operator <<(ostream& o, const bigint_matrix& X)
{
  unsigned long i, j;
  
  for (i = 0; i < X.num_row; i++)
  {
    for (j = 0; j < X.num_col; j++)
      o << X(i, j) << " ";
    
    o << endl;
  }
  
  return o;
}

