//  file:    bigrational_matrix.cc
//  update:  09/25/02

#include <bigrational_matrix.h>

bigrational_matrix :: bigrational_matrix(const unsigned long r,
                                         const unsigned long c)
  : num_row(r), num_col(c)
{
  unsigned long n;
  
  if (n = num_row * num_col)
    rep = new bigrational [n];
  else
    rep = 0;
}

bigrational_matrix :: bigrational_matrix(const unsigned long r,
                                         const unsigned long c,
                                         const bigrational* X)
  : num_row(r), num_col(c)
{
  unsigned long i, j;
  unsigned long n;
  
  if (n = num_row * num_col)
  {
    rep = new bigrational [n];
    
    for (i = 0; i < num_row; i++)
      for (j = 0; j < num_col; j++)
        rep[i * num_col + j] = X[i * num_col + j];
  }
  else
    rep = 0;
}

bigrational_matrix :: bigrational_matrix(const bigrational_matrix& X)
  : num_row(X.num_row), num_col(X.num_col)
{
  unsigned long i, j;
  unsigned long n;
  
  if (n = num_row * num_col)
  {
    rep = new bigrational [n];
    
    for (i = 0; i < num_row; i++)
      for (j = 0; j < num_col; j++)
        rep[i * num_col + j] = X(i, j);
  }
  else
    rep = 0;
}

//bigrational_matrix :: bigrational_matrix(const bigint_matrix& X)
//  : num_row(X.num_row), num_col(X.num_col)
//{
//  unsigned long i, j;
//  unsigned long n;
//  
//  if (n = num_row * num_col)
//  {
//    rep = new bigrational [n];
//    
//    for (i = 0; i < num_row; i++)
//      for (j = 0; j < num_col; j++)
//        rep[i * num_col + j] = X(i, j);
//  }
//  else
//    rep = 0;
//}

bigrational_matrix&
bigrational_matrix :: operator =(const bigrational_matrix& X)
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
        rep = new bigrational [n];
      else
        rep = 0;
    }
    
    for (i = 0; i < num_row; i++)
      for (j = 0; j < num_col; j++)
        rep[i * num_col + j] = X(i, j);
  }
  
  return *this;
}

bigrational_matrix :: ~bigrational_matrix()
{
  if (rep)
    delete [] rep;
}

bigrational_matrix bigrational_matrix :: add(const bigrational_matrix& X) const
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long      i, j;
  bigrational_matrix Y(num_row, num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      Y(i, j) = this->operator ()(i, j) + X(i, j);
  
  return Y;
}

bigrational_matrix operator +(const bigrational_matrix& X,
                              const bigrational_matrix& Y)
{
  return X.add(Y);
}

bigrational_matrix bigrational_matrix :: sub(const bigrational_matrix& X) const
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long      i, j;
  bigrational_matrix Y(num_row, num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      Y(i, j) = this->operator ()(i, j) - X(i, j);
  
  return Y;
}

bigrational_matrix operator -(const bigrational_matrix& X,
                              const bigrational_matrix& Y)
{
  return X.sub(Y);
}

bigrational_matrix bigrational_matrix :: mul(const bigrational_matrix& X) const
{
  assert(num_col == X.num_row);
  
  unsigned long      i, j, k;
  bigrational_matrix Y(num_row, X.num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < X.num_col; j++)
    {
      Y(i, j) = 0;
      
      for (k = 0; k < num_col; k++)
        Y(i, j) += this->operator ()(i, k) * X(k, j);
    }
  
  return Y;
}

bigrational_matrix operator *(const bigrational_matrix& X,
                              const bigrational_matrix& Y)
{
  return X.mul(Y);
}

bigrational_matrix bigrational_matrix :: neg() const
{
  unsigned long      i, j;
  bigrational_matrix X(num_row, num_col);
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      X(i, j) = - this->operator ()(i, j);
  
  return X;
}

bigrational_matrix operator -(const bigrational_matrix& X)
{
  return X.neg();
}

//bigrational_matrix bigrational_matrix :: scalar_mul(const bigrational& x) const
//{
//  unsigned long      i, j;
//  bigrational_matrix X(num_row, num_col);
//  
//  for (i = 0; i < num_row; i++)
//    for (j = 0; j < num_col; j++)
//      X(i, j) = this->operator ()(i, j) * x;
//  
//  return X;
//}
//
//bigrational_matrix scalar_mul(const bigrational& x,
//                              const bigrational_matrix& Y)
//{
//  return Y.scalr_mul(x);
//}

int bigrational_matrix :: cmp(const bigrational_matrix& X) const
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

int operator ==(const bigrational_matrix& X, const bigrational_matrix& Y)
{
  return X.cmp(Y) == 0;
}

int operator !=(const bigrational_matrix& X, const bigrational_matrix& Y)
{
  return X.cmp(Y) != 0;
}

bigrational_matrix&
bigrational_matrix :: operator +=(const bigrational_matrix& X)
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long i, j;
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      this->operator ()(i, j) += X(i, j);
  
  return *this;
}

bigrational_matrix&
bigrational_matrix :: operator -=(const bigrational_matrix& X)
{
  assert(num_row == X.num_row && num_col == X.num_col);
  
  unsigned long i, j;
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      this->operator ()(i, j) -= X(i, j);
  
  return *this;
}

//bigrational_matrix&
//bigrational_matrix :: scalar_mul_assign(const bigrational& x)
//{
//  unsigned long i, j;
//  
//  for (i = 0; i < num_row; i++)
//    for (j = 0; j < num_col; j++)
//      this->operator ()(i, j) *= x;
//  
//  return *this;
//}

int bigrational_matrix :: is_identity() const
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

bigrational bigrational_matrix :: det() const
{
  assert(num_row == num_col);
  
  long               i;
  unsigned long      j, k;
  bigrational_matrix X = *this;
  bigrational        p, x;
  int                s;
  bigrational        d;
  
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
    d = bigrational(s) * X(num_row - 1, num_col - 1);
  else
    d = 0;
  
  return d;
}

bigrational det(const bigrational_matrix& X)
{
  return X.det();
}

bigrational bigrational_matrix :: Bareiss(bigrational_matrix& Y) const
{
  assert(num_row == num_col && num_row == Y.num_row && num_col == Y.num_col);
  
  long               i, j, k;
  bigrational_matrix X = *this;
  bigrational        p, q, x, y;
  int                s;
  bigrational        d;
  
  for (i = 0; i < num_row; i++)
    for (j = 0; j < num_col; j++)
      if (i == j)
        Y(i, j) = 1;
      else
        Y(i, j) = 0;
  
  cerr << X << endl;
  cerr << Y << endl;
  cerr << " upper triangulate " << endl;
  
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
        
        for (k = 0; k < num_col; k++)
        {
          y       = Y(i, k);
          Y(i, k) = Y(j, k);
          Y(j, k) = y;
        }
        
        for (k = i; k < num_col; k++)
        {
          x       = X(i, k);
          X(i, k) = X(j, k);
          X(j, k) = x;
        }
      }
      
      for (j = i + 1; j < num_row; j++)
      {
        for (k = 0; k < num_col; k++)
        {
          Y(j, k) = X(i, i) * Y(j, k) - X(j, i) * Y(i, k);
          Y(j, k) /= p;
        }
        
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
    
    cerr << X << endl;
    cerr << Y << endl;
  }
  
  if (s)
    d = bigrational(s) * X(num_row - 1, num_col - 1);
  else
    d = 0;
  
  if (s)
  {
    cerr << " diagonalize " << endl;
    
    for (i = num_col - 1, q = 1; i >= 0; i--)
    {
      for (j = i - 1; j >= 0; j--)
      {
        for (k = num_col - 1; k >= 0; k--)
        {
          Y(j, k) = X(i, i) * Y(j, k) - X(j, i) * Y(i, k);
          Y(j, k) /= q;
        }
        
        for (k = i - 1; k >= 0; k--)
        {
          X(j, k) = X(i, i) * X(j, k) - X(j, i) * X(i, k);
          X(j, k) /= q;
        }
        
        X(j, i) = 0;
      }
      
      q = X(i, i);
      
      cerr << X << endl;
      cerr << Y << endl;
    }
    
    for (i = 0; i < num_row; i++)
      for (j = 0; j < num_col; j++)
        Y(i, j) /= X(i, i);
    
    cerr << Y << endl;
    
    bigrational_matrix Z = this->mul(Y);
    cerr << *this << endl;
    cerr << Y << endl;
    cerr << Z << endl;
  }
  
  return d;
}

bigrational_matrix bigrational_matrix :: inverse() const
{
  bigrational        d;
  bigrational_matrix X(num_row, num_col);
  
  d = this->Bareiss(X);
  
  return X;
}

unsigned long bigrational_matrix :: rank() const
{
  unsigned long      i, j, k;
  bigrational_matrix X = *this;
  bigrational        x;
  unsigned long      r;
  
//  cerr << " before: " << endl << flush;
//  cerr << X << endl << flush;
  
  r = 0;
  i = 0;
  
  while (i < num_col)
  {
    j = r;
    
    while (j < num_row && !sgn(X(j, i)))
      j++;
    
    if (j < num_row)
    {
      if (j != r)
        for (k = r; k < num_col; k++)
        {
          x = X(r, k);
          X(r, k) = X(j, k);
          X(j, k) = x;
        }
      
      for (j = r + 1; j < num_row; j++)
      {
        for (k = i + 1; k < num_col; k++)
          X(j, k) = X(j, k) - X(j, i) * X(r, k) / X(r, i);
        
        X(j, i) = 0;
      }
      
      r++;
    }
    
    i++;
  }
  
//  cerr << " after: " << endl << flush;
//  cerr << X << endl << flush;
  
  return r;
}

unsigned long rank(const bigrational_matrix& X)
{
  return X.rank();
}

ostream& operator <<(ostream& o, const bigrational_matrix& X)
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

