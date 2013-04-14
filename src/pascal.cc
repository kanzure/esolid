//  file:    pascal.cc
//  update:  12/07/02

#include <pascal.h>

Pascal :: Pascal()
{
  num_row = 0;
  num_col = 0;
  tab     = 0;
}

Pascal :: Pascal(const unsigned long r, const unsigned long c)
{
  unsigned long i, j;
  
  num_row = r;
  num_col = c;
  
  if (num_row * num_col > 0)
  {
    tab = new long [num_row * num_col];
    
    for (i = 0; i < num_col; i++)
      tab[i] = 1;
    
    for (i = 1; i < num_row; i++)
      tab[i * num_col] = 1;
    
    for (i = 1; i < num_row; i++)
      for (j = 1; j < num_col; j++)
        tab[i * num_col + j] =
          tab[i * num_col + j - 1] + tab[(i - 1) * num_col + j];
  }
  else
    tab = 0;
}

Pascal :: ~Pascal()
{
  if (num_row * num_col > 0)
    delete [] tab;
}

int Pascal :: resize(const unsigned long r, const unsigned long c)
{
  unsigned long i, j;
  
  if (num_row * num_col > 0)
    delete [] tab;
  
  if (num_row < r)
    num_row = r;
  
  if (num_col < c)
    num_col = c;
  
  if (num_row * num_col > 0)
  {
    tab = new long [num_row * num_col];
    
    for (i = 0; i < num_col; i++)
      tab[i] = 1;
    
    for (i = 1; i < num_row; i++)
      tab[i * num_col] = 1;
    
    for (i = 1; i < num_row; i++)
      for (j = 1; j < num_col; j++)
        tab[i * num_col + j] =
          tab[i * num_col + j - 1] + tab[(i - 1) * num_col + j];
  }
  else  //  if (!num_row && !num_col)
    tab = 0;
  
  return 0;
}

long Pascal :: get_Pascal(const unsigned long i, const unsigned long j)
{
  if (i + 1 > num_row || j + 1 > num_col)
    resize(i + 1, j + 1);
  
  assert(num_row * num_col > 0);
  
  return tab[i * num_col + j];
}

