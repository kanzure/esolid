#include <kboxco2.h>

K_BOXCO2 :: K_BOXCO2()
{
  unsigned long i;
  
  for (i = 0; i < 2; i++)
  {
    low[i]       = 0;
    high[i]      = 0;
    low_open[i]  = 1;
    high_open[i] = 1;
  }
}

K_BOXCO2 :: K_BOXCO2(const bigrational& l_s, const bigrational& h_s,
                     const bigrational& l_t, const bigrational& h_t,
                     const int l_open_s, const int h_open_s,
                     const int l_open_t, const int h_open_t)
{
  low[0]       = l_s;
  high[0]      = h_s;
  low[1]       = l_t;
  high[1]      = h_t;
  low_open[0]  = l_open_s;
  high_open[0] = h_open_s;
  low_open[1]  = l_open_t;
  high_open[1] = h_open_t;
}

K_BOXCO2 :: K_BOXCO2(const bigrational* const l, const bigrational* const h,
                     const int* const l_open, const int* const h_open)
{
  unsigned long i;
  
  for (i = 0; i < 2; i++)
  {
    low[i]       = l[i];
    high[i]      = h[i];
    low_open[i]  = l_open[i];
    high_open[i] = h_open[i];
  }
}

K_BOXCO2 :: K_BOXCO2(const K_BOXCO2& b)
{
  unsigned long i;
  
  for (i = 0; i < 2; i++)
  {
    low[i]       = b.low[i];
    high[i]      = b.high[i];
    low_open[i]  = b.low_open[i];
    high_open[i] = b.high_open[i];
  }
}

K_BOXCO2& K_BOXCO2 :: operator =(const K_BOXCO2& b)
{
  if (this != &b)
  {
    unsigned long i;
    
    for (i = 0; i < 2; i++)
    {
      low[i]       = b.low[i];
      high[i]      = b.high[i];
      low_open[i]  = b.low_open[i];
      high_open[i] = b.high_open[i];
    }
  }
  
  return *this;
}

K_BOXCO2 :: ~K_BOXCO2()
{ }

ostream& K_BOXCO2 :: output(ostream& o) const
{
  if (low_open[0])
    o << "( ";
  else  //  if (low_open[0] == 0)
    o << "[ ";
  
  o << low[0] << ", " << high[0];
  
  if (high_open[0])
    o << " )";
  else  //  if (high_open[0] == 0)
    o << " ]";
  
  o << " x ";
  
  if (low_open[1])
    o << "( ";
  else  //  if (low_open[1] == 0)
    o << "[ ";
  
  o << low[1] << ", " << high[1];
  
  if (high_open[1])
    o << " )";
  else  //  if (high_open[1] == 0)
    o << " ]";
  
  return o;
}

ostream& operator <<(ostream& o, const K_BOXCO2& b)
{
  return b.output(o);
}

bigrational K_BOXCO2 :: get_low_s() const
{
  return low[0];
}

bigrational K_BOXCO2 :: get_high_s() const
{
  return high[0];
}

bigrational K_BOXCO2 :: get_low_t() const
{
  return low[1];
}

bigrational K_BOXCO2 :: get_high_t() const
{
  return high[1];
}

int K_BOXCO2 :: is_low_s_open() const
{
  return low_open[0];
}

int K_BOXCO2 :: is_high_s_open() const
{
  return high_open[0];
}

int K_BOXCO2 :: is_low_t_open() const
{
  return low_open[1];
}

int K_BOXCO2 :: is_high_t_open() const
{
  return high_open[1];
}

K_BOXCO2 K_BOXCO2 :: merge(const K_BOXCO2& b) const
{
  unsigned long i;
  K_BOXCO2      c;
  
  for (i = 0; i < 2; i++)
  {
    if (low[i] < b.low[i])
    {
      c.low[i]      = low[i];
      c.low_open[i] = low_open[i];
    }
    else if (low[i] > b.low[i])
    {
      c.low[i]      = b.low[i];
      c.low_open[i] = b.low_open[i];
    }
    else  //  if (low[i] == b.low[i])
    {
      c.low[i]      = low[i];
      c.low_open[i] = low_open[i] * b.low_open[i];
    }
    
    if (high[i] > b.high[i])
    {
      c.high[i]      = high[i];
      c.high_open[i] = high_open[i];
    }
    else if (high[i] < b.high[i])
    {
      c.high[i]      = b.high[i];
      c.high_open[i] = b.high_open[i];
    }
    else  //  if (high[i] == b.high[i])
    {
      c.high[i]      = high[i];
      c.high_open[i] = high_open[i] * b.high_open[i];
    }
  }
  
  return c;
}

int K_BOXCO2 :: overlap(const K_BOXCO2& b) const
{
  unsigned long i;
  int           o;
  
  for (o = 1, i = 0; o && i < 2; i++)
    if (low[i] > b.high[i]
        ||
        high[i] < b.low[i]
        ||
        low[i] == b.high[i] && (low_open[i] || b.high_open[i])
        ||
        high[i] == b.low[i] && (high_open[i] || b.low_open[i]))
      o = 0;
  
  return o;
}

int K_BOXCO2 :: contains(const K_BOXCO2& b) const
{
  unsigned long i;
  int           c;
  
  for (i = 0, c = 1; c && i < 2; i++)
    if (low[i] > b.low[i]
        ||
        high[i] < b.high[i]
        ||
        low[i] == b.low[i] && low_open[i] && !b.low_open[i]
        ||
        high[i] == b.high[i] && high_open[i] && !b.high_open[i])
      c = 0;
  
  return c;
}

