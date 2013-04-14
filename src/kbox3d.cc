//  file:    kbox3d.cc
//  update:  11/30/02

#include <kbox3d.h>

K_BOX3D :: K_BOX3D()
{
  unsigned long i;
  
  for (i = 0; i < 3; i++)
  {
    low[i]        = 0;
    high[i]       = 0;
    low_infty[i]  = - 1;
    high_infty[i] = 1;
  }
}

K_BOX3D :: K_BOX3D(const bigrational& l_s, const bigrational& h_s,
                   const bigrational& l_t, const bigrational& h_t,
                   const bigrational& l_u, const bigrational& h_u)
{
  unsigned long i;
  
  low[0]  = l_s;
  high[0] = h_s;
  low[1]  = l_t;
  high[1] = h_t;
  low[2]  = l_u;
  high[2] = h_u;
  
  for (i = 0; i < 3; i++)
    low_infty[i] = high_infty[i] = 0;
}

K_BOX3D :: K_BOX3D(const bigrational* const l, const bigrational* const h)
{
  unsigned long i;
  
  for (i = 0; i < 3; i++)
  {
    low[i]        = l[i];
    high[i]       = h[i];
    low_infty[i]  = high_infty[i] = 0;
  }
}

K_BOX3D :: K_BOX3D(const K_BOX3D& b)
{
  unsigned long i;
  
  for (i = 0; i < 3; i++)
  {
    low[i]        = b.low[i];
    high[i]       = b.high[i];
    low_infty[i]  = b.low_infty[i];
    high_infty[i] = b.high_infty[i];
  }
}

K_BOX3D& K_BOX3D :: operator =(const K_BOX3D& b)
{
  if (this != &b)
  {
    unsigned long i;
    
    for (i = 0; i < 3; i++)
    {
      low[i]        = b.low[i];
      high[i]       = b.high[i];
      low_infty[i]  = b.low_infty[i];
      high_infty[i] = b.high_infty[i];
    }
  }
  
  return *this;
}

K_BOX3D :: ~K_BOX3D()
{ }

ostream& K_BOX3D :: output(ostream& o) const
{
  o << "[";
  
  if (!low_infty[0])
    o << low[0];
  else if (low_infty[0] < 0)
    o << "- infty";
  else  //  if (low_infty[0] > 0)
    o << "+ infty";
  
  o << ", ";
  
  if (!high_infty[0])
    o << high[0];
  else if (high_infty[0] < 0)
    o << "- infty";
  else  //  if (high_infty[0] > 0)
    o << "+ infty";
  
  o << "] x [";
  
  if (!low_infty[1])
    o << low[1];
  else if (low_infty[1] < 0)
    o << "- infty";
  else  //  if (low_infty[1] > 0)
    o << "+ infty";
  
  o << ", ";
  
  if (!high_infty[1])
    o << high[1];
  else if (high_infty[1] < 0)
    o << "- infty";
  else  //  if (high_infty[1] > 0)
    o << "+ infty";
  
  o << "] x [";
  
  if (!low_infty[2])
    o << low[2];
  else if (low_infty[2] < 0)
    o << "- infty";
  else  //  if (low_infty[2] > 0)
    o << "+ infty";
  
  o << ", ";
  
  if (!high_infty[2])
    o << high[2];
  else if (high_infty[2] < 0)
    o << "- infty";
  else  //  if (high_infty[2] > 0)
    o << "+ infty";
  
  o << "]";
  
  return o;
}

ostream& operator <<(ostream& o, const K_BOX3D& b)
{
  return b.output(o);
}

int K_BOX3D :: overlap(const K_BOX3D& b) const
{
  unsigned long i;
  int           o;
  
  for (o = 1, i = 0; o && i < 3; i++)
    if (low_infty[i] < 0)
    {
      if (high_infty[i] < 0 && b.low_infty[i] >= 0)
        o = 0;
      else if (!high_infty[i])
        if (b.low_infty[i] > 0)
          o = 0;
        else if (!b.low_infty[i] && high[i] < b.low[i])
          o = 0;
    }
    else if (!low_infty[i])
    {
      if (b.high_infty[i] < 0)
        o = 0;
      else if (!b.high_infty[i] && low[i] > b.high[i])
        o = 0;
      else if (!high_infty[i])
        if (b.low_infty[i] > 0)
          o = 0;
        else if (!b.low_infty[i] && high[i] < b.low[i])
          o = 0;
    }
    else  //  if (low_infty[i] > 0)
      if (b.high_infty[i] <= 0)
        o = 0;
  
  return o;
}

int K_BOX3D :: contains(const K_BOX3D& b) const
{
  unsigned long i;
  int           c;
  
  for (c = 1, i = 0; c && i < 3; i++)
    if (low_infty[i] > 0 && b.low_infty[i] <= 0)
      c = 0;
    else if (!low_infty[i])
    {
      if (b.low_infty[i] < 0)
        c = 0;
      else if (!b.low_infty[i] && low[i] > b.low[i])
        c = 0;
    }
    else if (high_infty[i] < 0 && b.high_infty[i] >= 0)
      c = 0;
    else if (!high_infty[i])
    {
      if (b.high_infty[i] > 0)
        c = 0;
      else if (!b.high_infty[i] && high[i] < b.high[i])
        c = 0;
    }
  
  return c;
}

