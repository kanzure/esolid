#include <config.h>
#ifdef _EXPERIMENT
#include <counter.h>
#endif

#include <kpoint2d.h>

#include <fpconversion.h>

//  K_POINT2D :: K_POINT2D()
//    constructs a type 0 instance.

K_POINT2D :: K_POINT2D()
{
  type = 0;
  
  PtRs = 0;
  PtRt = 0;
  PtBs = 0;
  PtBt = 0;
  
  poly1 = 0;
  poly2 = 0;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const bigrational& b_s, const bigrational& b_t)
//    construct a type 4 instance from bigrational b_s and b_t.

K_POINT2D :: K_POINT2D(const bigrational& b_s, const bigrational& b_t)
{
  type = 4;
  
  PtRs = 0;
  PtRt = 0;
  PtBs = b_s;
  PtBt = b_t;
  
  poly1 = new K_RATPOLY(1, 0, b_s);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(1, 0, b_t);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const ROOT1& r_s, const ROOT1& r_t,
//                         const K_RATPOLY& P1, const K_RATPOLY& P2)
//    construct a type 2 instance from
//      ROOT1 r_s and r_t and K_RATPOLY P1 and P2.
//    r_s.num_roots and r_t.num_roots must have been 1 and
//    P1 must intersect with P2 at (r_s, r_t).

K_POINT2D :: K_POINT2D(const ROOT1& r_s, const ROOT1& r_t,
                       const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  assert(r_s.num_roots == 1);
  assert(r_t.num_roots == 1);
  
  type = 1;
  
  PtRs = new ROOT1(r_s);
  PtRs->ref_count++;
  PtRt = new ROOT1(r_t);
  PtRt->ref_count++;
  
  poly1 = new K_RATPOLY(P1);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(P2);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const ROOT1& r_s, const bigrational& b_t,
//                         const K_RATPOLY& P1, const K_RATPOLY& P2)
//    construct a type 2 instance from
//      bigrational b_s, ROOT1 r_t and K_RATPOLY P1 and P2.
//    r_s.num_roots must have been 1 and
//    P1 must intersect with P2 at (r_s, b_t).

K_POINT2D :: K_POINT2D(const ROOT1& r_s, const bigrational& b_t,
                       const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  assert(r_s.num_roots == 1);
  
  type = 2;
  
  PtRs = new ROOT1(r_s);
  PtRs->ref_count++;
  PtRt = 0;
  PtBt = b_t;
  
  poly1 = new K_RATPOLY(P1);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(P2);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const bigrational& b_s, const ROOT1& r_t,
//                         const K_RATPOLY& P1, const K_RATPOLY& P2)
//    construct a type 3 instance from
//      bigrational b_s, ROOT1 r_t and K_RATPOLY P1 and P2.
//    r_t.num_roots must have been 1 and
//    P1 must intersect with P2 at (b_s, r_t).

K_POINT2D :: K_POINT2D(const bigrational& b_s, const ROOT1& r_t,
                       const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  assert(r_t.num_roots == 1);
  
  type = 3;
  
  PtRs = 0;
  PtRt = new ROOT1(r_t);
  PtRt->ref_count++;
  PtBs = b_s;
  
  poly1 = new K_RATPOLY(P1);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(P2);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const bigrational& b_s, const bigrational& b_t,
//                         const K_RATPOLY& P1, const K_RATPOLY& P2)
//    construct a type 4 instance from
//      bigrational b_s and b_t and K_RATPOLY P1 and P2.
//    P1 must intersect with P2 at (b_s, b_t).

K_POINT2D :: K_POINT2D(const bigrational& b_s, const bigrational& b_t,
                       const K_RATPOLY& P1, const K_RATPOLY& P2)
{
//  bigrational b[2];
//  
//  b[0] = b_s;
//  b[1] = b_t;
//  
//  assert(!sgn(P1.evaluate(b)) && !sgn(P2.evaluate(b)));
//  
  type = 4;
  
  PtRs = 0;
  PtRt = 0;
  PtBs = b_s;
  PtBt = b_t;
  
  poly1 = new K_RATPOLY(P1);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(P2);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const bigrational& b_s, const K_POINT1D& x_t,
//                         const K_RATPOLY& P1, const K_RATPOLY& P2)
//    construct a type 3 or 4 instance from
//      bigrational b_s, K_POINT1D x_t and K_RATPOLY P1 and P2.
//    P1 must intersect with P2 at (b_s, x_t).

K_POINT2D :: K_POINT2D(const bigrational& b_s, const K_POINT1D& x_t,
                       const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  assert(x_t.type > 0);
  
  if (x_t.type == 1)
  {
    type = 3;
    
    PtRs = 0;
    PtRt = x_t.PtR;
    PtRt->ref_count++;
    PtBs = b_s;
  }
  else  //  if (x_t.type == 2)
  {
    type = 4;
    
    PtRs = 0;
    PtRt = 0;
    PtBs = b_s;
    PtBt = x_t.PtB;
  }
  
  poly1 = new K_RATPOLY(P1);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(P2);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const K_POINT1D& x_s, const bigrational& b_t,
//                         const K_RATPOLY& P1, const K_RATPOLY& P2)
//    construct a type 2 or 4 instance from
//      K_POINT1D x_s, bigrational b_t and K_RATPOLY P1 and P2.
//    P1 must intersect with P2 at (x_s, b_t).

K_POINT2D :: K_POINT2D(const K_POINT1D& x_s, const bigrational& b_t,
                       const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  assert(x_s.type > 0);
  
  if (x_s.type == 1)
  {
    type = 2;
    
    PtRs = x_s.PtR;
    PtRs->ref_count++;
    PtRt = 0;
    PtBt = b_t;
  }
  else  //  if (x_s.type == 2)
  {
    type = 4;
    
    PtRs = 0;
    PtRt = 0;
    PtBs = x_s.PtB;
    PtBt = b_t;
  }
  
  poly1 = new K_RATPOLY(P1);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(P2);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const K_POINT1D& x_s, const K_POINT1D& x_t,
//                         const K_RATPOLY& P1, const K_RATPOLY& P2)
//    construct a type 1 or 2 or 3 or 4 instance from
//      K_POINT1D x_s and x_t and K_RATPOLY P1 and P2.
//    P1 must intersect with P2 at (x_s, x_t).

K_POINT2D :: K_POINT2D(const K_POINT1D& x_s, const K_POINT1D& x_t,
                       const K_RATPOLY& P1, const K_RATPOLY& P2)
{
  assert(x_s.type > 0);
  assert(x_t.type > 0);
  
  if (x_s.type == 1 && x_t.type == 1)
  {
    type = 1;
    
    PtRs = x_s.PtR;
    PtRs->ref_count++;
    PtRt = x_t.PtR;
    PtRt->ref_count++;
  }
  else if (x_s.type == 1 && x_t.type == 2)
  {
    type = 2;
    
    PtRs = x_s.PtR;
    PtRs->ref_count++;
    PtRt = 0;
    PtBt = x_t.PtB;
  }
  else if (x_s.type == 2 && x_t.type == 1)
  {
    type = 3;
    
    PtRs = 0;
    PtRt = x_t.PtR;
    PtRt->ref_count++;
    PtBs = x_s.PtB;
  }
  else  //  if (x_s.type == 2 && x_t.type == 2)
  {
    type = 4;
    
    PtRs = 0;
    PtRt = 0;
    PtBs = x_s.PtB;
    PtBt = x_t.PtB;
  }
  
  poly1 = new K_RATPOLY(P1);
  poly1->ref_count++;
  poly2 = new K_RATPOLY(P2);
  poly2->ref_count++;
  
  prev = next = this;
  
  pt_in_other_dom = 0;
  
  ref_count = 0;
}

//  K_POINT2D :: K_POINT2D(const K_POINT2D& x)
//    the copy constructor

K_POINT2D :: K_POINT2D(const K_POINT2D& x)
{
  type = x.type;
  
  if (PtRs = x.PtRs)
    PtRs->ref_count++;
  
  if (PtRt = x.PtRt)
    PtRt->ref_count++;
  
  PtBs = x.PtBs;
  PtBt = x.PtBt;
  
  if (poly1 = x.poly1)
    poly1->ref_count++;
  
  if (poly2 = x.poly2)
    poly2->ref_count++;
  
//  next         = x.next;
//  prev         = x.next->prev;
//  x.next->prev = this;
//  x.next       = this;
  prev = next = this;
  
  pt_in_other_dom = x.pt_in_other_dom;
  
  ref_count = 0;
}

//  K_POINT2D& K_POINT2D :: operator =(const K_POINT2D& x)
//    the assignment operator

K_POINT2D& K_POINT2D :: operator =(const K_POINT2D& x)
{
  if (this != &x)
  {
    if (PtRs && !--PtRs->ref_count)
      delete PtRs;
    
    if (PtRt && !--PtRt->ref_count)
      delete PtRt;
    
    if (poly1 && !--poly1->ref_count)
      delete poly1;
    
    if (poly2 && !--poly2->ref_count)
      delete poly2;
    
    next->prev = prev;
    prev->next = next;
    
    if (pt_in_other_dom)
      pt_in_other_dom->pt_in_other_dom = 0;
    
    type = x.type;
    
    if (PtRs = x.PtRs)
      PtRs->ref_count++;
    
    if (PtRt = x.PtRt)
      PtRt->ref_count++;
    
    PtBs = x.PtBs;
    PtBt = x.PtBt;
    
    if (poly1 = x.poly1)
      poly1->ref_count++;
    
    if (poly2 = x.poly2)
      poly2->ref_count++;
    
//    next         = x.next;
//    prev         = x.next->prev;
//    x.next->prev = this;
//    x.next       = this;
    prev = next = this;
    
    pt_in_other_dom = x.pt_in_other_dom;
  }
  
  return *this;
}

//  K_POINT2D :: ~K_POINT2D()
//    the destructor

K_POINT2D :: ~K_POINT2D()
{
  if (PtRs && !--PtRs->ref_count)
    delete PtRs;
  
  if (PtRt && !--PtRt->ref_count)
    delete PtRt;
  
  if (poly1 && !--poly1->ref_count)
    delete poly1;
  
  if (poly2 && !--poly2->ref_count)
    delete poly2;
  
  next->prev = prev;
  prev->next = next;
  
  if (pt_in_other_dom)
    pt_in_other_dom->pt_in_other_dom = 0;
}

ostream& K_POINT2D :: output(ostream& o) const
{
  if (type == 1)
    o << "(" << *PtRs << ", " << *PtRt << ")" << flush;
  else if (type == 2)
    o << "(" << *PtRs << ", " << PtBt << ")" << flush;
  else if (type == 3)
    o << "(" << PtBs << ", " << *PtRt << ")" << flush;
  else if (type == 4)
    o << "(" << PtBs << ", " << PtBt << ")" << flush;
  else  //  if type == 0
    o << " NULL " << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const K_POINT2D& x)
{
  return x.output(o);
}

bigrational K_POINT2D :: get_low_s() const
{
  assert(type);
  
  bigrational l;
  
  if (type == 1 || type == 2)
    l = PtRs->low;
  else  //  if (type == 3 || type == 4)
    l = PtBs;
  
  return l;
}

bigrational K_POINT2D :: get_high_s() const
{
  assert(type);
  
  bigrational h;
  
  if (type == 1 || type == 2)
    h = PtRs->high;
  else  //  if (type == 3 || type == 4)
    h = PtBs;
  
  return h;
}

bigrational K_POINT2D :: get_low_t() const
{
  assert(type);
  
  bigrational l;
  
  if (type == 1 || type == 3)
    l = PtRt->low;
  else  //  if (type == 2 || type == 4)
    l = PtBt;
  
  return l;
}

bigrational K_POINT2D :: get_high_t() const
{
  assert(type);
  
  bigrational h;
  
  if (type == 1 || type == 3)
    h = PtRt->high;
  else  //  if (type == 2 || type == 4)
    h = PtBt;
  
  return h;
}

bigrational_vector K_POINT2D :: get_low() const
{
  assert(type);
  
  bigrational_vector l(2);
  
  if (type == 1)
  {
    l[0] = PtRs->low;
    l[1] = PtRt->low;
  }
  else if (type == 2)
  {
    l[0] = PtRs->low;
    l[1] = PtBt;
  }
  else if (type == 3)
  {
    l[0] = PtBs;
    l[1] = PtRt->low;
  }
  else  //  if (type == 4)
  {
    l[0] = PtBs;
    l[1] = PtBt;
  }
  
  return l;
}

bigrational_vector K_POINT2D :: get_high() const
{
  assert(type);
  
  bigrational_vector h(2);
  
  if (type == 1)
  {
    h[0] = PtRs->high;
    h[1] = PtRt->high;
  }
  else if (type == 2)
  {
    h[0] = PtRs->high;
    h[1] = PtBt;
  }
  else if (type == 3)
  {
    h[0] = PtBs;
    h[1] = PtRt->high;
  }
  else  //  if (type == 4)
  {
    h[0] = PtBs;
    h[1] = PtBt;
  }
  
  return h;
}

//  int K_POINT2D :: cut_s(const bigrational& b_s) const
//    cut *this by the line s = b_s, i.e.,
//    refine the box for *this
//      by setting get_low_s() or get_high_s() to be b_s.

int K_POINT2D :: cut_s(const bigrational& b_s) const
{
  assert(type > 0);
  
  if ((type == 1 || type == 2) && b_s > PtRs->low && b_s < PtRs->high)
  {
    assert(PtRs->num_roots == 1);
    
    if (!PtRs->poly->sgn_at(b_s))
    {
      if (type == 1)
        type = 3;
      else  //  if (type == 2)
        type = 4;
      
      if (!--PtRs->ref_count)
        delete PtRs;
      
      PtRs = 0;
      PtBs = b_s;
    }
    else  //  if (PtRs->poly->sgn_at(b_s))
      PtRs->cut(b_s);
  }
  
  return 0;
}

//  int K_POINT2D :: cut_t(const bigrational& b_t) const
//    cut *this by the line t = b_t, i.e.,
//    refine the box for *this
//      by setting get_low_t() or get_high_t() to be b_t.

int K_POINT2D :: cut_t(const bigrational& b_t) const
{
  assert(type);
  
  if ((type == 1 || type == 3) && b_t > PtRt->low && b_t < PtRt->high)
  {
    assert(PtRt->num_roots == 1);
    
    if (!PtRt->poly->sgn_at(b_t))
    {
      if (type == 1)
        type = 2;
      else  //  if (type == 3)
        type = 4;
      
      if (!--PtRt->ref_count)
        delete PtRt;
      
      PtRt = 0;
      PtBt = b_t;
    }
    else  //  if (PtRt->poly->sgn_at(b_t))
      PtRt->cut(b_t);
  }
  
  return 0;
}

//  int K_POINT2D :: halve_s() const
//    cut *this by the line s = half that halves the box for *this, i.e.,
//    refine the box for *this
//      by setting get_low_s() or get_high_s() to be half.

int K_POINT2D :: halve_s() const
{
  assert(type);
  
  if (type == 1 || type == 2)
  {
    assert(PtRs->num_roots == 1);
    
    bigrational half = (PtRs->low + PtRs->high) / 2;
    
    if (!PtRs->poly->sgn_at(half))
    {
      if (type == 1)
        type = 3;
      else  //  if (type == 2)
        type = 4;
      
      if (!--PtRs->ref_count)
        delete PtRs;
      
      PtRs = 0;
      PtBs = half;
    }
    else  //  if (PtRs->poly->sgn_at(half))
      PtRs->cut(half);
  }
  
  return 0;
}

//  int K_POINT2D :: halve_t() const
//    cut *this by the line t = half that halves the box for *this, i.e.,
//    refine the box for *this
//      by setting get_low_t() or get_high_t() to be half.

int K_POINT2D :: halve_t() const
{
  assert(type);
  
  if (type == 1 || type == 3)
  {
    assert(PtRt->num_roots == 1);
    
    bigrational half = (PtRt->low + PtRt->high) / 2;
    
    if (!PtRt->poly->sgn_at(half))
    {
      if (type == 1)
        type = 2;
      else  //  if (type == 3)
        type = 4;
      
      if (!--PtRt->ref_count)
        delete PtRt;
      
      PtRt = 0;
      PtBt = half;
    }
    else  //  if (PtRt->poly->sgn_at(half))
      PtRt->cut(half);
  }
  
  return 0;
}

//  int K_POINT2D :: reduce_s(const unsigned long num_bits) const
//    reduce *this s.t.
//      [get_low_s(), get_high_s()] will contain all the numbers
//        approximated by PtRs->float_est to num_bits precision.
//    return 1 if some reduction occurs and
//           0 otherwise.

int K_POINT2D :: reduce_s(const unsigned long num_bits) const
{
  assert(type);
  
  int reduced;
  
  if (type == 1 || type == 2)
  {
    assert(PtRs->num_roots == 1);
    assert(PtRs->ok_float);
    
    bigrational l, h;
    
    //  Compute l and h s.t. [l, h] containing all the numbers approximated
    //  by float_est to num_bits precision.
    
    double2bigrational_interval(PtRs->float_est, num_bits, l, h);
    
    //  If [low, high] \sebseteq [l, h], or
    //     (h >=) l > high (>= low) or (l <=) h < low (<= high) then
    //  nothing to do.
    //  Otherwise, see if reducible.
    
    if (l <= PtRs->low && h >= PtRs->high)
      reduced = 1;
    else if (l <= PtRs->high && h >= PtRs->low)
    //  if (PtRs->low < l <= PtRs->high && PtRs->low <= h < PtRs->high)
    {
      int sgn_l, sgn_h;
      
      reduced = 1;
      
      if (l > PtRs->low)
        if (sgn_l = PtRs->poly->sgn_at(l))
          if (PtRs->sig_low < 0 && sgn_l < 0
              ||
              PtRs->sig_low > 0 && sgn_l > 0
              ||
              PtRs->sig_low == 0 &&
              PtRs->poly->num_Sturm_seq_perm(l) == PtRs->num_perm_low)
          //  Keep reduced == 1.
            PtRs->low = l;
          else
            reduced = 0;
        else  //  if l is a root of PtRs->poly
        //  Keep reduced == 1.
        {
          if (type == 1)
            type = 3;
          else  //  if (type == 2)
            type = 4;
          
          if (!--PtRs->ref_count)
            delete PtRs;
          
          PtRs = 0;
          PtBs = l;
        }
      
      //  Do the following if PtRs != 0 (and even if reduced == 0)
      //  since better high or low might be obtained.
      
      if (PtRs && h < PtRs->high)
        if (sgn_h = PtRs->poly->sgn_at(h))
          if (PtRs->sig_low < 0 && sgn_h > 0
              ||
              PtRs->sig_low > 0 && sgn_h < 0
              ||
              PtRs->sig_low == 0 &&
              PtRs->poly->num_Sturm_seq_perm(h) == PtRs->num_perm_high)
          //  Keep reduced as it is.
            PtRs->high = h;
          else
            reduced = 0;
        else  //  if h is a root of PtRs->poly
        {
          if (type == 1)
            type = 3;
          else  //  if (type == 2)
            type = 4;
          
          if (!--PtRs->ref_count)
            delete PtRs;
          
          PtRs = 0;
          PtBs = h;
          
          reduced = 1;
        }
    }
    else  //  if (h >= l) > high (>= low) or (l <=) h < low (<= high)
      reduced = 0;
  }
  else  //  if (type == 3 || type == 4)
    reduced = 1;
  
  return reduced;
}

//  int K_POINT2D :: reduce_t(const unsigned long num_bits) const
//    reduce *this s.t.
//      [get_low_t(), get_high_t()] will contain all the numbers
//        approximated by PtRt->float_est to num_bits precision.
//    return 1 if some reduction occurs and
//           0 otherwise.

int K_POINT2D :: reduce_t(const unsigned long num_bits) const
{
  assert(type);
  
  int reduced;
  
  if (type == 1 || type == 3)
  {
    assert(PtRt->num_roots == 1);
    assert(PtRt->ok_float);
    
    bigrational l, h;
    
    //  Compute l and h s.t. [l, h] containing all the numbers approximated
    //  by float_est to num_bits precision.
    
    double2bigrational_interval(PtRt->float_est, num_bits, l, h);
    
    //  If [low, high] \sebseteq [l, h], or
    //     (h >=) l > high (>= low) or (l <=) h < low (<= high) then
    //  nothing to do.
    //  Otherwise, see if reducible.
    
    if (l <= PtRt->low && h >= PtRt->high)
      reduced = 1;
    else if (l <= PtRt->high && h >= PtRt->low)
    //  if (PtRt->low < l <= PtRt->high) or (PtRt->low <= h < PtRt->high)
    {
      int sgn_l, sgn_h;
      
      reduced = 1;
      
      if (l > PtRt->low)
        if (sgn_l = PtRt->poly->sgn_at(l))
          if (PtRt->sig_low < 0 && sgn_l < 0
              ||
              PtRt->sig_low > 0 && sgn_l > 0
              ||
              PtRt->sig_low == 0 &&
              PtRt->poly->num_Sturm_seq_perm(l) == PtRt->num_perm_low)
          //  Keep reduced == 1.
            PtRt->low = l;
          else
            reduced = 0;
        else  //  if l is a root of PtRt->poly
        //  Keep reduced == 1.
        {
          if (type == 1)
            type = 2;
          else  //  if (type == 3)
            type = 4;
          
          if (!--PtRt->ref_count)
            delete PtRt;
          
          PtRt = 0;
          PtBt = l;
        }
      
      //  Do the following if PtRt != 0 (and even if reduced == 0)
      //  since better high or low might be obtained.
      
      if (PtRt && h < PtRt->high)
        if (sgn_h = PtRt->poly->sgn_at(h))
          if (PtRt->sig_low < 0 && sgn_h > 0
              ||
              PtRt->sig_low > 0 && sgn_h < 0
              ||
              PtRt->sig_low == 0 &&
              PtRt->poly->num_Sturm_seq_perm(h) == PtRt->num_perm_high)
          //  Keep reduced as it is.
            PtRt->high = h;
          else
            reduced = 0;
        else  //  if h is a root of PtRt->poly
        {
          if (type == 1)
            type = 2;
          else  //  if (type == 3)
            type = 4;
          
          if (!--PtRt->ref_count)
            delete PtRt;
          
          PtRt = 0;
          PtBt = h;
          
          reduced = 1;
        }
    }
    else  //  if (h >= l) > high (>= low) or (l <=) h < low (<= high)
      reduced = 0;
  }
  else  //  if (type == 2 || type == 4)
    reduced = 1;
  
  return reduced;
}

//  int K_POINT2D :: reduce(const unsigned long num_bits_s,
//                          const unsigned long num_bits_t) const
//    reduce *this s.t.
//      [get_low_s(), get_high_s()] will contain all the numbers
//        approximated by PtRs->float_est to num_bits_s precision and
//      [get_low_t(), get_high_t()] will contain all the numbers
//        approximated by PtRt->float_est to num_bits_t precision.
//    return 1 if some reduction occurs and
//           0 otherwise.

int K_POINT2D :: reduce(const unsigned long num_bits_s,
                        const unsigned long num_bits_t) const
{
  assert(type);
  
  int reduced;

  if (type == 1)
  {
    assert(PtRs->num_roots == 1);
    assert(PtRs->ok_float);
    assert(PtRt->num_roots == 1);
    assert(PtRt->ok_float);
    
    reduced  = reduce_s(num_bits_s);
    reduced *= reduce_t(num_bits_t);
  }
  else if (type == 2)
  {
    assert(PtRs->num_roots == 1);
    assert(PtRs->ok_float);
    
    reduced = reduce_s(num_bits_s);
  }
  else if (type == 3)
  {
    assert(PtRt->num_roots == 1);
    assert(PtRt->ok_float);
    
    reduced = reduce_t(num_bits_t);
  }
  else  //  if (type == 4)
    reduced = 1;
  
  return reduced;
}

//  int K_POINT2D :: contract_s(const bigrational& tol) const
//    contract *this s.t.
//      [get_low_s(), get_high_s()] will be no larger than tol.

int K_POINT2D :: contract_s(const bigrational& tol) const
{
  assert(tol > 0);
  assert(type);
  
  if ((type == 1 || type == 2) && PtRs->high - PtRs->low > tol)
  {
    assert(PtRs->num_roots == 1);
    
    if (PtRs->ok_float)
    {
      double  d_tol;
      long    num_bits;
      
      d_tol    = tol.as_double();
      num_bits = 0;
      
      while (d_tol < 1.0 && num_bits < MAX_NUM_GOOD_FP_BITS)
      {
        d_tol *= 2.0;
        num_bits++;
      }
      
      while (num_bits > 0 && !reduce_s(num_bits))
        num_bits -= 3;
    }
    
    while ((type == 1 || type == 2) && PtRs->high - PtRs->low > tol)
      halve_s();
  }
  
  return 0;
}

//  int K_POINT2D :: contract_t(const bigrational& tol) const
//    contract *this s.t.
//      [get_low_t(), get_high_t()] will be no larger than tol.

int K_POINT2D :: contract_t(const bigrational& tol) const
{
  assert(tol > 0);
  assert(type);
  
  if ((type == 1 || type == 3) && PtRt->high - PtRt->low > tol)
  {
    assert(PtRt->num_roots == 1);
    
    if (PtRt->ok_float)
    {
      double  d_tol;
      long    num_bits;
      
      d_tol    = tol.as_double();
      num_bits = 0;
      
      while (d_tol < 1.0 && num_bits < MAX_NUM_GOOD_FP_BITS)
      {
        d_tol *= 2.0;
        num_bits++;
      }
      
      while (num_bits > 0 && !reduce_t(num_bits))
        num_bits -= 3;
    }
    
    while ((type == 1 || type == 3) && PtRt->high - PtRt->low > tol)
      halve_t();
  }
  
  return 0;
}

//  int K_POINT2D :: contract(const bigrational& tol_s,
//                            const bigrational& tol_t) const
//    contract *this s.t.
//      [get_low_s(), get_high_s()] will be no larger than tol_s and
//      [get_low_t(), get_high_t()] will be no larger than tol_t.

int K_POINT2D :: contract(const bigrational& tol_s,
                          const bigrational& tol_t) const
{
  assert(tol_s > 0);
  assert(tol_t > 0);
  assert(type);
  
  if (type == 1)
  {
    if (PtRs->high - PtRs->low > tol_s && PtRt->high - PtRt->low > tol_t)
    {
      assert(PtRs->num_roots == 1);
      assert(PtRt->num_roots == 1);
      
      if (PtRs->ok_float && PtRt->ok_float)
      {
        double  d_tol_s, d_tol_t;
        long    num_bits_s, num_bits_t;
        
        d_tol_s    = tol_s.as_double();
        num_bits_s = 0;
        
        while (d_tol_s < 1.0 && num_bits_s < MAX_NUM_GOOD_FP_BITS)
        {
          d_tol_s *= 2.0;
          num_bits_s++;
        }
        
        d_tol_t    = tol_t.as_double();
        num_bits_t = 0;
        
        while (d_tol_t < 1.0 && num_bits_t < MAX_NUM_GOOD_FP_BITS)
        {
          d_tol_t *= 2.0;
          num_bits_t++;
        }
        
        while (num_bits_s > 0 && num_bits_t > 0
                              && !reduce(num_bits_s, num_bits_t))
        {
          num_bits_s -= 3;
          num_bits_t -= 3;
        }
      }
    }
    
    if  (type == 1 && PtRs->high - PtRs->low > tol_s)
      contract_s(tol_s);
    
    if  (type == 1 && PtRt->high - PtRt->low > tol_t)
      contract_t(tol_t);
  }
  
  if (type == 2 && PtRs->high - PtRs->low > tol_s)
  {
    assert(PtRs->num_roots == 1);
    
    if (PtRs->ok_float)
    {
      double  d_tol;
      long    num_bits;
      
      d_tol    = tol_s.as_double();
      num_bits = 0;
      
      while (d_tol < 1.0 && num_bits < MAX_NUM_GOOD_FP_BITS)
      {
        d_tol *= 2.0;
        num_bits++;
      }
      
      while (num_bits > 0 && !reduce_s(num_bits))
        num_bits -= 3;
    }
    
    while (type == 2 && PtRs->high - PtRs->low > tol_s)
      halve_s();
  }
  else if (type == 3 && PtRt->high - PtRt->low > tol_t)
  {
    assert(PtRt->num_roots == 1);
    
    if (PtRt->ok_float)
    {
      double  d_tol;
      long    num_bits;
      
      d_tol    = tol_t.as_double();
      num_bits = 0;
      
      while (d_tol < 1.0 && num_bits < MAX_NUM_GOOD_FP_BITS)
      {
        d_tol *= 2.0;
        num_bits++;
      }
      
      while (num_bits > 0 && !reduce_t(num_bits))
        num_bits -= 3;
    }
    
    while (type == 3 && PtRt->high - PtRt->low > tol_t)
      halve_t();
  }
  
  return 0;
}

//  int K_POINT2D :: shrink(const bigrational& fac_s,
//                          const bigrational& fac_t) const
//    shrink *this s.t.
//      [get_low_s(), get_high_s()] will become smaller by fac_s and
//      [get_low_t(), get_high_t()] will become smaller by fac_t.

int K_POINT2D :: shrink(const bigrational& fac_s,
                        const bigrational& fac_t) const
{
  assert(type);
  
  if (type == 1)
  {
    assert(PtRs->num_roots == 1);
    assert(PtRt->num_roots == 1);
    
    contract(fac_s * (PtRs->high - PtRs->low),
             fac_t * (PtRt->high - PtRt->low));
  }
  else if (type == 2)
  {
    assert(PtRs->num_roots == 1);
    
    contract_s(fac_s * (PtRs->high - PtRs->low));
  }
  else if (type == 3)
  {
    assert(PtRt->num_roots == 1);
    
    contract_t(fac_t * (PtRt->high - PtRt->low));
  }
  
  return 0;
}

//  int K_POINT2D :: compare_s(const K_POINT2D& x) const
//    return   1 if *this >  x,
//             0 if *this == x, and
//           - 1 if *this <  x  in the s-coordinate.

int K_POINT2D :: compare_s(const K_POINT2D& x) const
{
  assert(type);
  assert(x.type);
  
  int c;
  
  if (type == 3 || type == 4)
    if (x.type == 3 || x.type == 4)
      if (PtBs > x.PtBs)
        c = 1;
      else if  (PtBs < x.PtBs)
        c = - 1;
      else  //  if (PtBs == x.PtBs)
        c = 0;
    else  //  if (x.type == 1 || x.type == 2)
      if (PtBs >= x.PtRs->high)
        c = 1;
      else if (PtBs <= x.PtRs->low)
        c = - 1;
      else  //  if (x.PtRs->low < PtBs < x.PtRs->high)
      {
        x.cut_s(PtBs);
        c = compare_s(x);
      }
  else  //  if (type == 1 || type == 2)
    if (x.type == 3 || x.type == 4)
      if (PtRs->low >= x.PtBs)
        c = 1;
      else if (PtRs->high <= x.PtBs)
        c = - 1;
      else  //  if (PtRs->low < x.PtBs < PtRs->high)
      {
        cut_s(x.PtBs);
        c = compare_s(x);
      }
    else  //  if (x.type == 1 || x.type == 2)
      if (PtRs->low > x.PtRs->low && PtRs->low < x.PtRs->high)
      {
        x.cut_s(PtRs->low);
        c = compare_s(x);
      }
      else if (PtRs->high > x.PtRs->low && PtRs->high < x.PtRs->high)
      {
        x.cut_s(PtRs->high);
        c = compare_s(x);
      }
      else if (PtRs->low < x.PtRs->low && PtRs->high > x.PtRs->low)
      {
        cut_s(x.PtRs->low);
        c = compare_s(x);
      }
      else if (PtRs->low < x.PtRs->high && PtRs->high > x.PtRs->high)
      {
        cut_s(x.PtRs->high);
        c = compare_s(x);
      }
      else if (PtRs->low >= x.PtRs->high)
        c = 1;
      else if (PtRs->high <= x.PtRs->low)
        c = - 1;
      else  //  if (PtRs->low == x.PtRs->low && PtRs->high == x.PtRs->high)
        c = 0;
  
  return c;
}

//  int K_POINT2D :: compare_t(const K_POINT2D& x) const
//    return   1 if *this >  x,
//             0 if *this == x, and
//           - 1 if *this <  x  in the t-coordinate.

int K_POINT2D :: compare_t(const K_POINT2D& x) const
{
  assert(type);
  assert(x.type);
  
  int c;
  
  if (type == 2 || type == 4)
    if (x.type == 2 || x.type == 4)
      if (PtBt > x.PtBt)
        c = 1;
      else if (PtBt < x.PtBt)
        c = - 1;
      else  //  if (PtBt == x.PtBt)
        c = 0;
    else  //  if (x.type == 1 || x.type == 3)
      if (PtBt >= x.PtRt->high)
        c = 1;
      else if (PtBt <= x.PtRt->low)
        c = - 1;
      else  //  if (x.PtRt->low < PtBt < x.PtRt->high)
      {
        x.cut_t(PtBt);
        c = compare_t(x);
      }
  else  //  if (type == 1 || type == 3)
    if (x.type == 2 || x.type == 4)
      if (PtRt->low >= x.PtBt)
        c = 1;
      else if (PtRt->high <= x.PtBt)
        c = - 1;
      else  //  if (PtRt->low < x.PtBt < PtRt->high)
      {
        cut_t(x.PtBt);
        c = compare_t(x);
      }
    else  //  if (x.type == 1 || x.type == 3)
      if (PtRt->low > x.PtRt->low && PtRt->low < x.PtRt->high)
      {
        x.cut_t(PtRt->low);
        c = compare_t(x);
      }
      else if (PtRt->high > x.PtRt->low && PtRt->high < x.PtRt->high)
      {
        x.cut_t(PtRt->high);
        c = compare_t(x);
      }
      else if (PtRt->low < x.PtRt->low && PtRt->high > x.PtRt->low)
      {
        cut_t(x.PtRt->low);
        c = compare_t(x);
      }
      else if (PtRt->low < x.PtRt->high && PtRt->high > x.PtRt->high)
      {
        cut_t(x.PtRt->high);
        c = compare_t(x);
      }
      else if (PtRt->low >= x.PtRt->high)
        c = 1;
      else if (PtRt->high <= x.PtRt->low)
        c = - 1;
      else  //  if (PtRt->low == x.PtRt->low && PtRt->high == x.PtRt->high)
        c = 0;
  
  return c;
}

////  int K_POINT2D :: compare(const K_POINT2D& x, const unsigned long i) const
////    return   1 if *this >  x,
////             0 if *this == x, and
////           - 1 if *this <  x  in the i-th coordinate.
//
//int K_POINT2D :: compare(const K_POINT2D& x, const unsigned long i) const
//{
//  assert(type);
//  assert(x.type);
//  assert(i == 0 || i == 1);
//  
//  int c;
//  
//  if (i == 0)  //  comparison in s
//    c = compare_s(x);
//  else  //  if (i == 1); comparison in t
//    c = compare_t(x);
//  
//  return c;
//}

//  int K_POINT2D :: compare_s(const bigrational& b_s) const
//    return   1 if *this >  b_s,
//             0 if *this == b_s, and
//           - 1 if *this <  b_s  in the s-coordinate.

int K_POINT2D :: compare_s(const bigrational& b_s) const
{
  assert(type);
  
  int c;
  
  if (type == 3 || type == 4)
    if (PtBs > b_s)
      c = 1;
    else if  (PtBs < b_s)
      c = - 1;
    else  //  if (PtBs == b_s)
      c = 0;
  else  //  if (typr == 1 || type == 2)
    if (PtRs->low >= b_s)
      c = 1;
    else if (PtRs->high <= b_s)
      c = - 1;
    else  //  if (PtRs->low < b_s < PtRs->high)
    {
      cut_s(b_s);
      c = compare_s(b_s);
    }
  
  return c;
}

//  int K_POINT2D :: compare_t(const bigrational& b_t) const
//    return   1 if *this >  b_t,
//             0 if *this == b_t, and
//           - 1 if *this <  b_t  in the t-coordinate.

int K_POINT2D :: compare_t(const bigrational& b_t) const
{
  assert(type);
  
  int c;
  
  if (type == 2 || type == 4)
    if (PtBt > b_t)
      c = 1;
    else if (PtBt < b_t)
      c = - 1;
    else  //  if (PtBt == b_t)
      c = 0;
  else  //  if (type == 1 || type == 3)
    if (PtRt->low >= b_t)
      c = 1;
    else if (PtRt->high <= b_t)
      c = - 1;
    else  //  if (PtRt->low < b_t < PtRt->high)
    {
      cut_t(b_t);
      c = compare_t(b_t);
    }
  
  return c;
}

////  int K_POINT2D :: compare(const bigrational& b, const unsigned long i) const
////    return   1 if *this >  b,
////             0 if *this == b, and
////           - 1 if *this <  b  in the i-th coordinate.
//
//int K_POINT2D :: compare(const bigrational& b, const unsigned long i) const
//{
//  assert(type);
//  assert(i == 0 || i == 1);
//  
//  int c;
//  
//  if (i == 0)  //  comparison in s
//    c = compare_s(b);
//  else  //  if (i == 1); comparison in t
//    c = compare_t(b);
//  
//  return c;
//}

//  int sort_s(K_POINT2D** const X, const unsigned long n)
//    insertion-sort the array X of length n in the s-coordinate.
//    return 1 if all the elements are distinct in the s-coordinate and
//           0 otherwise.

int sort_s(K_POINT2D** const X, const unsigned long n)
{
  long       i, j;
  K_POINT2D* x;
  int        c, distinct;
  
  distinct = 1;
  
  for (i = 1; i < n; i++)
  {
    x = X[i];
    
    for (j = i - 1; j >= 0 && (c = x->compare_s(*X[j])) < 0; j--)
      X[j + 1] = X[j];
    
    X[j + 1] = x;
    
    if (!c)
      distinct = 0;
  }
  
  return distinct;
}

//  int sort_t(K_POINT2D** const X, const unsigned long n)
//    insertion-sort the array X of length n in the t-coordinate.
//    return 1 if all the elements are distinct in the t-coordinate and
//           0 otherwise.

int sort_t(K_POINT2D** const X, const unsigned long n)
{
  long       i, j;
  K_POINT2D* x;
  int        c, distinct;
  
  distinct = 1;
  
  for (i = 1; i < n; i++)
  {
    x = X[i];
    
    for (j = i - 1; j >= 0 && (c = x->compare_t(*X[j])) < 0; j--)
      X[j + 1] = X[j];
    
    X[j + 1] = x;
    
    if (!c)
      distinct = 0;
  }
  
  return distinct;
}

////  int sort(K_POINT2D** const X, const unsigned long n, const unsigned long i)
////    insertion-sort the array X of length n in the i-th coordinate.
////    return 1 if all the elements are distinct in the i-th coordinate and
////           0 otherwise.
//
//int sort(K_POINT2D** const X, const unsigned long n, const unsigned long i)
//{
//  assert(i == 0 || i == 1);
//  
//  int distinct;
//  
//  if (i == 0)
//    distinct = sort_s(X, n);
//  else  //  if (i == 1)
//    distinct = sort_t(X, n);
//  
//  return distinct;
//}

//  unsigned long get_pts(const bigrational& l_s, const bigrational& h_s,
//                        const bigrational& l_t, const bigrational& h_t,
//                        const K_RATPOLY& P1, const K_RATPOLY& P2,
//                        K_POINT2D**& pts,
//                        const bigrational& tol, const int count_edges)
//    computes the intersections "pts" of 2 bivariate polynomials P1 and P2
//                              on/in the region [l_s, h_s] x [l_t, h_t].
//    returns the number of intersections.
//    if tol > 0 then
//      a box for each intersection is no larger than tol in any direction.
//    if tol = 0 then
//      there is no limit on the size of boxes for intersections.
//    if count_edges = 1 then
//      the intersections on the boundary edges of the region are counted.
//    if count_edges = 0 then
//      the intersections on the boundary edges of the region are ignored.

unsigned long get_pts(const bigrational& l_s, const bigrational& h_s,
                      const bigrational& l_t, const bigrational& h_t,
                      const K_RATPOLY& P1, const K_RATPOLY& P2,
                      K_POINT2D**& pts,
                      const bigrational& tol, const int count_edges)
{
  assert(P1.num_vars == 2);
  assert(P2.num_vars == 2);
  
#ifdef _EXPERIMENT
  num_kpoint2d_get_pts++;
#endif
  
  unsigned long num_pts;
  
  if (!P1.deg[0] && !P1.deg[1] || !P2.deg[0] && !P2.deg[1])
  //  0.1.  when P1 or P2 is a constant polynomial
  {
    pts     = 0;
    num_pts = 0;
  }
  else if (P1.eq_upto_const(P2))
  //  0.2.  when P1 and P2 are the same (up to constant multiple)
  {
    pts     = 0;
    num_pts = 0;
  }
  else if (l_s > h_s || l_t > h_t)
  //  0.3.  when the region [l_s, l_t] x [h_s, h_t] is not well-defined
  {
    pts     = 0;
    num_pts = 0;
  }
  else if (l_s == h_s && l_t == h_t)
  //  1.  when the region [l_s, l_t] x [h_s, h_t] is
  //        shrunken to the point (l_s, l_t).
    if (count_edges)
    {
      bigrational b[2];
      
      b[0] = l_s;
      b[1] = l_t;
      
      if (!P1.sgn_at(b) && !P2.sgn_at(b))
      //  If (l_s, l_t) is a common root of P1 and P2 then
      //    pts consists of a single point, namely (l_s, l_t).
      {
        pts    = new K_POINT2D* [num_pts = 1];
        pts[0] = new K_POINT2D(l_s, h_s, P1, P2);
        pts[0]->ref_count++;
      }
      else  //  if (P1.sgn_at(b) || P2.sgn_at(b))
      {
        pts     = 0;
        num_pts = 0;
      }
    }
    else  //  if (!count_edges)
    {
      pts     = 0;
      num_pts = 0;
    }
  else if (l_s == h_s && l_t < h_t)
  //  2.  when the region [l_s, l_t] x [h_s, h_t] is
  //        shrunken to the vertical line segment [(l_s, l_t), (l_s, h_t)]
    if (count_edges)
    //  See if 2 univariate polynomials P1 | s = l_s and P2 | s = l_s
    //        have common roots in the interval [l_t, h_t].
      num_pts = get_pts_proto(l_s,
                              l_t, h_t,
                              gcd(P1.subst_val(0, l_s), P2.subst_val(0, l_s)),
                              P1, P2,
                              pts,
                              tol, 1);
    else  //  if (!count_edges)
    {
      pts     = 0;
      num_pts = 0;
    }
  else if (l_s < h_s && l_t == h_t)
  //  3.  when the region [l_s, l_t] x [h_s, h_t] is
  //        shrunken to the horizontal line segment [(l_s, l_t), (h_s, l_t)]
    if (count_edges)
    //  See if 2 univariate polynomials P1 | t = l_t and P2 | t = l_t
    //        have common roots on the interval [l_s, h_s].
      num_pts = get_pts_proto(l_s, h_s,
                              gcd(P1.subst_val(1, l_t), P2.subst_val(1, l_t)),
                              l_t,
                              P1, P2,
                              pts,
                              tol, 1);
    else  //  if (!count_edges)
    {
      pts     = 0;
      num_pts = 0;
    }
  else if (P1.deg[0] == 1 && !P1.deg[1])
  //  4.  when P1 = 0 is a vertical line s = b1_s
  {
    bigrational b1_s;
    
    b1_s = - P1.coeffs[1] / P1.coeffs[0];
    //  P1.deg[0] == 1 => P1.coeffs[0] != 0
    
    if (b1_s > l_s && b1_s < h_s
        ||
        count_edges && b1_s >= l_s && b1_s <= h_s)
    //  If b1_s is in/on the interval [l_s, h_s] then
    //    see if a univariate polynomial P2 | s = b1_s
    //          has roots in/on the interval [l_t, h_t].
      num_pts = get_pts_proto(b1_s,
                              l_t, h_t, P2.subst_val(0, b1_s),
                              P1, P2,
                              pts,
                              tol, count_edges);
    else  //  if b1_s is not in/on the interval [l_s, h_s]
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else if (!P1.deg[0] && P1.deg[1] == 1)
  //  5.  when P1 is a horizontal line t = b1_t.
  {
    bigrational b1_t;
    
    b1_t = - P1.coeffs[1] / P1.coeffs[0];
    //  P1.deg[1] == 1 => P1.coeffs[0] != 0
    
    if (b1_t > l_t && b1_t < h_t
        ||
        count_edges && b1_t >= l_t && b1_t <= h_t)
    //  If b1_t is in/on the interval [l_t, h_t] then
    //    see if a univariate polynomial P2 | t = b1_t
    //          has roots in/on the interval [l_s, h_s].
      num_pts = get_pts_proto(l_s, h_s, P2.subst_val(1, b1_t),
                              b1_t,
                              P1, P2,
                              pts,
                              tol, count_edges);
    else  //  if b1_t is not in/on the interval [l_t, h_t]
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else if (P2.deg[0] == 1 && !P2.deg[1])
  //  6.  when P2 is a vertical line s = b2_s
  {
    bigrational b2_s;
    
    b2_s = - P2.coeffs[1] / P2.coeffs[0];
    //  P2.deg[0] == 1 => P2.coeffs[0] != 0
    
    if (b2_s > l_s && b2_s < h_s
        ||
        count_edges && b2_s >= l_s && b2_s <= h_s)
    //  If b2_s is in/on the interval [l_s, h_s] then
    //    see if a univariate polynomial P2 | s = b2_s
    //          has roots in/on the interval [l_t, h_t].
      num_pts = get_pts_proto(b2_s,
                              l_t, h_t, P1.subst_val(0, b2_s),
                              P1, P2,
                              pts,
                              tol, count_edges);
    else  //  if b2_s is not in/on the interval [l_s, h_s]
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else if (!P2.deg[0] && P2.deg[1] == 1)
  //  7.  when P2 is a horizontal line t = b2_t
  {
    bigrational b2_t;
    
    b2_t = - P2.coeffs[1] / P2.coeffs[0];
    //  P2.deg[1] == 1 => P2.coeffs[0] != 0
    
    if (b2_t > l_t && b2_t < h_t
        ||
        count_edges && b2_t >= l_t && b2_t <= h_t)
    //  If b2_t is in/on the interval [l_t, h_t] then
    //    see if a univariate polynomial P2 | t = b2_t
    //          has roots in/on the interval [l_s, h_s].
      num_pts = get_pts_proto(l_s, h_s, P1.subst_val(1, b2_t),
                              b2_t,
                              P1, P2,
                              pts,
                              tol, count_edges);
    else  //  if b2_t is not in/on the interval [l_t, h_t]
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else if (P1.get_total_deg() == 1 && P2.get_total_deg() == 1)
  //  8.  When P1 and P2 are linear, say,
  //        P1 = c11 s + c12 t + c13 and P2 = c21 s + c22 t + c23,
  //      use Cramer's formula to solve the system P1 = P2 = 0 and
  //      see if it has a solution in/on the region [l_s, h_s] x [l_t, h_t].
  {
    bigrational c11, c12, c13, c21, c22, c23, d;
    
    c11 = P1.coeffs[1];
    c12 = P1.coeffs[2];
    c13 = P1.coeffs[3];
    c21 = P2.coeffs[1];
    c22 = P2.coeffs[2];
    c23 = P2.coeffs[3];
    d   = c11 * c22 - c12 * c21;
    
    if (sgn(d))
    {
      bigrational b_s, b_t;
      
      b_s = (- c13 * c22 + c12 * c23) / d;
      b_t = (- c11 * c23 + c13 * c21) / d;
      
      if (b_s > l_s && b_s < h_s && b_t > l_t && b_t < h_t
          ||
          count_edges && b_s >= l_s && b_s <= h_s && b_t >= l_t && b_t <= h_t)
      {
        pts    = new K_POINT2D* [num_pts = 1];
        pts[0] = new K_POINT2D(b_s, b_t, P1, P2);
        pts[0]->ref_count++;
      }
      else
      //  if the solution (b_s, b_t) is not
      //    in/on the region [l_s, h_s] x [l_t, h_t]
      {
        pts     = 0;
        num_pts = 0;
      }
    }
    else  //  if (d == 0), i.e., P1 = 0 and P2 = 0 are parallel
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else
  //  9.  Otherwise, perform 2-D root finding.
  //      We look the following objects:
  //
  //        (l_s, h_t)  edge h_t  (h_s, h_t)
  //                 \            /
  //                  *----------*
  //                  |          |
  //            edge  |          |  edge
  //             l_s  |          |   h_s
  //                  |          |
  //                  |          |
  //                  *----------*
  //                 /            \
  //        (l_s, l_t)  edge l_t  (h_s, l_t)
  //
  //
  //      CAUTION: Corners and boundary edges are looked at
  //               only if count_edges == 1.
  {
    unsigned long i, j;
    K_POINT2D     c_l_l, c_l_h, c_h_h, c_h_l;
    unsigned long num_c_l_l, num_c_l_h, num_c_h_h, num_c_h_l;
    K_POINT2D**   e_l_t;
    K_POINT2D**   e_l_s;
    K_POINT2D**   e_h_t;
    K_POINT2D**   e_h_s;
    unsigned long num_e_l_t, num_e_l_s, num_e_h_t, num_e_h_s;
    K_POINT2D**   interior;
    unsigned long num_interior;
    
    if (count_edges)
    {
      //  9.1.  See if corners are intersections of P1 and P2.
      
      bigrational b[2];
      
      //  9.1.1.  See if (l_s, l_t) is an intersection of P1 and P2.
      
      b[0] = l_s;
      b[1] = l_t;
      
      if (!P1.sgn_at(b) && !P2.sgn_at(b))
        num_c_l_l = 1;
      else  //  if (P1.sgn_at(b) || P2.sgn_at(b))
        num_c_l_l = 0;
      
      //  9.1.2.  See if (l_s, h_t) is an intersection of P1 and P2.
      
      b[0] = l_s;
      b[1] = h_t;
      
      if (!P1.sgn_at(b) && !P2.sgn_at(b))
        num_c_l_h = 1;
      else  //  if (P1.sgn_at(b) || P2.sgn_at(b))
        num_c_l_h = 0;
      
      //  9.1.3.  See if (h_s, h_t) is an intersection of P1 and P2.
      
      b[0] = h_s;
      b[1] = h_t;
      
      if (!P1.sgn_at(b) && !P2.sgn_at(b))
        num_c_h_h = 1;
      else  //  if (P1.sgn_at(b) || P2.sgn_at(b))
        num_c_h_h = 0;
      
      //  9.1.4.  See if (h_s, l_t) is an intersection of P1 and P2.
      
      b[0] = h_s;
      b[1] = l_t;
      
      if (!P1.sgn_at(b) && !P2.sgn_at(b))
        num_c_h_l = 1;
      else  //  if (P1.sgn_at(b) || P2.sgn_at(b))
        num_c_h_l = 0;
      
      //  9.2.  See if there are intersections of P1 and P2 on boundary edges.
      
      //  9.2.1.  See if there are intersections of P1 and P2 on edge l_t.
      
      num_e_l_t = get_pts_proto(l_s, h_s,
                                gcd(P1.subst_val(1, l_t),
                                    P2.subst_val(1, l_t)),
                                l_t,
                                P1, P2,
                                e_l_t,
                                tol, 0);
      
      //  9.2.2.  See if there are intersections of P1 and P2 on edge l_s.
      
      num_e_l_s = get_pts_proto(l_s,
                                l_t, h_t,
                                gcd(P1.subst_val(0, l_s),
                                    P2.subst_val(0, l_s)),
                                P1, P2,
                                e_l_s,
                                tol, 0);
      
      //  9.2.3.  See if there are intersections of P1 and P2 on edge h_t.
      
      num_e_h_t = get_pts_proto(l_s, h_s,
                                gcd(P1.subst_val(1, h_t),
                                    P2.subst_val(1, h_t)),
                                h_t,
                                P1, P2,
                                e_h_t,
                                tol, 0);
      
      //  9.2.4.  See if there are intersections of P1 and P2 on edge h_s.
      
      num_e_h_s = get_pts_proto(h_s,
                                l_t, h_t,
                                gcd(P1.subst_val(0, h_s),
                                    P2.subst_val(0, h_s)),
                                P1, P2,
                                e_h_s,
                                tol, 0);
    }
    else  //  if (!count_edegs)
      num_c_l_l = num_c_l_h = num_c_h_h = num_c_h_l =
      num_e_l_t = num_e_l_s = num_e_h_t = num_e_h_s = 0;
    
    //  9.3.  See if there are intersections of P1 and P2
    //              in the open region (l_s, l_t) x (h_s, h_t).
    
    num_interior = get_pts_interior(l_s, l_t, h_s, h_t, P1, P2, interior, tol);
    
    //  9.4.  Gather all together.
    
    if ((num_pts = num_c_l_l + num_c_l_h + num_c_h_h + num_c_h_l +
                   num_e_l_t + num_e_l_s + num_e_h_t + num_e_h_s +
                   num_interior) > 0)
    {
      pts = new K_POINT2D* [num_pts];
      i   = 0;
      
      if (num_c_l_l > 0)
      {
        pts[i] = new K_POINT2D(l_s, l_t, P1, P2);
        pts[i++]->ref_count++;
      }
      
      if (num_c_l_h > 0)
      {
        pts[i] = new K_POINT2D(l_s, h_t, P1, P2);
        pts[i++]->ref_count++;
      }
      
      if (num_c_h_h > 0)
      {
        pts[i] = new K_POINT2D(h_s, h_t, P1, P2);
        pts[i++]->ref_count++;
      }
      
      if (num_c_h_l > 0)
      {
        pts[i] = new K_POINT2D(h_s, l_t, P1, P2);
        pts[i++]->ref_count++;
      }
      
      for (j = 0; j < num_e_l_t; j++)
      {
        pts[i] = e_l_t[j];
        pts[i++]->ref_count++;
      }
      
      for (j = 0; j < num_e_l_s; j++)
      {
        pts[i] = e_l_s[j];
        pts[i++]->ref_count++;
      }
      
      for (j = 0; j < num_e_h_t; j++)
      {
        pts[i] = e_h_t[j];
        pts[i++]->ref_count++;
      }
      
      for (j = 0; j < num_e_h_s; j++)
      {
        pts[i] = e_h_s[j];
        pts[i++]->ref_count++;
      }
      
      for (j = 0; j < num_interior; j++)
      {
        pts[i] = interior[j];
        pts[i++]->ref_count++;
      }
      
      assert(i == num_pts);
      
      for (j = 0; j < num_e_l_t; j++)
        if (!--e_l_t[j]->ref_count)
          delete e_l_t[j];
      
      if (num_e_l_t > 0)
        delete [] e_l_t;
      
      for (j = 0; j < num_e_l_s; j++)
        if (!--e_l_s[j]->ref_count)
          delete e_l_s[j];
      
      if (num_e_l_s > 0)
        delete [] e_l_s;
      
      for (j = 0; j < num_e_h_t; j++)
        if (!--e_h_t[j]->ref_count)
          delete e_h_t[j];
      
      if (num_e_h_t > 0)
        delete [] e_h_t;
      
      for (j = 0; j < num_e_h_s; j++)
        if (!--e_h_s[j]->ref_count)
          delete e_h_s[j];
      
      if (num_e_h_s > 0)
        delete [] e_h_s;
      
      for (j = 0; j < num_interior; j++)
        if (!--interior[j]->ref_count)
          delete interior[j];
      
      if (num_interior > 0)
        delete [] interior;
    }
    else  //  if (num_pts == 0)
      pts = 0;
  }
  
  return num_pts;
}

#define MAX_NUM_PTS_HIT 16

//  unsigned long get_pts_interior(const bigrational& l_s,
//                                 const bigrational& l_t,
//                                 const bigrational& h_s,
//                                 const bigrational& h_t,
//                                 const K_RATPOLY& P1, const K_RATPOLY& P2,
//                                 K_POINT2D**& pts,
//                                 const bigrational& tol)
//    computes the intersections "pts" of 2 bivariate polynomials P1 and P2
//                              in the open region (l_s, h_s) x (l_t, h_t)
//    by using box-hit analysis.
//    returns the number of intersections.
//    if tol > 0 then
//      a box for each intersection is no larger than tol in any direction.
//    if tol = 0 then
//      there is no limit on the size of boxes for intersections.

unsigned long get_pts_interior(const bigrational& l_s, const bigrational& l_t,
                               const bigrational& h_s, const bigrational& h_t,
                               const K_RATPOLY& P1, const K_RATPOLY& P2,
                               K_POINT2D**& pts,
                               const bigrational& tol)
{
  assert(P1.num_vars == 2);
  assert(P2.num_vars == 2);
  
  unsigned long num_pts;
  
  K_RATPOLY s_only = GoodSylvester(P1, P2, 1);
  K_RATPOLY t_only = GoodSylvester(P1, P2, 0);
  
  //  1.  Compute all the roots pts_s and pts_t of s_only and t_only
  //        that are supersets of the s- and t- coordinate
  //                                 of the intersections pts of P1 and P2.
  
  unsigned long num_pts_s, num_pts_t;
  K_POINT1D**   pts_s;
  K_POINT1D**   pts_t;
  
  if (!(num_pts_s = get_pts(l_s, h_s, s_only, pts_s, tol, 0)))
  {
    pts     = 0;
    num_pts = 0;
  }
  else if (!(num_pts_t = get_pts(l_t, h_t, t_only, pts_t, tol, 0)))
  {
    pts     = 0;
    num_pts = 0;
  }
  else  //  if (num_pts_s > 0 && num_pts_t > 0)
  {
    long        i, j, k, l;
    bigrational b[2];
    
    //  2.  Divide pts_s into pts_s1 and pts_s2
    //        that are the sets of type 1 and 2 points.
    //      Also, divide pts_t into pts_t1 and pts_t2
    //              that are the sets of type 1 and 2 points.
    
    unsigned long num_pts_s1, num_pts_s2, num_pts_t1, num_pts_t2;
    K_POINT1D**   pts_s1;
    K_POINT1D**   pts_s2;
    K_POINT1D**   pts_t1;
    K_POINT1D**   pts_t2;
    K_POINT2D**   pts_proto;
    
    pts_s1     = new K_POINT1D* [num_pts_s];
    pts_s2     = new K_POINT1D* [num_pts_s];
    num_pts_s1 = num_pts_s2 = 0;
    
    for (i = 0; i < num_pts_s; i++)
      if (pts_s[i]->type == 1)
      {
        pts_s1[num_pts_s1] = pts_s[i];
        pts_s1[num_pts_s1++]->ref_count++;
      }
      else  //  if (pts_s[i]->type == 2)
      {
        pts_s2[num_pts_s2] = pts_s[i];
        pts_s2[num_pts_s2++]->ref_count++;
      }
    
    assert(num_pts_s == num_pts_s1 + num_pts_s2);    
    assert(sort(pts_s1, num_pts_s1));  //  Make sure that the points in pts_s1
                                       //                   are separable.
    
    for (i = 0; i < num_pts_s; i++)
      if (!--pts_s[i]->ref_count)
        delete pts_s[i];
    
    delete [] pts_s;  //  num_pts_s > 0 => pts_s != 0
    
    pts_t1 = new K_POINT1D* [num_pts_t];
    pts_t2 = new K_POINT1D* [num_pts_t];
    num_pts_t1 = num_pts_t2 = 0;
    
    for (i = 0; i < num_pts_t; i++)
      if (pts_t[i]->type == 1)
      {
        pts_t1[num_pts_t1] = pts_t[i];
        pts_t1[num_pts_t1++]->ref_count++;
      }
      else  //  if (pts_t[i]->type == 2)
      {
        pts_t2[num_pts_t2] = pts_t[i];
        pts_t2[num_pts_t2++]->ref_count++;
      }
    
    assert(num_pts_t == num_pts_t1 + num_pts_t2);
    assert(sort(pts_t1, num_pts_t1));  //  Make sure that the points in pts_t1
                                       //                   are separable.
    
    for (i = 0; i < num_pts_t; i++)
      if (!--pts_t[i]->ref_count)
        delete pts_t[i];
    
    delete [] pts_t;  //  num_pts_t > 0 => pts_t != 0
    
    //  3.  See if a point in pts_s x pts_t
    //            is actually an intersection of P1 and P2.
    //      Let num_pts_all to be
    //        the number of possible intersections of P1 and P2.
    //      Then, num_pts_all =
    //              min { num_pts_s * num_pts_t
    //                    deg_s P1 * deg_t P2 + deg_t P1 * deg_s P2 }.
    
    unsigned long num_pts_all, num_pts_st;
    
    num_pts_all = P1.deg[0] * P2.deg[1] + P1.deg[1] * P2.deg[0];
    num_pts_st  = num_pts_s * num_pts_t;
    
    if (num_pts_all > num_pts_st)
      num_pts_all = num_pts_st;
    
    pts_proto = new K_POINT2D* [num_pts_all];
    num_pts   = 0;
    
    //  3-1.  type 2 v.s. type 2
    
    for (i = 0; i < num_pts_s2; i++)
    {
      b[0] = pts_s2[i]->PtB;
      
      for (j = 0; j < num_pts_t2; j++)
      {
        b[1] = pts_t2[j]->PtB;
        
        if (!P1.sgn_at(b) && !P2.sgn_at(b))
        {
          pts_proto[num_pts] = new K_POINT2D(b[0], b[1], P1, P2);
          pts_proto[num_pts++]->ref_count++;
        }
      }
    }
    
    //  3-2.  type 2 v.s. type 1
    
    bigrational   b_s;
    K_RATPOLY     G21;
    unsigned long num_pts_21;
    K_POINT1D**   pts_21;
    
    if (num_pts_t1 > 0)
      for (i = 0; i < num_pts_s2; i++)
      {
        b_s = pts_s2[i]->PtB;
        G21 = gcd(P1.subst_val(0, b_s), P2.subst_val(0, b_s));
        
        //  (b_s, y) is an intersection of P1 and P2            <=>
        //  y is a common root of P1 | s = b_s and P2 | s = b_s <=>
        //  y is a root of G21 = gcd(P1 | s = b_s, P2 | s = b_s).
        
        if (G21.deg[0] > 0)
        {
          num_pts_21 = get_pts(pts_t1[0]->get_low(),
                               pts_t1[num_pts_t1 - 1]->get_high(),
                               G21,
                               pts_21,
                               tol, 0);
          
          if (num_pts_21 > 0)
          {
            assert(sort(pts_21, num_pts_21));  //  Make sure that
                                               //    the points in pts_21
                                               //      are separable.
            
            for (j = 0; j < num_pts_21; j++)
            {
              pts_proto[num_pts] = new K_POINT2D(b_s, *pts_21[j], P1, P2);
              pts_proto[num_pts++]->ref_count++;
              
              if (!--pts_21[j]->ref_count)
                delete pts_21[j];
            }
            
            delete [] pts_21;
          }
        }
      }
    
    //  3-3.  type 1 v.s. type 2
    
    bigrational   b_t;
    K_RATPOLY     G12;
    unsigned long num_pts_12;
    K_POINT1D**   pts_12;
    
    if (num_pts_s1 > 0)
      for (i = 0; i < num_pts_t2; i++)
      {
        b_t = pts_t2[i]->PtB;
        G12 = gcd(P1.subst_val(1, b_t), P2.subst_val(1, b_t));
        
        //  (x, b_t) is an intersection of P1 and P2            <=>
        //  x is a common root of P1 | t = b_t and P2 | t = b_t <=>
        //  x is a root of G12 = gcd(P1 | t = b_t, P2 | t = b_t).
        
        if (G12.deg[0] > 0)
        {
          num_pts_12 = get_pts(pts_s1[0]->get_low(),
                               pts_s1[num_pts_s1 - 1]->get_high(),
                               G12,
                               pts_12,
                               tol, 0);
          
          if (num_pts_12 > 0)
          {
            assert(sort(pts_12, num_pts_12));  //  Make sure that
                                               //    the points in pts_12
                                               //      are separable.
            
            for (j = 0; j < num_pts_12; j++)
            {
              pts_proto[num_pts] = new K_POINT2D(*pts_12[j], b_t, P1, P2);
              pts_proto[num_pts++]->ref_count++;
              
              if (!--pts_12[j]->ref_count)
                delete pts_12[j];
            }
            
            delete [] pts_12;
          }
        }
      }
    
    //  3-4.  type 1 v.s. type 1
    //        Perform box hit analysis.
    //        For each side of the box for a point in pts_s1 x pts_t1,
    //          count the number of intersections with P1 and P2.
    
    bigrational*     s_l;
    bigrational*     s_h;
    bigrational*     t_l;
    bigrational*     t_h;
    unsigned long*** boxhits;
    unsigned long**  boxhits1;
    unsigned long**  boxhits2;
    
    if (num_pts_s1 > 0 && num_pts_t1 > 0)
    {
      s_l = new bigrational [num_pts_s1];
      s_h = new bigrational [num_pts_s1];
      t_l = new bigrational [num_pts_t1];
      t_h = new bigrational [num_pts_t1];
      
      for (i = 0; i < num_pts_s1; i++)
      {
        s_l[i] = pts_s1[i]->get_low();
        s_h[i] = pts_s1[i]->get_high();
      }
      
      for (i = 0; i < num_pts_t1; i++)
      {
        t_l[i] = pts_t1[i]->get_low();
        t_h[i] = pts_t1[i]->get_high();
      }
      
      boxhits  = new unsigned long** [num_pts_s1];
      boxhits1 = new unsigned long* [num_pts_s1];
      boxhits2 = new unsigned long* [num_pts_s1];
      
      for (i = 0; i < num_pts_s1; i++)
      {
        boxhits[i]  = new unsigned long* [num_pts_t1];
        boxhits1[i] = new unsigned long [num_pts_t1];
        boxhits2[i] = new unsigned long [num_pts_t1];
        
        for (j = 0; j < num_pts_t1; j++)
        {
          boxhits[i][j]  = new unsigned long [MAX_NUM_PTS_HIT];
          boxhits1[i][j] = 0;
          boxhits2[i][j] = 0;
          
          for (k = 0; k < MAX_NUM_PTS_HIT; k++)
            boxhits[i][j][k] = 0;
        }
      }
      
      K_RATPOLY     Q1, Q2;
      unsigned long num_pts1, num_pts2;
      K_POINT1D**   pts1;
      K_POINT1D**   pts2;
      K_RATPOLY*    poly_hit;
      long*         match1;
      long*         match2;
      int           sgn_dtds;
      unsigned long count1, count2;
      unsigned long num_pts_hit;
      K_POINT1D**   pts_hit;
      
      //  low_s: the left side of a box
      
      for (i = 0; i < num_pts_s1; i++)
      {
        Q1       = P1.subst_val(0, s_l[i]);
        num_pts1 = get_pts(t_l[0], t_h[num_pts_t1 - 1], Q1, pts1, init_tol, 0);
        
        if (num_pts1 > 0)
        {
          poly_hit = pts1[0]->poly;
          sort(pts1, num_pts1);
          match_intervals(t_l, t_h, num_pts_t1, pts1, match1, num_pts1);
        }
        else  //  if (num_pts1 == 0)
          poly_hit = 0;
        
        Q2       = P2.subst_val(0, s_l[i]);
        num_pts2 = get_pts(t_l[0], t_h[num_pts_t1 - 1], Q2, pts2, init_tol, 0);
        
        if (num_pts2 > 0)
        {
          sort(pts2, num_pts2);
          match_intervals(t_l, t_h, num_pts_t1, pts2, match2, num_pts2);
        }
        
        count1 = count2 = 0;
        
        for (j = 0; j < num_pts_t1; j++)
        {
          if (!Q1.sgn_at(t_l[j]))
          {
            b[0]     = s_l[i];
            b[1]     = t_l[j];
            sgn_dtds = - sgn(P1.derivative(0).evaluate(b)) *
                         sgn(P1.derivative(1).evaluate(b));
            
            if (sgn_dtds > 0)
            {
              boxhits[i][j][boxhits1[i][j] + boxhits2[i][j]] = 1;
              boxhits1[i][j]++;
            }
            else if (!sgn_dtds)
              boxhits1[i][j] += 5;
          }
          
          if (!Q2.sgn_at(t_l[j]))
          {
            b[0]     = s_l[i];
            b[1]     = t_l[j];
            sgn_dtds = - sgn(P2.derivative(0).evaluate(b)) *
                         sgn(P2.derivative(1).evaluate(b));
            
            if (sgn_dtds > 0)
            {
              boxhits[i][j][boxhits1[i][j] + boxhits2[i][j]] = 2;
              boxhits2[i][j]++;
            }
            else if (!sgn_dtds)
              boxhits2[i][j] += 5;
          }
          
          pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
          num_pts_hit = 0;
          
          while (count1 < num_pts1 && match1[count1] < j)
            count1++;
          
          while (count1 < num_pts1 && match1[count1] == j)
          {
            boxhits1[i][j]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts1[count1++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          while (count2 < num_pts2 && match2[count2] < j)
            count2++;
          
          while (count2 < num_pts2 && match2[count2] == j)
          {
            boxhits2[i][j]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts2[count2++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          if (!sort(pts_hit, num_pts_hit))
          {
            for (k = 0; k < num_pts_hit - 1; k++)
              for (l = k + 1; l < num_pts_hit; l++)
                pts_hit[k]->separate(*pts_hit[l]);
            
            assert(sort(pts_hit, num_pts_hit));
          }
          
          for (k = 0; k < num_pts_hit; k++)
          {
            l = boxhits1[i][j] + boxhits2[i][j] - num_pts_hit + k;
            
            if (poly_hit && *pts_hit[k]->poly == *poly_hit)
              boxhits[i][j][l] = 1;
            else  //  if (!poly_hit || *pts_hit[k]->poly != *poly_hit)
              boxhits[i][j][l] = 2;
            
            if (!--pts_hit[k]->ref_count)
              delete pts_hit[k];
          }
          
          delete [] pts_hit;
        }
        
        for (j = 0; j < num_pts1; j++)
          if (!--pts1[j]->ref_count)
            delete pts1[j];
        
        if (num_pts1 > 0)
        {
          delete [] pts1;
          delete [] match1;
        }
        
        for (j = 0; j < num_pts2; j++)
          if (!--pts2[j]->ref_count)
            delete pts2[j];
        
        if (num_pts2 > 0)
        {
          delete [] pts2;
          delete [] match2;
        }
      }
      
      //  high_t: the top side of a box
      
      for (i = 0; i < num_pts_t1; i++)
      {
        Q1       = P1.subst_val(1, t_h[i]);
        num_pts1 = get_pts(s_l[0], s_h[num_pts_s1 - 1], Q1, pts1, init_tol, 0);
        
        if (num_pts1 > 0)
        {
          poly_hit = pts1[0]->poly;
          sort(pts1, num_pts1);
          match_intervals(s_l, s_h, num_pts_s1, pts1, match1, num_pts1);
        }
        else  //  if (num_pts1 == 0)
          poly_hit = 0;
        
        Q2       = P2.subst_val(1, t_h[i]);
        num_pts2 = get_pts(s_l[0], s_h[num_pts_s1 - 1], Q2, pts2, init_tol, 0);
        
        if (num_pts2 > 0)
        {
          sort(pts2, num_pts2);
          match_intervals(s_l, s_h, num_pts_s1, pts2, match2, num_pts2);
        }
        
        count1 = count2 = 0;
        
        for (j = 0; j < num_pts_s1; j++)
        {
          if (!Q1.sgn_at(s_l[j]))
          {
            b[0]     = s_l[j];
            b[1]     = t_h[i];
            sgn_dtds = - sgn(P1.derivative(0).evaluate(b)) *
                         sgn(P1.derivative(1).evaluate(b));
            
            if (sgn_dtds < 0)
            {
              boxhits[j][i][boxhits1[j][i] + boxhits2[j][i]] = 1;
              boxhits1[j][i]++;
            }
            else if (!sgn_dtds)
              boxhits1[j][i] += 5;
          }
          
          if (!Q2.sgn_at(s_l[j]))
          {
            b[0]     = s_l[j];
            b[1]     = t_h[i];
            sgn_dtds = - sgn(P2.derivative(0).evaluate(b)) *
                         sgn(P2.derivative(1).evaluate(b));
            
            if (sgn_dtds < 0)
            {
              boxhits[j][i][boxhits1[j][i] + boxhits2[j][i]] = 2;
              boxhits2[j][i]++;
            }
            else if (!sgn_dtds)
              boxhits2[j][i] += 5;
          }
          
          pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
          num_pts_hit = 0;
          
          while (count1 < num_pts1 && match1[count1] < j)
            count1++;
          
          while (count1 < num_pts1 && match1[count1] == j)
          {
            boxhits1[j][i]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts1[count1++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          while (count2 < num_pts2 && match2[count2] < j)
            count2++;
          
          while (count2 < num_pts2 && match2[count2] == j)
          {
            boxhits2[j][i]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts2[count2++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          if (!sort(pts_hit, num_pts_hit))
          {
            for (k = 0; k < num_pts_hit - 1; k++)
              for (l = k + 1; l < num_pts_hit; l++)
                pts_hit[k]->separate(*pts_hit[l]);
            
            assert(sort(pts_hit, num_pts_hit));
          }
          
          for (k = 0; k < num_pts_hit; k++)
          {
            l = boxhits1[j][i] + boxhits2[j][i] - num_pts_hit + k;
            
            if (poly_hit && *pts_hit[k]->poly == *poly_hit)
              boxhits[j][i][l] = 1;
            else  //  if (!poly_hit || *pts_hit[k]->poly != *poly_hit)
              boxhits[j][i][l] = 2;
            
            if (!--pts_hit[k]->ref_count)
              delete pts_hit[k];
          }
          
          delete [] pts_hit;
        }
        
        for (j = 0; j < num_pts1; j++)
          if (!--pts1[j]->ref_count)
            delete pts1[j];
        
        if (num_pts1 > 0)
        {
          delete [] pts1;
          delete [] match1;
        }
        
        for (j = 0; j < num_pts2; j++)
          if (!--pts2[j]->ref_count)
            delete pts2[j];
        
        if (num_pts2 > 0)
        {
          delete [] pts2;
          delete [] match2;
        }
      }
      
      //  high_s: the right side of a box: (counted in reverse order)
      
      for (i = 0; i < num_pts_s1; i++)
      {
        Q1       = P1.subst_val(0, s_h[i]);
        num_pts1 = get_pts(t_l[0], t_h[num_pts_t1 - 1], Q1, pts1, init_tol, 0);
        
        if (num_pts1 > 0)
        {
          poly_hit = pts1[0]->poly;
          sort(pts1, num_pts1);
          match_intervals(t_l, t_h, num_pts_t1, pts1, match1, num_pts1);
        }
        else  //  if (num_pts1 == 0)
          poly_hit = 0;
        
        Q2       = P2.subst_val(0, s_h[i]);
        num_pts2 = get_pts(t_l[0], t_h[num_pts_t1 - 1], Q2, pts2, init_tol, 0);
        
        if (num_pts2 > 0)
        {
          sort(pts2, num_pts2);
          match_intervals(t_l, t_h, num_pts_t1, pts2, match2, num_pts2);
        }
        
        count1 = count2 = 0;
        
        for (j = 0; j < num_pts_t1; j++)
        {
          if (!Q1.sgn_at(t_h[j]))
          {
            b[0]     = s_h[i];
            b[1]     = t_h[j];
            sgn_dtds = - sgn(P1.derivative(0).evaluate(b)) *
                         sgn(P1.derivative(1).evaluate(b));
            
            if (sgn_dtds > 0)
            {
              boxhits[i][j][boxhits1[i][j] + boxhits2[i][j]] = 1;
              boxhits1[i][j]++;
            }
            else if (!sgn_dtds)
              boxhits1[i][j] += 5;
          }
          
          if (!Q2.sgn_at(t_h[j]))
          {
            b[0]     = s_h[i];
            b[1]     = t_h[j];
            sgn_dtds = - sgn(P2.derivative(0).evaluate(b)) *
                         sgn(P2.derivative(1).evaluate(b));
            
            if (sgn_dtds > 0)
            {
              boxhits[i][j][boxhits1[i][j] + boxhits2[i][j]] = 2;
              boxhits2[i][j]++;
            }
            else if (!sgn_dtds)
              boxhits2[i][j] += 5;
          }
          
          pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
          num_pts_hit = 0;
          
          while (count1 < num_pts1 && match1[count1] < j)
            count1++;
          
          while (count1 < num_pts1 && match1[count1] == j)
          {
            boxhits1[i][j]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts1[count1++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          while (count2 < num_pts2 && match2[count2] < j)
            count2++;
          
          while (count2 < num_pts2 && match2[count2] == j)
          {
            boxhits2[i][j]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts2[count2++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          if (!sort(pts_hit, num_pts_hit))
          {
            for (k = 0; k < num_pts_hit - 1; k++)
              for (l = k + 1; l < num_pts_hit; l++)
                pts_hit[k]->separate(*pts_hit[l]);
            
            assert(sort(pts_hit, num_pts_hit));
          }
          
          for (k = num_pts_hit - 1; k >= 0; k--)  //  in reverse order
          {
            l = boxhits1[i][j] + boxhits2[i][j] - k - 1;
            
            if (poly_hit && *pts_hit[k]->poly == *poly_hit)
              boxhits[i][j][l] = 1;
            else  //  if (!poly_hit || *pts_hit[k]->poly != *poly_hit)
              boxhits[i][j][l] = 2;
            
            if (!--pts_hit[k]->ref_count)
              delete pts_hit[k];
          }
          
          delete [] pts_hit;
        }
        
        for (j = 0; j < num_pts1; j++)
          if (!--pts1[j]->ref_count)
            delete pts1[j];
        
        if (num_pts1 > 0)
        {
          delete [] pts1;
          delete [] match1;
        }
        
        for (j = 0; j < num_pts2; j++)
          if (!--pts2[j]->ref_count)
            delete pts2[j];
        
        if (num_pts2 > 0)
        {
          delete [] pts2;
          delete [] match2;
        }
      }
      
      //  low_t: the bottom side of a box (counted in reverse order)
      
      for (i = 0; i < num_pts_t1; i++)
      {
        Q1       = P1.subst_val(1, t_l[i]);
        num_pts1 = get_pts(s_l[0], s_h[num_pts_s1 - 1], Q1, pts1, init_tol, 0);
        
        if (num_pts1 > 0)
        {
          poly_hit = pts1[0]->poly;
          sort(pts1, num_pts1);
          match_intervals(s_l, s_h, num_pts_s1, pts1, match1, num_pts1);
        }
        else  //  if (num_pts1 > 0)
          poly_hit = 0;
        
        Q2       = P2.subst_val(1, t_l[i]);
        num_pts2 = get_pts(s_l[0], s_h[num_pts_s1 - 1], Q2, pts2, init_tol, 0);
        
        if (num_pts2 > 0)
        {
          sort(pts2, num_pts2);
          match_intervals(s_l, s_h, num_pts_s1, pts2, match2, num_pts2);
        }
        
        count1 = count2 = 0;
        
        for (j = 0; j < num_pts_s1; j++)
        {
          if (!Q1.sgn_at(s_h[j]))
          {
            b[0]     = s_h[j];
            b[1]     = t_l[i];
            sgn_dtds = - sgn(P1.derivative(0).evaluate(b)) *
                         sgn(P1.derivative(1).evaluate(b));
            
            if (sgn_dtds < 0)
            {
              boxhits[j][i][boxhits1[j][i] + boxhits2[j][i]] = 1;
              boxhits1[j][i]++;
            }
            else if (!sgn_dtds)
              boxhits1[j][i] += 5;
          }
          
          if (!Q2.sgn_at(s_h[j]))
          {
            b[0]     = s_h[j];
            b[1]     = t_l[i];
            sgn_dtds = - sgn(P2.derivative(0).evaluate(b)) *
                         sgn(P2.derivative(1).evaluate(b));
            
            if (sgn_dtds < 0)
            {
              boxhits[j][i][boxhits1[j][i] + boxhits2[j][i]] = 2;
              boxhits2[j][i]++;
            }
            else if (!sgn_dtds)
              boxhits2[j][i] += 5;
          }
          
          pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
          num_pts_hit = 0;
          
          while (count1 < num_pts1 && match1[count1] < j)
            count1++;
          
          while (count1 < num_pts1 && match1[count1] == j)
          {
            boxhits1[j][i]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts1[count1++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          while (count2 < num_pts2 && match2[count2] < j)
            count2++;
          
          while (count2 < num_pts2 && match2[count2] == j)
          {
            boxhits2[j][i]++;
//            assert(num_pts_hit < MAX_NUM_PTS_HIT);
            pts_hit[num_pts_hit] = pts2[count2++];
            pts_hit[num_pts_hit++]->ref_count++;
          }
          
          if (!sort(pts_hit, num_pts_hit))
          {
            for (k = 0; k < num_pts_hit - 1; k++)
              for (l = k + 1; l < num_pts_hit; l++)
                pts_hit[k]->separate(*pts_hit[l]);
            
            assert(sort(pts_hit, num_pts_hit));
          }
          
          for (k = num_pts_hit - 1; k >= 0; k--)  //  in reverse order
          {
            l = boxhits1[j][i] + boxhits2[j][i] - k - 1;
            
            if (poly_hit && *pts_hit[k]->poly == *poly_hit)
              boxhits[j][i][l] = 1;
            else  //  if (!poly_hit || *pts_hit[k]->poly != *poly_hit)
              boxhits[j][i][l] = 2;
            
            if (!--pts_hit[k]->ref_count)
              delete pts_hit[k];
          }
          
          delete [] pts_hit;
        }
        
        for (j = 0; j < num_pts1; j++)
          if (!--pts1[j]->ref_count)
            delete pts1[j];
        
        if (num_pts1 > 0)
        {
          delete [] pts1;
          delete [] match1;
        }
        
        for (j = 0; j < num_pts2; j++)
          if (!--pts2[j]->ref_count)
            delete pts2[j];
        
        if (num_pts2 > 0)
        {
          delete [] pts2;
          delete [] match2;
        }
      }
      
      //  See which of those boxes have roots in them.
      
      for (i = 0; i < num_pts_s1; i++)
        for (j = 0; j < num_pts_t1; j++)
          if (boxhits1[i][j] >= 2 && boxhits2[i][j] >= 2)
            if (boxhits1[i][j] + boxhits2[i][j] > 4)
            {
              if (refine_interior(pts_s1[i], pts_t1[j],
                                  P1, P2,
                                  pts_proto[num_pts],
                                  tol))
                num_pts++;
            }
            else if (boxhits[i][j][0] == 1 && boxhits[i][j][1] == 2 &&
                     boxhits[i][j][2] == 1 && boxhits[i][j][3] == 2
                     ||
                     boxhits[i][j][0] == 2 && boxhits[i][j][1] == 1 &&
                     boxhits[i][j][2] == 2 && boxhits[i][j][3] == 1)
            {
              pts_proto[num_pts] = new K_POINT2D(*pts_s1[i]->PtR,
                                                 *pts_t1[j]->PtR,
                                                 P1,
                                                 P2);
              pts_proto[num_pts++]->ref_count++;
            }
      
      for (i = 0; i < num_pts_s1; i++)
      {
        for (j = 0; j < num_pts_t1; j++)
          delete [] boxhits[i][j];
        
        delete [] boxhits[i];
        delete [] boxhits1[i];
        delete [] boxhits2[i];
      }
      
      delete [] boxhits;
      delete [] boxhits1;
      delete [] boxhits2;
      
      delete [] s_l;
      delete [] s_h;
      delete [] t_l;
      delete [] t_h;
    }
    
    for (i = 0; i < num_pts_s1; i++)
      if (!--pts_s1[i]->ref_count)
        delete pts_s1[i];
    
    delete [] pts_s1;  //  num_pts_s > 0 => pts_s1 != 0
    
    for (i = 0; i < num_pts_s2; i++)
      if (!--pts_s2[i]->ref_count)
        delete pts_s2[i];
    
    delete [] pts_s2;  //  num_pts_s > 0 => pts_s2 != 0
    
    for (i = 0; i < num_pts_t1; i++)
      if (!--pts_t1[i]->ref_count)
        delete pts_t1[i];
    
    delete [] pts_t1;  //  num_pts_t > 0 => pts_t1 != 0
    
    for (i = 0; i < num_pts_t2; i++)
      if (!--pts_t2[i]->ref_count)
        delete pts_t2[i];
    
    delete [] pts_t2;  //  num_pts_t > 0 => pts_t2 != 0
    
    if (num_pts > 0)
    {
      pts = new K_POINT2D* [num_pts];
      
      for (i = 0; i < num_pts; i++)
      {
        pts[i] = pts_proto[i];
        pts[i]->ref_count++;
      }
    }
    else  //  if (!num_pts)
      pts = 0;
    
    for (i = 0; i < num_pts; i++)
      if (!--pts_proto[i]->ref_count)
        delete pts_proto[i];
    
    delete [] pts_proto;  //  num_pts_all >= 1 => pts_proto != 0
  }
  
  return num_pts;
}

int refine_interior(K_POINT1D* const x_s, K_POINT1D* const x_t,
                    const K_RATPOLY& P1, const K_RATPOLY& P2,
                    K_POINT2D*& y,
                    const bigrational& tol)
{
  assert(x_s->type);
  assert(x_t->type);
  assert(P1.num_vars == 2);
  assert(P2.num_vars == 2);
  
  long i, j;
  int  r;
  
  r = 0;
  
  if (x_s->type == 1)
    x_s->shrink(shrink_step);
  
  if (x_t->type == 1)
    x_t->shrink(shrink_step);
  
  if (x_s->type == 2 && x_t->type == 2)
  {
    y = new K_POINT2D(x_s->PtB, x_t->PtB, P1, P2);
    y->ref_count++;
    r = 1;
  }
  else if (x_s->type == 2 && x_t->type == 1)
  {
    bigrational   b_s;
    K_RATPOLY     G21;
    unsigned long num_pts_21;
    K_POINT1D**   pts_21;
    
    b_s = x_s->PtB;
    G21 = gcd(P1.subst_val(0, b_s), P2.subst_val(0, b_s));
    
    if (G21.deg[0] > 0)
    {
      num_pts_21 =
        get_pts(x_t->get_low(), x_t->get_high(), G21, pts_21, tol, 0);
      
      assert(num_pts_21 <= 1);
      
      if (num_pts_21 == 1)
      {
        y = new K_POINT2D(b_s, *pts_21[0], P1, P2);
        y->ref_count++;
        r = 1;
        
        if (!--pts_21[0]->ref_count)
          delete pts_21[0];
        
        delete [] pts_21;
      }
    }
  }
  else if (x_s->type == 1 && x_t->type == 2)
  {
    bigrational   b_t;
    K_RATPOLY     G12;
    unsigned long num_pts_12;
    K_POINT1D**   pts_12;
    
    b_t = x_t->PtB;
    G12 = gcd(P1.subst_val(1, b_t), P2.subst_val(1, b_t));
    
    if (G12.deg[0] > 0)
    {
      num_pts_12 =
        get_pts(x_s->get_low(), x_s->get_high(), G12, pts_12, tol, 0);
      
      assert(num_pts_12 <= 1);
      
      if (num_pts_12 == 1)
      {
        y = new K_POINT2D(*pts_12[0], b_t, P1, P2);
        y->ref_count++;
        r = 1;
        
        if (!--pts_12[0]->ref_count)
          delete pts_12[0];
        
        delete [] pts_12;
      }
    }
  }
  else  //  if (x_s->type == 1 && x_t->type == 1)
  {
    unsigned long* boxhits;
    unsigned long  boxhits1;
    unsigned long  boxhits2;
    
    boxhits  = new unsigned long [MAX_NUM_PTS_HIT];
    boxhits1 = boxhits2 = 0;
    
    for (i = 0; i < MAX_NUM_PTS_HIT; i++)
      boxhits[i] = 0;
    
    K_RATPOLY     Q1, Q2;
    unsigned long num_pts1, num_pts2;
    K_POINT1D**   pts1;
    K_POINT1D**   pts2;
    K_RATPOLY*    poly_hit;
    bigrational   v[2];
    int           sgn_dtds;
    unsigned long count1, count2;
    unsigned long num_pts_hit;
    K_POINT1D**   pts_hit;
    
    //  low_s
    
    Q1       = P1.subst_val(0, x_s->get_low());
    num_pts1 = get_pts(x_t->get_low(), x_t->get_high(), Q1, pts1, init_tol, 0);
    
    if (num_pts1 > 0)
    {
      poly_hit = pts1[0]->poly;
      sort(pts1, num_pts1);
    }
    else  //  if (num_pts1 == 0)
      poly_hit = 0;
    
    Q2       = P2.subst_val(0, x_s->get_low());
    num_pts2 = get_pts(x_t->get_low(), x_t->get_high(), Q2, pts2, init_tol, 0);
    
    if (num_pts2 > 0)
      sort(pts2, num_pts2);
    
    if (!Q1.sgn_at(x_t->get_low()))
    {
      v[0]     = x_s->get_low();
      v[1]     = x_t->get_low();
      sgn_dtds = - sgn(P1.derivative(0).evaluate(v)) *
                   sgn(P1.derivative(1).evaluate(v));
      
      if (sgn_dtds > 0)
      {
        boxhits[boxhits1 + boxhits2] = 1;
        boxhits1++;
      }
      else if (!sgn_dtds)
        boxhits1 += 5;
    }
    
    if (!Q2.sgn_at(x_t->get_low()))
    {
      v[0]     = x_s->get_low();
      v[1]     = x_t->get_low();
      sgn_dtds = - sgn(P2.derivative(0).evaluate(v)) *
                   sgn(P2.derivative(1).evaluate(v));
      
      if (sgn_dtds > 0)
      {
        boxhits[boxhits1 + boxhits2] = 2;
        boxhits2++;
      }
      else if (!sgn_dtds)
        boxhits2 += 5;
    }
    
    pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
    num_pts_hit = 0;
    
    for (count1 = 0; count1 < num_pts1; count1++)
    {
      boxhits1++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts1[count1];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    for (count2 = 0; count2 < num_pts2; count2++)
    {
      boxhits2++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts2[count2];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    if (!sort(pts_hit, num_pts_hit))
    {
      for (i = 0; i < num_pts_hit - 1; i++)
        for (j = i + 1; j < num_pts_hit; j++)
          pts_hit[i]->separate(*pts_hit[j]);
      
      assert(sort(pts_hit, num_pts_hit));
    }
    
    for (i = 0; i < num_pts_hit; i++)
    {
      j = boxhits1 + boxhits2 - num_pts_hit + i;
      
      if (poly_hit && *pts_hit[i]->poly == *poly_hit)
        boxhits[j] = 1;
      else  //  if (!poly_hit || *pts_hit[i]->poly != *poly_hit)
        boxhits[j] = 2;
      
      if (!--pts_hit[i]->ref_count)
        delete pts_hit[i];
    }
    
    delete [] pts_hit;
    
    for (i = 0; i < num_pts1; i++)
      if (!--pts1[i]->ref_count)
        delete pts1[i];
    
    if (num_pts1 > 0)
      delete [] pts1;
    
    for (i = 0; i < num_pts2; i++)
      if (!--pts2[i]->ref_count)
        delete pts2[i];
    
    if (num_pts2 > 0)
      delete [] pts2;
    
    //  high_t
    
    Q1       = P1.subst_val(1, x_t->get_high());
    num_pts1 = get_pts(x_s->get_low(), x_s->get_high(), Q1, pts1, init_tol, 0);
    
    if (num_pts1 > 0)
    {
      poly_hit = pts1[0]->poly;
      sort(pts1, num_pts1);
    }
    else  //  if (num_pts1 == 0)
      poly_hit = 0;
    
    Q2       = P2.subst_val(1, x_t->get_high());
    num_pts2 = get_pts(x_s->get_low(), x_s->get_high(), Q2, pts2, init_tol, 0);
    
    if (num_pts2 > 0)
      sort(pts2, num_pts2);
    
    if (!Q1.sgn_at(x_s->get_low()))
    {
      v[0]     = x_s->get_low();
      v[1]     = x_t->get_high();
      sgn_dtds = - sgn(P1.derivative(0).evaluate(v)) *
                   sgn(P1.derivative(1).evaluate(v));
      
      if (sgn_dtds < 0)
      {
        boxhits[boxhits1 + boxhits2] = 1;
        boxhits1++;
      }
      else if (!sgn_dtds)
        boxhits1 += 5;
    }
    
    if (!Q2.sgn_at(x_s->get_low()))
    {
      v[0]     = x_s->get_low();
      v[1]     = x_t->get_high();
      sgn_dtds = - sgn(P2.derivative(0).evaluate(v)) *
                   sgn(P2.derivative(1).evaluate(v));
      
      if (sgn_dtds < 0)
      {
        boxhits[boxhits1 + boxhits2] = 2;
        boxhits2++;
      }
      else if (!sgn_dtds)
        boxhits2 += 5;
    }
    
    pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
    num_pts_hit = 0;
    
    for (count1 = 0; count1 < num_pts1; count1++)
    {
      boxhits1++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts1[count1];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    for (count2 = 0; count2 < num_pts2; count2++)
    {
      boxhits2++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts2[count2];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    if (!sort(pts_hit, num_pts_hit))
    {
      for (i = 0; i < num_pts_hit - 1; i++)
        for (j = i + 1; j < num_pts_hit; j++)
          pts_hit[i]->separate(*pts_hit[j]);
      
      assert(sort(pts_hit, num_pts_hit));
    }
    
    for (i = 0; i < num_pts_hit; i++)
    {
      j = boxhits1 + boxhits2 - num_pts_hit + i;
      
      if (poly_hit && *pts_hit[i]->poly == *poly_hit)
        boxhits[j] = 1;
      else  //  if (!poly_hit || *pts_hit[i]->poly != *poly_hit)
        boxhits[j] = 2;
      
      if (!--pts_hit[i]->ref_count)
        delete pts_hit[i];
    }
    
    delete [] pts_hit;
    
    for (i = 0; i < num_pts1; i++)
      if (!--pts1[i]->ref_count)
        delete pts1[i];
    
    if (num_pts1 > 0)
      delete [] pts1;
    
    for (i = 0; i < num_pts2; i++)
      if (!--pts2[i]->ref_count)
        delete pts2[i];
    
    if (num_pts2 > 0)
      delete [] pts2;
    
    //  high_s (in reverse order)
    
    Q1       = P1.subst_val(0, x_s->get_high());
    num_pts1 = get_pts(x_t->get_low(), x_t->get_high(), Q1, pts1, init_tol, 0);
    
    if (num_pts1 > 0)
    {
      poly_hit = pts1[0]->poly;
      sort(pts1, num_pts1);
    }
    else  //  if (num_pts1 == 0)
      poly_hit = 0;
    
    Q2       = P2.subst_val(0, x_s->get_high());
    num_pts2 = get_pts(x_t->get_low(), x_t->get_high(), Q2, pts2, init_tol, 0);
    
    if (num_pts2 > 0)
      sort(pts2, num_pts2);
    
    if (!Q1.sgn_at(x_t->get_high()))
    {
      v[0]     = x_s->get_high();
      v[1]     = x_t->get_high();
      sgn_dtds = - sgn(P1.derivative(0).evaluate(v)) *
                   sgn(P1.derivative(1).evaluate(v));
      
      if (sgn_dtds > 0)
      {
        boxhits[boxhits1 + boxhits2] = 1;
        boxhits1++;
      }
      else if (!sgn_dtds)
        boxhits1 += 5;
    }
    
    if (!Q2.sgn_at(x_t->get_high()))
    {
      v[0]     = x_s->get_high();
      v[1]     = x_t->get_high();
      sgn_dtds = - sgn(P2.derivative(0).evaluate(v)) *
                   sgn(P2.derivative(1).evaluate(v));
      
      if (sgn_dtds > 0)
      {
        boxhits[boxhits1 + boxhits2] = 2;
        boxhits2++;
      }
      else if (!sgn_dtds)
        boxhits2 += 5;
    }
    
    pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
    num_pts_hit = 0;
    
    for (count1 = 0; count1 < num_pts1; count1++)
    {
      boxhits1++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts1[count1];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    for (count2 = 0; count2 < num_pts2; count2++)
    {
      boxhits2++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts2[count2];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    if (!sort(pts_hit, num_pts_hit))
    {
      for (i = 0; i < num_pts_hit - 1; i++)
        for (j = i + 1; j < num_pts_hit; j++)
          pts_hit[i]->separate(*pts_hit[j]);
      
      assert(sort(pts_hit, num_pts_hit));
    }
    
    for (i = num_pts_hit - 1; i >= 0; i--)  //  in reverse order
    {
      j = boxhits1 + boxhits2 - i - 1;
      
      if (poly_hit && *pts_hit[i]->poly == *poly_hit)
        boxhits[j] = 1;
      else  //  if (!poly_hit || *pts_hit[i]->poly != *poly_hit)
        boxhits[j] = 2;
      
      if (!--pts_hit[i]->ref_count)
        delete pts_hit[i];
    }
    
    delete [] pts_hit;
    
    for (i = 0; i < num_pts1; i++)
      if (!--pts1[i]->ref_count)
        delete pts1[i];
    
    if (num_pts1 > 0)
      delete [] pts1;
    
    for (i = 0; i < num_pts2; i++)
      if (!--pts2[i]->ref_count)
        delete pts2[i];
    
    if (num_pts2 > 0)
      delete [] pts2;
    
    //  low_t (in reverse order)
    
    Q1       = P1.subst_val(1, x_t->get_low());
    num_pts1 = get_pts(x_s->get_low(), x_s->get_high(), Q1, pts1, init_tol, 0);
    
    if (num_pts1 > 0)
    {
      poly_hit = pts1[0]->poly;
      sort(pts1, num_pts1);
    }
    else  //  if (num_pts1 == 0)
      poly_hit = 0;
    
    Q2       = P2.subst_val(1, x_t->get_low());
    num_pts2 = get_pts(x_s->get_low(), x_s->get_high(), Q2, pts2, init_tol, 0);
    
    if (num_pts2 > 0)
      sort(pts2, num_pts2);
    
    if (!Q1.sgn_at(x_s->get_high()))
    {
      v[0]     = x_s->get_high();
      v[1]     = x_t->get_low();
      sgn_dtds = - sgn(P1.derivative(0).evaluate(v)) *
                   sgn(P1.derivative(1).evaluate(v));
      
      if (sgn_dtds < 0)
      {
        boxhits[boxhits1 + boxhits2] = 1;
        boxhits1++;
      }
      else if (!sgn_dtds)
        boxhits1 += 5;
    }
    
    if (!Q2.sgn_at(x_s->get_high()))
    {
      v[0]     = x_s->get_high();
      v[1]     = x_t->get_low();
      sgn_dtds = - sgn(P2.derivative(0).evaluate(v)) *
                   sgn(P2.derivative(1).evaluate(v));
      
      if (sgn_dtds < 0)
      {
        boxhits[boxhits1 + boxhits2] = 2;
        boxhits2++;
      }
      else if (!sgn_dtds)
        boxhits2 += 5;
    }
    
    pts_hit     = new K_POINT1D* [MAX_NUM_PTS_HIT];
    num_pts_hit = 0;
    
    for (count1 = 0; count1 < num_pts1; count1++)
    {
      boxhits1++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts1[count1];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    for (count2 = 0; count2 < num_pts2; count2++)
    {
      boxhits2++;
//      assert(num_pts_hit < MAX_NUM_PTS_HIT);
      pts_hit[num_pts_hit] = pts2[count2];
      pts_hit[num_pts_hit++]->ref_count++;
    }
    
    if (!sort(pts_hit, num_pts_hit))
    {
      for (i = 0; i < num_pts_hit - 1; i++)
        for (j = i + 1; j < num_pts_hit; j++)
          pts_hit[i]->separate(*pts_hit[j]);
      
      assert(sort(pts_hit, num_pts_hit));
    }
    
    for (i = num_pts_hit - 1; i >= 0; i--)
    {
      j = boxhits1 + boxhits2 - i - 1;
      
      if (poly_hit && *pts_hit[i]->poly == *poly_hit)
        boxhits[j] = 1;
      else  //  if (!poly_hit || *pts_hit[i]->poly != *poly_hit)
        boxhits[j] = 2;
      
      if (!--pts_hit[i]->ref_count)
        delete pts_hit[i];
    }
    
    delete [] pts_hit;
    
    for (i = 0; i < num_pts1; i++)
      if (!--pts1[i]->ref_count)
        delete pts1[i];
    
    if (num_pts1 > 0)
      delete [] pts1;
    
    for (i = 0; i < num_pts2; i++)
      if (!--pts2[i]->ref_count)
        delete pts2[i];
    
    if (num_pts2 > 0)
      delete [] pts2;
    
    //  See which of those boxes have roots in them.
    
    if (boxhits1 >= 2 &&  boxhits2 >= 2)
      if (boxhits1 + boxhits2 > 4)
        r = refine_interior(x_s, x_t, P1, P2, y, tol);
      else if (boxhits[0] == 1 && boxhits[1] == 2 &&
               boxhits[2] == 1 && boxhits[3] == 2
               ||
               boxhits[0] == 2 && boxhits[1] == 1 &&
               boxhits[2] == 2 && boxhits[3] == 1)
      {
        y = new K_POINT2D(*x_s->PtR, *x_t->PtR, P1, P2);
        y->ref_count++;
        r = 1;
      }
    
    delete [] boxhits;
  }
  
  return r;
}

//  unsigned long get_pts_proto(const bigrational& l_s, const bigrational& h_s,
//                              const K_RATPOLY& s_only,
//                              const bigrational& b_t,
//                              const K_RATPOLY& P1, const K_RATPOLY& P2,
//                              K_POINT2D**& pts,
//                              const bigrational& tol, const int count_edges)
//    computes the intersections "pts" of 2 bivariate polynomials P1 and P2
//                on/in the horizontal line segment [(l_s, b_t), (h_s, b_t)].
//    by computing the roots of a univariate polynomial s_only
//                 on/in the interval [l_s, h_s].
//    returns the number of intersections.
//    if tol >  0 then
//      a box for each intersection is no larger than tol in the s-direction.
//    if tol = 0 then there is no limit on the size of boxes for intersections.
//    if count_edges = 1 then
//      the intersections at the endpoints of the line segment are counted.
//    if count_edges = 0 then
//      the intersections at the endpoints of the line segment are ignored.

unsigned long get_pts_proto(const bigrational& l_s, const bigrational& h_s,
                            const K_RATPOLY& s_only,
                            const bigrational& b_t,
                            const K_RATPOLY& P1, const K_RATPOLY& P2,
                            K_POINT2D**& pts,
                            const bigrational& tol, const int count_edges)
{
  assert(s_only.num_vars == 1);
  
  unsigned long num_pts;
  
  if (!s_only.deg[0])
  //  0.1.  when s_only is a constant polynomial
  {
    pts     = 0;
    num_pts = 0;
  }
//  else if (!P1.deg[0] && !P1.deg[1] || !P2.deg[0] && !P2.deg[1])
//  //  0.2.  when P1 or P2 is a constant polynomial
//  {
//    pts     = 0;
//    num_pts = 0;
//  }
//  else if (P1.eq_upto_const(P2))
//  //  0.3.  when P1 and P2 are the same (up to constant multiple)
//  {
//    pts     = 0;
//    num_pts = 0;
//  }
  else if (l_s > h_s)
  //  0.4.  when the horizontal line segment [(l_s, b_t), (h_s, b_t)] is
  //          not well-defined
  {
    pts     = 0;
    num_pts = 0;
  }
  else if (l_s == h_s)
  //  1.  when the horizontal line segment [(l_s, b_t), (h_s, b_t)] is
  //        shrunken to a signle point (l_s, b_t).
    if (count_edges && !s_only.sgn_at(l_s))
    {
      pts    = new K_POINT2D* [num_pts = 1];
      pts[0] = new K_POINT2D(l_s, b_t, P1, P2);
      pts[0]->ref_count++;
    }
    else  //  if (!count_edges || s_only.sgn_at(l_s))
    {
      pts     = 0;
      num_pts = 0;
    }
  else if (s_only.deg[0] == 1)
  //  2.  when s_only is linear
  {
    bigrational b_s;
    
    b_s = - s_only.coeffs[1] / s_only.coeffs[0];
    //  s_only.deg[0] == 1 => s_only.coeffs[0] != 0
    
    if (b_s > l_s && b_s < h_s || count_edges && b_s >= l_s && b_s <= h_s)
    {
      pts    = new K_POINT2D* [num_pts = 1];
      pts[0] = new K_POINT2D(b_s, b_t, P1, P2);
      pts[0]->ref_count++;
    }
    else  //  if b_s is not in/on the interval [l_s, h_s]
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else
  //  3.  Otherwise, perform 1-D root finding.
  {
    unsigned long i;
    unsigned long num_all, num_each;
    ROOT1*        each;
    K_POINT2D**   pts_proto;
    bigrational   e_l;
    
    //  3.1  Find all roots of s_only on/in the half open interval [l_s, h_s)
    //       and
    //       isolate them.
    
    ROOT1 all(s_only, l_s, h_s);
    
    num_all  = all.num_roots;
    num_each = all.isolate_roots(each, tol);
    
    assert(num_each == num_all);
    
    pts_proto = new K_POINT2D* [num_all + 1];
    num_pts   = 0;
    
    for (i = 0; i < num_all; i++)
    {
      e_l = each[i].low;
      
      if (!s_only.sgn_at(e_l))
      //  3.1.1  A root e_l = each[i].low of s_only will belong to pts_proto
      //         iff it is not l_s unless count_edges == 1.
      {
        if (count_edges || e_l != l_s)
        {
          pts_proto[num_pts] = new K_POINT2D(e_l, b_t, P1, P2);
          pts_proto[num_pts++]->ref_count++;
        }
      }
      else  //  if (s_only.sgn_at(e_l))
      {
        pts_proto[num_pts] = new K_POINT2D(each[i], b_t, P1, P2);
        pts_proto[num_pts++]->ref_count++;
      }
    }
    
    if (num_each > 0)
      delete [] each;
    
    //  3.2.  See if h_s is a root of s_only.
    
    if (count_edges && !s_only.sgn_at(h_s))
    {
      pts_proto[num_pts] = new K_POINT2D(h_s, b_t, P1, P2);
      pts_proto[num_pts++]->ref_count++;
    }
    
    //  3.3.  Gather all together.
    
    if (num_pts > 0)
    {
      pts = new K_POINT2D* [num_pts];
      
      for (i = 0; i < num_pts; i++)
      {
        pts[i] = pts_proto[i];
        pts[i]->ref_count++;
      }
    }
    else  //  if (!num_pts)
      pts = 0;
    
    for (i = 0; i < num_pts; i++)
      if (!--pts_proto[i]->ref_count)
        delete pts_proto[i];
    
    delete [] pts_proto;  //  pts_proto != 0
  }
  
  return num_pts;
}

//  unsigned long get_pts_proto(const bigrational& b_s,
//                              const bigrational& l_t, const bigrational& h_t,
//                              const K_RATPOLY& t_only,
//                              const K_RATPOLY& P1, const K_RATPOLY& P2,
//                              K_POINT2D**& pts,
//                              const bigrational& tol, const int count_edges)
//    computes the intersections "pts" of 2 bivariate polynomials P1 and P2
//                on/in the vertical line segment [(b_s, l_t), (b_s, h_t)].
//    by computing the roots of a univariate polynomial t_only
//                 on/in the interval [l_t, h_t].
//    returns the number of intersections.
//    if tol >  0 then
//      a box for each intersection is no larger than tol in the t-direction.
//    if tol = 0 then there is no limit on the size of boxes for intersections.
//    if count_edges = 1 then
//      the intersections at the endpoints of the line segment are counted.
//    if count_edges = 0 then
//      the intersections at the endpoints of the line segment are ignored.

unsigned long get_pts_proto(const bigrational& b_s,
                            const bigrational& l_t, const bigrational& h_t,
                            const K_RATPOLY& t_only,
                            const K_RATPOLY& P1, const K_RATPOLY& P2,
                            K_POINT2D**& pts,
                            const bigrational& tol, const int count_edges)
{
  assert(t_only.num_vars == 1);
  
  unsigned long num_pts;
  
  if (!t_only.deg[0])
  //  0.1.  when t_only is a constant polynomial
  {
    pts     = 0;
    num_pts = 0;
  }
//  else if (!P1.deg[0] && !P1.deg[1] || !P2.deg[0] && !P2.deg[1])
//  //  0.2.  when P1 or P2 is a constant polynomial
//  {
//    pts     = 0;
//    num_pts = 0;
//  }
//  else if (P1.eq_upto_const(P2))
//  //  0.3.  when P1 and P2 are the same (up to constant multiple)
//  {
//    pts     = 0;
//    num_pts = 0;
//  }
  else if (l_t > h_t)
  //  0.4.  when the vertical line segment [(b_s, l_t), (b_s, h_t)] is
  //          not well-defined
  {
    pts     = 0;
    num_pts = 0;
  }
  else if (l_t == h_t)
  //  1.  when the vertical line segment [(b_s, l_t), (b_s, h_t)] is
  //        shrunken to a signle point (l_t, b_t).
    if (count_edges && !t_only.sgn_at(l_t))
    {
      pts    = new K_POINT2D* [num_pts = 1];
      pts[0] = new K_POINT2D(b_s, l_t, P1, P2);
      pts[0]->ref_count++;
    }
    else  //  if (!count_edges || t_only.sgn_at(l_t))
    {
      pts     = 0;
      num_pts = 0;
    }
  else if (t_only.deg[0] == 1)
  //  2.  when t_only is linear
  {
    bigrational b_t;
    
    b_t = - t_only.coeffs[1] / t_only.coeffs[0];
    //  t_only.deg[0] == 1 => t_only.coeffs[0] != 0
    
    if (b_t > l_t && b_t < h_t || count_edges && b_t >= l_t && b_t <= h_t)
    {
      pts    = new K_POINT2D* [num_pts = 1];
      pts[0] = new K_POINT2D(b_s, b_t, P1, P2);
      pts[0]->ref_count++;
    }
    else  //  if b_t is not in/on the interval [l_t, h_t]
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else
  //  3.  Otherwise, perform 1-D root finding.
  {
    unsigned long i;
    unsigned long num_all, num_each;
    ROOT1*        each;
    K_POINT2D**   pts_proto;
    bigrational   e_l;
    
    //  3.1  Find all roots of t_only on/in the half open intrval [l_t, h_t)
    //       and
    //       isolate them.
    
    ROOT1 all(t_only, l_t, h_t);
    
    num_all  = all.num_roots;
    num_each = all.isolate_roots(each, tol);
    assert(num_each == num_all);
    
    pts_proto = new K_POINT2D* [num_all + 1];
    num_pts   = 0;
    
    for (i = 0; i < num_all; i++)
    {
      e_l = each[i].low;
      
      if (!t_only.sgn_at(e_l))
      //  3.1.1  A root e_l = each[i].low of t_only will belong to pts_proto
      //         iff it is not l_t unless count_edges == 1.
      {
        if (count_edges || e_l != l_t)
        {
          pts_proto[num_pts] = new K_POINT2D(b_s, e_l, P1, P2);
          pts_proto[num_pts++]->ref_count++;
        }
      }
      else  //  if (t_only.sgn_at(e_l))
      {
        pts_proto[num_pts] = new K_POINT2D(b_s, each[i], P1, P2);
        pts_proto[num_pts++]->ref_count++;
      }
    }
    
    if (num_each > 0)
      delete [] each;
    
    //  3.2.  See if h_t is a root of t_only.
    
    if (count_edges && !t_only.sgn_at(h_t))
    {
      pts_proto[num_pts] = new K_POINT2D(b_s, h_t, P1, P2);
      pts_proto[num_pts++]->ref_count++;
    }
    
    //  3.3.  Gather all together.
    
    if (num_pts > 0)
    {
      pts = new K_POINT2D* [num_pts];
      
      for (i = 0; i < num_pts; i++)
      {
        pts[i] = pts_proto[i];
        pts[i]->ref_count++;
      }
    }
    else  //  if (!num_pts)
      pts = 0;
    
    for (i = 0; i < num_pts; i++)
      if (!--pts_proto[i]->ref_count)
        delete pts_proto[i];
    
    delete [] pts_proto;  //  pts_proto != 0
  }
  
  return num_pts;
}

//  K_BOXCO2 K_POINT2D :: bbox() const
//    return a bounding box for *this.

K_BOXCO2 K_POINT2D :: bbox() const
{
  assert(type);
  
  K_BOXCO2 b2;
  
  if (type == 1)
    b2 = K_BOXCO2(PtRs->low, PtRs->high, PtRt->low, PtRt->high, 1, 1, 1, 1);
  else if (type == 2)
    b2 = K_BOXCO2(PtRs->low, PtRs->high, PtBt, PtBt, 1, 1, 0, 0);
  else if (type == 3)
    b2 = K_BOXCO2(PtBs, PtBs, PtRt->low, PtRt->high, 0, 0, 1, 1);
  else  //  if (type == 4)
    b2 = K_BOXCO2(PtBs, PtBs, PtBt, PtBt, 0, 0, 0, 0);
  
  return b2;
}

//  int K_POINT2D :: overlap(const K_POINT2D& x) const
//    return 1 if *this and x overlap, and
//           0 otherwise.

int K_POINT2D :: overlap(const K_POINT2D& x) const
{
  assert(type);
  assert(x.type);
  
  int o;
  
  if (type == 1 && x.type == 1)
    if ((PtRs->low >= x.PtRs->low && PtRs->low < x.PtRs->high
         ||
         PtRs->high > x.PtRs->low && PtRs->high <= x.PtRs->high
         ||
         PtRs->low <= x.PtRs->low && PtRs->high > x.PtRs->low
         ||
         PtRs->low < x.PtRs->high && PtRs->high >= x.PtRs->high)
         &&
        (PtRt->low >= x.PtRt->low && PtRt->low < x.PtRt->high
         ||
         PtRt->high > x.PtRt->low && PtRt->high <= x.PtRt->high
         ||
         PtRt->low <= x.PtRt->low && PtRt->high > x.PtRt->low
         ||
         PtRt->low < x.PtRt->high && PtRt->high >= x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 1 && x.type == 2)
    if ((PtRs->low >= x.PtRs->low && PtRs->low < x.PtRs->high
         ||
         PtRs->high > x.PtRs->low && PtRs->high <= x.PtRs->high
         ||
         PtRs->low <= x.PtRs->low && PtRs->high > x.PtRs->low
         ||
         PtRs->low < x.PtRs->high && PtRs->high >= x.PtRs->high)
         &&
        (PtRt->low < x.PtBt && PtRt->high > x.PtBt))
      o = 1;
    else
      o = 0;
  else if (type == 1 && x.type == 3)
    if ((PtRs->low < x.PtBs && PtRs->high > x.PtBs)
        &&
        (PtRt->low >= x.PtRt->low && PtRt->low < x.PtRt->high
         ||
         PtRt->high > x.PtRt->low && PtRt->high <= x.PtRt->high
         ||
         PtRt->low <= x.PtRt->low && PtRt->high > x.PtRt->low
         ||
         PtRt->low < x.PtRt->high && PtRt->high >= x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 1 && x.type == 4)
    if ((PtRs->low < x.PtBs && PtRs->high > x.PtBs)
        &&
        (PtRt->low < x.PtBt && PtRt->high > x.PtBt))
      o = 1;
    else
      o = 0;
  else if (type == 2 && x.type == 1)
    if ((PtRs->low >= x.PtRs->low && PtRs->low < x.PtRs->high
         ||
         PtRs->high > x.PtRs->low && PtRs->high <= x.PtRs->high
         ||
         PtRs->low <= x.PtRs->low && PtRs->high > x.PtRs->low
         ||
         PtRs->low < x.PtRs->high && PtRs->high >= x.PtRs->high)
         &&
        (PtBt > x.PtRt->low && PtBt < x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 2 && x.type == 2)
    if ((PtRs->low >= x.PtRs->low && PtRs->low < x.PtRs->high
         ||
         PtRs->high > x.PtRs->low && PtRs->high <= x.PtRs->high
         ||
         PtRs->low <= x.PtRs->low && PtRs->high > x.PtRs->low
         ||
         PtRs->low < x.PtRs->high && PtRs->high >= x.PtRs->high)
         &&
        (PtBt == x.PtBt))
      o = 1;
    else
      o = 0;
  else if (type == 2 && x.type == 3)
    if ((PtRs->low < x.PtBs && PtRs->high > x.PtBs)
        &&
        (PtBt > x.PtRt->low && PtBt < x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 2 && x.type == 4)
    if ((PtRs->low < x.PtBs && PtRs->high > x.PtBs)
        &&
        (PtBt == x.PtBt))
      o = 1;
    else
      o = 0;
  else if (type == 3 && x.type == 1)
    if ((PtBs > x.PtRs->low && PtBs < x.PtRs->high)
        &&
        (PtRt->low >= x.PtRt->low && PtRt->low < x.PtRt->high
         ||
         PtRt->high > x.PtRt->low && PtRt->high <= x.PtRt->high
         ||
         PtRt->low <= x.PtRt->low && PtRt->high > x.PtRt->low
         ||
         PtRt->low < x.PtRt->high && PtRt->high >= x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 3 && x.type == 2)
    if ((PtBs > x.PtRs->low && PtBs < x.PtRs->high)
        &&
        (PtRt->low < x.PtBt && PtRt->high > x.PtBt))
      o = 1;
    else
      o = 0;
  else if (type == 3 && x.type == 3)
    if ((PtBs == x.PtBs)
        &&
        (PtRt->low >= x.PtRt->low && PtRt->low < x.PtRt->high
         ||
         PtRt->high > x.PtRt->low && PtRt->high <= x.PtRt->high
         ||
         PtRt->low <= x.PtRt->low && PtRt->high > x.PtRt->low
         ||
         PtRt->low < x.PtRt->high && PtRt->high >= x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 3 && x.type == 4)
    if ((PtBs == x.PtBs)
        &&
        (PtRt->low < x.PtBt && PtRt->high > x.PtBt))
      o = 1;
    else
      o = 0;
  else if (type == 4 && x.type == 1)
    if ((PtBs > x.PtRs->low && PtBs < x.PtRs->high)
        &&
        (PtBt > x.PtRt->low && PtBt < x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 4 && x.type == 2)
    if ((PtBs > x.PtRs->low && PtBs < x.PtRs->high)
        &&
        (PtBt == x.PtBt))
      o = 1;
    else
      o = 0;
  else if (type == 4 && x.type == 3)
    if ((PtBs == x.PtBs)
        &&
        (PtBt > x.PtRt->low && PtBt < x.PtRt->high))
      o = 1;
    else
      o = 0;
  else if (type == 4 && x.type == 4)
    if ((PtBs == x.PtBs)
        &&
        (PtBt == x.PtBt))
      o = 1;
    else
      o = 0;
  
  return o;
}

//  int K_POINT2D :: equiv(const K_POINT2D& x) const
//    return 1 if *this and x are identical or
//                            have already been identified, and
//            0 otherwise

int K_POINT2D :: equiv(const K_POINT2D& x) const
{
  K_POINT2D* curr;
  int        e;
  
  if (this == &x)
    e = 1;
  else  //  if (this != &x)
    for (curr = this->next, e = 0; curr != this && !e; curr = curr->next)
      if (curr == &x)
        e = 1;
  
  return e;
}

//  int K_POINT2D :: equal(const K_POINT2D& x) const
//    return 1 if the boxes for *this and x are equal, and
//           0 otherwise

int K_POINT2D :: equal(K_POINT2D& x)
{
  assert(type);
  assert(x.type);
  
  int e;
  
  //  See if *this and x are identical or have already been identified.
  
  if (equiv(x))
    e = 1;
  else  //  if (!equiv(x)) => this != &x
  {
    //  Cut *this and x s.t. one is entirely contained in the other.
    
    cut_s(x.get_low_s());
    cut_s(x.get_high_s());
    cut_t(x.get_low_t());
    cut_t(x.get_high_t());
    x.cut_s(get_low_s());
    x.cut_s(get_high_s());
    x.cut_t(get_low_t());
    x.cut_t(get_high_t());
    
    if (type != x.type)
      e = 0;
    else if (!overlap(x))
      e = 0;
    else  //  if (type == x.type && overlap(x))
    {
      unsigned long i;
      unsigned long num_pts_s, num_pts_t;
      K_POINT1D**   pts_s;
      K_POINT1D**   pts_t;
      
      if (type == 1)  //  type  == 1 => x.type == 1
      {
//        assert(PtRs->low  == x.PtRs->low);
//        assert(PtRs->high == x.PtRs->high);
//        assert(PtRt->low  == x.PtRt->low);
//        assert(PtRt->high == x.PtRt->high);
        
        //  See if 2 univariate polynomials
        //           that define the s-coordinate of *this (and x)
        //        have 1 and only 1 common root
        //          in the interval for the s-coordinate of *this (and x).
        
        if ((num_pts_s = get_pts(PtRs->low, PtRs->high,
                                 gcd(*PtRs->poly, *x.PtRs->poly),
                                 pts_s,
                                 0, 0)) == 1)
        {
          //  See if 2 univariate polynomials
          //           that define the t-coordinate of *this (and x)
          //        have 1 and only 1 common root
          //          in the interval for the t-coordinate of *this (and x).
          
          if ((num_pts_t = get_pts(PtRt->low, PtRt->high,
                                   gcd(*PtRt->poly, *x.PtRt->poly),
                                   pts_t,
                                   0, 0)) == 1)
          {
            prev->next   = x.next;
            x.next->prev = prev;
            prev         = &x;
            x.next       = this;
            e            = 1;
          }
          else  //  if (num_pts_t == 0 || num_pts_t > 1)
            e = 0;
          
          for (i = 0; i < num_pts_t; i++)
            if (!--pts_t[i]->ref_count)
              delete pts_t[i];
          
          if (num_pts_t > 0)
            delete [] pts_t;
        }
        else  //  if (num_pts_s == 0 || num_pts_s > 1)
          e = 0;
        
        for (i = 0; i < num_pts_s; i++)
          if (!--pts_s[i]->ref_count)
            delete pts_s[i];
        
        if (num_pts_s > 0)
          delete [] pts_s;
      }
      else if (type == 2)  //  type  == 2 => x.type == 2
      {
//        assert(PtRs->low  == x.PtRs->low);
//        assert(PtRs->high == x.PtRs->high);
//        assert(PtBt       == x.PtBt);
        
        //  See if 2 univariate polynomials
        //           that define the s-coordinate of *this (and x)
        //        have 1 and only 1 common root
        //          in the interval for the s-coordinate of *this (and x).
        
        if ((num_pts_s = get_pts(PtRs->low, PtRs->high,
                                 gcd(*PtRs->poly, *x.PtRs->poly),
                                 pts_s,
                                 0, 0)) == 1)
        {
          prev->next   = x.next;
          x.next->prev = prev;
          prev         = &x;
          x.next       = this;
          e            = 1;
        }
        else  //  if (num_pts_s == 0 || num_pts_s > 1)
          e = 0;
        
        for (i = 0; i < num_pts_s; i++)
          if (!--pts_s[i]->ref_count)
            delete pts_s[i];
        
        if (num_pts_s > 0)
          delete [] pts_s;
      }
      else if (type == 3)  //  type == 3 => x.type == 3
      {
//        assert(PtBs       == x.PtBs);
//        assert(PtRt->low  == x.PtRt->low);
//        assert(PtRt->high == x.PtRt->high);
        
        //  See if 2 univariate polynomials
        //           that define the t-coordinate of *this (and x)
        //        have 1 and only 1 common root
        //          in the interval for the t-coordinate of *this (and x).
        
        if ((num_pts_t = get_pts(PtRt->low, PtRt->high,
                                 gcd(*PtRt->poly, *x.PtRt->poly),
                                 pts_t,
                                 0, 0)) == 1)
        {
          prev->next   = x.next;
          x.next->prev = prev;
          prev         = &x;
          x.next       = this;
          e            = 1;
        }
        else  //  if (num_pts_t == 0 || num_pts_t > 1)
          e = 0;
        
        for (i = 0; i < num_pts_t; i++)
          if (!--pts_t[i]->ref_count)
            delete pts_t[i];
        
        if (num_pts_t > 0)
          delete [] pts_t;
      }
      else  //  if (type == x.type == 4)
      {
//        assert(PtBs == x.PtBs);
//        assert(PtBt == x.PtBt);
        
        prev->next   = x.next;
        x.next->prev = prev;
        prev         = &x;
        x.next       = this;
        e            = 1;
      }
    }
  }
  
  return e;
}

//  int K_POINT2D :: separate(const K_POINT2D& x) const
//    separate *this and x s.t. they will not overlap.
//    POSSIBLY DOES NOT TERMINATE!

int K_POINT2D :: separate(const K_POINT2D& x) const
{
  assert(type);
  assert(x.type);
  
  if (type == 4 && x.type == 3)
    x.cut_t(PtBt);
  else if (type == 4 && x.type == 2)
    x.cut_s(PtBs);
  else if (type == 4 && x.type == 1)
  {
    x.cut_s(PtBs);
    x.cut_t(PtBt);
  }
  else if (type == 3 && x.type == 4)
    cut_t(x.PtBt);
  else if (type == 3 && x.type == 3)
  {
    cut_t(x.PtRt->low);
    cut_t(x.PtRt->high);
    x.cut_t(PtRt->low);
    x.cut_t(PtRt->high);
  }
  else if (type == 3 && x.type == 2)
  {
    cut_t(x.PtBt);
    x.cut_s(PtBs);
  }
  else if (type == 3 && x.type == 1)
  {
    cut_t(x.PtRt->low);
    cut_t(x.PtRt->high);
    x.cut_s(PtBs);
    x.cut_t(PtRt->low);
    x.cut_t(PtRt->high);
  }
  else if (type == 2 && x.type == 4)
    cut_s(x.PtBs);
  else if (type == 2 && x.type == 3)
  {
    cut_s(x.PtBs);
    x.cut_t(PtBt);
  }
  else if (type == 2 && x.type == 2)
  {
    cut_s(x.PtRs->low);
    cut_s(x.PtRs->high);
    x.cut_s(PtRs->low);
    x.cut_s(PtRs->high);
  }
  else if (type == 2 && x.type == 1)
  {
    cut_s(x.PtRs->low);
    cut_s(x.PtRs->high);
    x.cut_s(PtRs->low);
    x.cut_s(PtRs->high);
    x.cut_t(PtBt);
  }
  else if (type == 1 && x.type == 4)
  {
    cut_s(x.PtBs);
    cut_t(x.PtBt);
  }
  else if (type == 1 && x.type == 3)
  {
    cut_s(x.PtBs);
    cut_t(x.PtRt->low);
    cut_t(x.PtRt->high);
    x.cut_t(PtRt->low);
    x.cut_t(PtRt->high);
  }
  else if (type == 1 && x.type == 2)
  {
    cut_s(x.PtRs->low);
    cut_s(x.PtRs->high);
    cut_t(x.PtBs);
    x.cut_s(PtRs->low);
    x.cut_s(PtRs->high);
  }
  else if (type == 1 && x.type == 1)
  {
    cut_s(x.PtRs->low);
    cut_s(x.PtRs->high);
    cut_t(x.PtRt->low);
    cut_t(x.PtRt->high);
    x.cut_s(PtRs->low);
    x.cut_s(PtRs->high);
    x.cut_t(PtRt->low);
    x.cut_t(PtRt->high);
  }
  
  while (overlap(x))
  {
    shrink(shrink_step, shrink_step);
    x.shrink(shrink_step, shrink_step);
    separate(x);
  }
  
  return 0;
}

//  K_POINT2D :: equal_s(K_POINT2D& x)
//    returns 1 if *this and x are equal in the s-coordinate, and
//            0 otherwise.

int K_POINT2D :: equal_s(K_POINT2D& x)
{
  assert(type);
  assert(x.type);
  
  unsigned long i;
  K_POINT1D**   pts;
  unsigned long num_pts;
  int           e;
  
  //  See if *this and x are identical or have already been identified.
  
  if (equiv(x))
    e = 1;
  else  //  if (!equiv(x))
  {
    //  Cut *this and x s.t. one is entirely contained in the other.
    
    cut_s(x.get_low_s());
    cut_s(x.get_high_s());
    x.cut_s(get_low_s());
    x.cut_s(get_high_s());
    
    if (compare_s(x))
      e = 0;
    else if ((type == 1 || type == 2) && (x.type == 1 || x.type == 2))
    {
      //  See if 2 univariate polynomials
      //           that define the s-coordinate of *this (and x)
      //        have 1 and only 1 common root
      //          in the interval for the s-coordinate of *this (and x).
      
      if ((num_pts = get_pts(PtRs->low, PtRs->high,
                             gcd(*PtRs->poly, *x.PtRs->poly),
                             pts,
                             0, 0)) == 1)
        e = 1;
      else  //  if (num_pts == 0 || num_pts > 1)
        e = 0;
      
      for (i = 0; i < num_pts; i++)
        if (!--pts[i]->ref_count)
          delete pts[i];
      
      if (num_pts > 0)
        delete [] pts;
    }
    else if ((type == 3 || type == 4) && (x.type == 3 || x.type == 4))
      if (PtBs == x.PtBs)
        e = 1;
      else  //  if (PtBs != x.PtBs)
        e = 0;
    else  //  if ((type == 1 || type == 2) && (x.type == 3 || x.type == 4)
          //      ||
          //      (type == 3 || type == 4) && (x.type == 1 || x.type == 2))
      e = 0;
  }
  
  return e;
}

//  K_POINT2D :: equal_t(K_POINT2D& x)
//    returns 1 if *this and x are equal in the t-coordinate, and
//            0 otherwise.

int K_POINT2D :: equal_t(K_POINT2D& x)
{
  assert(type);
  assert(x.type);
  
  unsigned long i;
  K_POINT1D**   pts;
  unsigned long num_pts;
  int           e;
  
  //  See if *this and x are identical or have already been identified.
  
  if (equiv(x))
    e = 1;
  else  //  if (!equiv(x))
  {
    //  Cut *this and x s.t. one is entirely contained in the other.
    
    cut_t(x.get_low_t());
    cut_t(x.get_high_t());
    x.cut_t(get_low_t());
    x.cut_t(get_high_t());
    
    if (compare_t(x))
      e = 0;
    else if ((type == 1 || type == 3) && (x.type == 1 || x.type == 3))
    {
      //  See if 2 univariate polynomials
      //           that define the t-coordinate of *this (and x)
      //        have 1 and only 1 common root
      //          in the interval for the t-coordinate of *this (and x).
      
      if ((num_pts = get_pts(PtRt->low, PtRt->high,
                             gcd(*PtRt->poly, *x.PtRt->poly),
                             pts,
                             0, 0)) == 1)
        e = 1;
      else  //  if (num_pts == 0 || num_pts > 1)
        e = 0;
      
      for (i = 0; i < num_pts; i++)
        if (!--pts[i]->ref_count)
          delete pts[i];
      
      if (num_pts > 0)
        delete [] pts;
    }
    else if ((type == 2 || type == 4) && (x.type == 2 || x.type == 4))
      if (PtBt == x.PtBt)
        e = 1;
      else  //  if (PtBt != x.PtBt)
        e = 0;
    else  //  if ((type == 1 || type == 3) && (x.type == 2 || x.type == 4)
          //      ||
          //      (type == 2 || type == 4) && (x.type == 1 || x.type == 3))
      e = 0;
  }
  
  return e;
}

////  K_POINT2D :: equal_dir(K_POINT2D& x, const unsigned long i)
////    returns 1 if *this and x are equal in the i-th coordinate, and
////            0 otherwise.
//
//int K_POINT2D :: equal_dir(K_POINT2D& x, const unsigned long i)
//{
//  assert(type);
//  assert(x.type);
//  assert(i == 0 || i == 1);
//  
//  int e;
//  
//  //  See if *this and x are identical or have already been identified.
//  
//  if (equiv(x))
//    e = 1;
//  else  //  if (!equiv(x))
//    if (i == 0)
//      e = equal_s(x);
//    else  //  if (i == 1)
//      e = equal_t(x);
//  
//  return e;
//}

//  unsigned long get_all_pts(const K_RATPOLY& P1, const K_RATPOLY& P2,
//                            K_POINT2D**& pts,
//                            const bigrational& tol)
//    computes all the intersections "pts" of
//      2 bivariate polynomials P1 and P2.
//    returns the number of intersections.
//    if tol >  0 then
//      a box for each intersection is no larger than tol in any direction.
//    if tol = 0 then
//      there is no limit on the size of boxes for intersections.

unsigned long get_all_pts(const K_RATPOLY& P1, const K_RATPOLY& P2,
                          K_POINT2D**& pts,
                          const bigrational& tol)
{
  assert(P1.get_num_vars() == 2);
  assert(P2.get_num_vars() == 2);
  
  unsigned long num_int_pts;
  bigrational   l_s, h_s, l_t, h_t;
  
  h_s = GoodSylvester(P1, P2, 1).get_Mignotte_bd() + s_epsilon;
  l_s = - h_s - s_epsilon;
  h_t = GoodSylvester(P1, P2, 0).get_Mignotte_bd() + t_epsilon;
  l_t = - h_t - t_epsilon;
  
  return get_pts(l_s, h_s, l_t, h_t, P1, P2, pts, tol, 0);
}

int K_POINT2D :: get_fp_approx(double*& a) const
{
  if (type == 1 || type == 2)
  {
    if (PtRs->ok_float)
      a[0] = PtRs->float_est;
    else  //  if (!PtRs->ok_float)
      a[0] = ((PtRs->low + PtRs->high) / 2).as_double();
  }
  else  //  if (type == 3 || type == 4)
    a[0] = PtBs.as_double();
  
  if (type == 1 || type == 3)
  {
    if (PtRt->ok_float)
      a[1] = PtRt->float_est;
    else  //  if (!PtRt->ok_float)
      a[1] = ((PtRt->low + PtRt->high) / 2).as_double();
  }
  else  //  if (type == 2 || type == 4)
    a[1] = PtBt.as_double();
  
  return 0;
}

