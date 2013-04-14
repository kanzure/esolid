//  file:    ksurf.cc
//  update:  11/18/02

#include <config.h>

#include <cstdlib>
#include <ctime>

#include <ksurf.h>
#include <bigrational_matrix.h>
#include <pascal.h>

K_SURF :: K_SURF()
  : Impl(0), Impl_ok(0),
    X(0), Y(0), Z(0), W(0), mon_ok(0),
    ref_count(0)
{ }

K_SURF :: K_SURF(K_RATPOLY* const impl,
                 K_RATPOLY* const x, K_RATPOLY* const y, K_RATPOLY* const z,
                 K_RATPOLY* const w)
{
  if (Impl = impl)
  {
    Impl->ref_count++;
    Impl_ok = 1;
  }
  else
    Impl_ok = 0;
  
  if (x && y && z && w)
  {
    X      = x;
    X->ref_count++;
    Y      = y;
    Y->ref_count++;
    Z      = z;
    Z->ref_count++;
    W      = w;
    W->ref_count++;
    mon_ok = 1;
  }
  else
  {
    X = 0;
    Y = 0;
    Z = 0;
    W = 0;
    mon_ok = 0;
  }
  
  ref_count = 0;
}

K_SURF :: K_SURF(const K_RATPOLY& impl)
{
  Impl    = new K_RATPOLY(impl);
  Impl->ref_count++;
  Impl_ok = 1;
  
  X      = 0;
  Y      = 0;
  Z      = 0;
  W      = 0;
  mon_ok = 0;
  
  ref_count = 0;
}

K_SURF :: K_SURF(const K_SURF& s)
{
  if (Impl_ok = s.Impl_ok)
  {
    Impl = s.Impl;
    Impl->ref_count++;
  }
  else
    Impl = 0;
  
  if (mon_ok = s.mon_ok)
  {
    X = s.X;
    X->ref_count++;
    Y = s.Y;
    Y->ref_count++;
    Z = s.Z;
    Z->ref_count++;
    W = s.W;
    W->ref_count++;
  }
  else
    X = Y = Z = W = 0;
  
  ref_count = 0;
}

K_SURF& K_SURF :: operator =(const K_SURF& s)
{
  if (this != &s)
  {
    if (Impl && !--Impl->ref_count)
      delete Impl;
    
    if (X && !--X->ref_count)
      delete X;
    
    if (Y && !--Y->ref_count)
      delete Y;
    
    if (Z && !--Z->ref_count)
      delete Z;
    
    if (W && !--W->ref_count)
      delete W;
    
    if (Impl_ok = s.Impl_ok)
    {
      Impl = s.Impl;
      Impl->ref_count++;
    }
    else
      Impl = 0;
    
    if (mon_ok = s.mon_ok)
    {
      X = s.X;
      X->ref_count++;
      Y = s.Y;
      Y->ref_count++;
      Z = s.Z;
      Z->ref_count++;
      W = s.W;
      W->ref_count++;
    }
    else
      X = Y = Z = W = 0;
  }
  
  return *this;
}

K_SURF :: ~K_SURF()
{
  if (Impl && !--Impl->ref_count)
    delete Impl;
  
  if (X && !--X->ref_count)
    delete X;
  
  if (Y && !--Y->ref_count)
    delete Y;
  
  if (Z && !--Z->ref_count)
    delete Z;
  
  if (W && !--W->ref_count)
    delete W;  
}

ostream& K_SURF :: output(ostream& o) const
{
  if (Impl_ok)
    o << " implicit: " << endl << *Impl << flush;
  
  if (mon_ok)
  {
    o << " parametrized: X = " << endl << *X << endl << flush;
    o << " parametrized: Y = " << endl << *Y << endl << flush;
    o << " parametrized: Z = " << endl << *Z << endl << flush;
    o << " parametrized: W = " << endl << *W << flush;
  }
  
  if (!Impl_ok && !mon_ok)
    o << " NULL " << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const K_SURF& s)
{
  return s.output(o);
}

//  int K_SURF :: param_to_coord(const bigrational& s,
//                             const bigrational& t,
//                             bigrational& x,
//                             bigrational& y,
//                             bigrational& z) const
//    compute
//      (x, y, z) = (X(s, t) / W(s, t), Y(s, t) / W(s, t), Z(s, t) / W(s, t)).

int K_SURF :: param_to_coord(const bigrational& s,
                             const bigrational& t,
                             bigrational& x,
                             bigrational& y,
                             bigrational& z) const
{
  assert(mon_ok);
  assert(X);
  assert(Y);
  assert(Z);
  assert(W);
  
  bigrational p[2];
  bigrational w;
  
  p[0] = s;
  p[1] = t;
  assert(sgn(w = W->evaluate(p)));
  x = X->evaluate(p) / w;
  y = Y->evaluate(p) / w;
  z = Z->evaluate(p) / w;
  
  return 0;
}

//  int K_SURF :: param_to_coord(const bigrational_vector& p,
//                               bigrational_vector& c) const
//    compute
//      c = (X(p) / W(p), Y(p) / W(p), Z(p) / W(p)).

int K_SURF :: param_to_coord(const bigrational_vector& p,
                             bigrational_vector& c) const
{
  assert(mon_ok);
  assert(X);
  assert(Y);
  assert(Z);
  assert(W);
  assert(p.get_dim() == 2);
  assert(c.get_dim() == 3);
  
  bigrational w;
  
  assert(sgn(w = W->evaluate(p)));
  c[0] = X->evaluate(p) / w;
  c[1] = Y->evaluate(p) / w;
  c[2] = Z->evaluate(p) / w;
  
  return 0;
}

K_BOX3D K_SURF :: get_range(const bigrational& l_s,
                            const bigrational& h_s,
                            const bigrational& l_t,
                            const bigrational& h_t) const
{
  assert(mon_ok);
  
  unsigned long i;
  K_BOX3D       b;
  bigrational   v_l[2];
  bigrational   v_h[2];
  bigrational   X_l, X_h, Y_l, Y_h, Z_l, Z_h, W_l, W_h;
  
  v_l[0] = l_s;
  v_l[1] = l_t;
  v_h[0] = h_s;
  v_h[1] = h_t;
  
  X->eval_range(v_l, v_h, X_l, X_h);
  Y->eval_range(v_l, v_h, Y_l, Y_h);
  Z->eval_range(v_l, v_h, Z_l, Z_h);
  W->eval_range(v_l, v_h, W_l, W_h);
  
  if (sgn(W_h) < 0)  //  if (W_l <= W_h < 0)
  {
    if (sgn(X_l) > 0)  //  if (0 < X_l <= X_h)
    {
      b.low[0]  = X_h / W_h;
      b.high[0] = X_l / W_l;
    }
    else if (sgn(X_h) < 0)  //  if (X_l <= X_h < 0)
    {
      b.low[0]  = X_h / W_l;
      b.high[0] = X_l / W_h;
    }
    else  //  if (X_l <= 0 <= X_h)
    {
      b.low[0]  = X_h / W_h;
      b.high[0] = X_l / W_h;
    }
    
    if (sgn(Y_l) > 0)  //  if (0 < Y_l <= Y_h)
    {
      b.low[1]  = Y_h / W_h;
      b.high[1] = Y_l / W_l;
    }
    else if (sgn(Y_h) < 0)  //  if (Y_l <= Y_h < 0)
    {
      b.low[1]  = Y_h / W_l;
      b.high[1] = Y_l / W_h;
    }
    else  //  if (Y_l <= 0 <= Y_h)
    {
      b.low[1]  = Y_h / W_h;
      b.high[1] = Y_l / W_h;
    }
    
    if (sgn(Z_l) > 0)  //  if (0 < Z_l <= Z_h)
    {
      b.low[2]  = Z_h / W_h;
      b.high[2] = Z_l / W_l;
    }
    else if (sgn(Z_h) < 0)  //  if (Z_l <= Z_h < 0)
    {
      b.low[2]  = Z_h / W_l;
      b.high[2] = Z_l / W_h;
    }
    else  //  if (Z_l <= 0 <= Z_h)
    {
      b.low[2]  = Z_h / W_h;
      b.high[2] = Z_l / W_h;
    }
    
    for (i = 0; i < 3; i++)
      b.low_infty[i] = b.high_infty[i] = 0;
  }
  else if (sgn(W_l) < 0 && !sgn(W_h))
  {
    if (sgn(X_l) > 0)  //  if (0 < X_l <= X_h)
    {
      b.low[0]       = X_h / W_l;
      b.low_infty[0] = 0;
    }
    else if (sgn(X_h) < 0)  //  if (X_l <= X_h < 0)
    {
      b.high[0]       = X_l / W_l;
      b.high_infty[0] = 0;
    }
    
    if (sgn(Y_l) > 0)  //  if (0 < Y_l <= Y_h)
    {
      b.low[1]       = Y_h / W_l;
      b.low_infty[1] = 0;
    }
    else if (sgn(Y_h) < 0)  //  if (Y_l <= Y_h < 0)
    {
      b.high[1]       = Y_l / W_l;
      b.high_infty[1] = 0;
    }
    
    if (sgn(Z_l) > 0)  //  if (0 < Z_l <= Z_h)
    {
      b.low[2]       = Z_h / W_l;
      b.low_infty[2] = 0;
    }
    else if (sgn(Z_h) < 0)  //  if (Z_l <= Z_h < 0)
    {
      b.high[2]       = Z_l / W_l;
      b.high_infty[2] = 0;
    }
  }
  else if (!sgn(W_l) && !sgn(W_h))
  {
    if (sgn(X_l) > 0)  //  if (0 < X_l <= X_h)
      b.low_infty[0] = 1;
    else if (sgn(X_h) < 0)  //  if (X_l <= X_h < 0)
      b.high_infty[0] = - 1;
    
    if (sgn(Y_l) > 0)  //  if (0 < Y_l <= Y_h)
      b.low_infty[1] = 1;
    else if (sgn(Y_h) < 0)  //  if (Y_l <= Y_h < 0)
      b.high_infty[1] = - 1;
    
    if (sgn(Z_l) > 0)  //  if (0 < Z_l <= Z_h)
      b.low_infty[2] = 1;
    else if (sgn(Z_h) < 0)  //  if (Z_l <= Z_h < 0)
      b.high_infty[2] = - 1;
  }
  else if (!sgn(W_l) && sgn(W_h) > 0)
  {
    if (sgn(X_l) > 0)  //  if (0 < X_l <= X_h)
    {
      b.low[0]       = X_l / W_h;
      b.low_infty[0] = 0;
    }
    else if (sgn(X_h) < 0)  //  if (X_l <= X_h < 0)
    {
      b.high[0]       = X_h / W_h;
      b.high_infty[0] = 0;
    }
    
    if (sgn(Y_l) > 0)  //  if (0 < Y_l <= Y_h)
    {
      b.low[1]       = Y_l / W_h;
      b.low_infty[1] = 0;
    }
    else if (sgn(Y_h) < 0)  //  if (Y_l <= Y_h < 0)
    {
      b.high[1]       = Y_h / W_h;
      b.high_infty[1] = 0;
    }
    
    if (sgn(Z_l) > 0)  //  if (0 < Z_l <= Z_h)
    {
      b.low[2]       = Z_l / W_h;
      b.low_infty[2] = 0;
    }
    else if (sgn(Z_h) < 0)  //  if (Z_l <= Z_h < 0)
    {
      b.high[2]       = Z_h / W_h;
      b.high_infty[2] = 0;
    }
  }
  else if (sgn(W_l) > 0)  //  if (0 < W_l <= W_h)
  {
    if (sgn(X_l) > 0)  //  if (0 < X_l <= X_h)
    {
      b.low[0]  = X_l / W_h;
      b.high[0] = X_h / W_l;
    }
    else if (sgn(X_h) < 0)  //  if (X_l <= X_h < 0)
    {
      b.low[0]  = X_l / W_l;
      b.high[0] = X_h / W_h;
    }
    else  //  if (X_l <= 0 <= X_h)
    {
      b.low[0]  = X_l / W_l;
      b.high[0] = X_h / W_l;
    }
    
    if (sgn(Y_l) > 0)  //  if (0 < Y_l <= Y_h)
    {
      b.low[1]  = Y_l / W_h;
      b.high[1] = Y_h / W_l;
    }
    else if (sgn(Y_h) < 0)  //  if (Y_l <= Y_h < 0)
    {
      b.low[1]  = Y_l / W_l;
      b.high[1] = Y_h / W_h;
    }
    else  //  if (Y_l <= 0 <= Y_h)
    {
      b.low[1]  = Y_l / W_l;
      b.high[1] = Y_h / W_l;
    }
    
    if (sgn(Z_l) > 0)  //  if (0 < Z_l <= Z_h)
    {
      b.low[2]  = Z_l / W_h;
      b.high[2] = Z_h / W_l;
    }
    else if (sgn(Z_h) < 0)  //  if (Z_l <= Z_h < 0)
    {
      b.low[2]  = Z_l / W_l;
      b.high[2] = Z_h / W_h;
    }
    else  //  if (Z_l <= 0 <= Z_h)
    {
      b.low[2]  = Z_l / W_l;
      b.high[2] = Z_h / W_l;
    }
    
    for (i = 0; i < 3; i++)
      b.low_infty[i] = b.high_infty[i] = 0;
  }
  
  return b;
}

int match_pts(K_POINT2D** const pts1, const unsigned long num_pts1,
              K_SURF* const s1,
              K_POINT2D** const pts2, const unsigned long num_pts2,
              K_SURF* const s2)
{
  assert(num_pts1 <= num_pts2);
  
  unsigned long i, j;
  K_BOX3D*      boxes1;
  K_BOX3D*      boxes2;
  unsigned long num_match;
  unsigned long num_new_match, last_match;
  K_POINT2D*    p1;
  K_BOX3D       b1;
  K_POINT2D*    p2;
  K_BOX3D       b2;
  
  if (num_pts1 > 0)
  {
    boxes1 = new K_BOX3D [num_pts1];
    
    for (i = 0; i < num_pts1; i++)
      boxes1[i] = s1->get_range(pts1[i]->get_low_s(), pts1[i]->get_high_s(),
                                pts1[i]->get_low_t(), pts1[i]->get_high_t());
  }
  else  //  if (!num_pts1)
    boxes1 = 0;
  
  if (num_pts2 > 0)
  {
    boxes2 = new K_BOX3D [num_pts2];
    
    for (i = 0; i < num_pts2; i++)
      boxes2[i] = s2->get_range(pts2[i]->get_low_s(), pts2[i]->get_high_s(),
                                pts2[i]->get_low_t(), pts2[i]->get_high_t());
  }
  else  //  if(!num_pts2)
    boxes2 = 0;
  
//  cerr << " ksurf: match_pts: -------------------- " << endl << flush;
//  for (i = 0; i < num_pts1; i++)
//    cerr << " ksurf: match_pts: pts1[" << i << "] = " << *pts1[i] << ", boxes1[" << i << "] = " << boxes1[i] << endl << flush;
//  for (i = 0; i < num_pts2; i++)
//    cerr << " ksurf: match_pts: pts2[" << i << "] = " << *pts2[i] << ", boxes2[" << i << "] = " << boxes2[i] << endl << flush;
//  cerr << " ksurf: match_pts: 1: -------------------- " << endl << flush;
  
  num_match = 0;
  
  while (num_match < num_pts1)
  {
    //  Assume pts1[0], ..., pts1[num_match - 1] and
    //         pts2[0], ..., pts2[num_match - 1] are matched, respectively.
    
    for (i = num_match; i < num_pts1; i++)
    {
      num_new_match = 0;
      
      for (j = num_match; j < num_pts2; j++)
        if (boxes1[i].overlap(boxes2[j]))
        {
          last_match = j;
          num_new_match++;
        }
      
      assert(num_new_match > 0);
      
      //  If precisely one of
      //    pts2[num_match], ..., pts2[num_pts2 - 1]
      //  matches with pts1[i], say pts2[last_match], then exchange
      //    pts1[num_match]   & pts1[i],
      //    boxes1[num_match] & boxes1[i],
      //    pts2[num_match]   & pts2[last_match], and
      //    boxes2[num_match] & boxes2[last_match].
      
      if (num_new_match == 1)
      {
        if (i > num_match)
        {
          p1                = pts1[i];
          b1                = boxes1[i];
          pts1[i]           = pts1[num_match];
          boxes1[i]         = boxes1[num_match];
          pts1[num_match]   = p1;
          boxes1[num_match] = b1;
        }
        
        if (last_match > num_match)
        {
          p2                 = pts2[last_match];
          b2                 = boxes2[last_match];
          pts2[last_match]   = pts2[num_match];
          boxes2[last_match] = boxes2[num_match];
          pts2[num_match]    = p2;
          boxes2[num_match]  = b2;
        }
        
        num_match++;
      }
      
//      for (unsigned long ii = 0; ii < num_match; ii++)
//      {
//        cerr << " ksurf: match_pts: pts1[" << ii << "] = " << *pts1[ii] << ", boxes1[" << ii << "] = " << boxes1[ii] << endl << flush;
//        cerr << " ksurf: match_pts: pts2[" << ii << "] = " << *pts2[ii] << ", boxes2[" << ii << "] = " << boxes2[ii] << endl << flush;
//      }
//      cerr << " ksurf: match_pts: 2: -------------------- " << endl << flush;
    }
    
    //  Shrink pts1[num_match], ..., pts1[num_pts1 - 1] &
    //         pts2[num_match], ..., pts2[num_pts2 - 1].
    
    for (i = num_match; i < num_pts1; i++)
    {
      pts1[i]->shrink(shrink_step, shrink_step);
      boxes1[i] = s1->get_range(pts1[i]->get_low_s(), pts1[i]->get_high_s(),
                                pts1[i]->get_low_t(), pts1[i]->get_high_t());
    }
    
    if (num_match < num_pts1)
      for (i = num_match; i < num_pts2; i++)
      {
        pts2[i]->shrink(shrink_step, shrink_step);
        boxes2[i] = s2->get_range(pts2[i]->get_low_s(), pts2[i]->get_high_s(),
                                  pts2[i]->get_low_t(), pts2[i]->get_high_t());
      }
  }
  
//  for (i = 0; i < num_match; i++)
//  {
//    cerr << " ksurf: match_pts: pts1[" << i << "] = " << *pts1[i] << ", boxes1[" << i << "] = " << boxes1[i] << endl << flush;
//    cerr << " ksurf: match_pts: pts2[" << i << "] = " << *pts2[i] << ", boxes2[" << i << "] = " << boxes2[i] << endl << flush;
//  }
//  cerr << " ksurf: match_pts: 3: -------------------- " << endl << flush;
  
  if (num_pts1 > 0)
    delete [] boxes1;
  
  if (num_pts2 > 0)
    delete [] boxes2;
  
  return 0;
}

//  unsigned long get_plane_coeffs(const bigrational_vector v0,
//                                 const bigrational_vector v1,
//                                 const bigrational_vector v2,
//                                 bigrational*& p)
//    computes the plane p[0] X + p[1] Y + p[2] Z + p[3]
//      that passes 3 pts v0, v1 and v2.
//    returns 4
//            0 otherwise (i.e. pts v0, v1 and v2 are collinear.)

unsigned long get_plane_coeffs(const bigrational_vector& v0,
                               const bigrational_vector& v1,
                               const bigrational_vector& v2,
                               bigrational*& p)
{
  assert(v0.get_dim() == 3);
  assert(v1.get_dim() == 3);
  assert(v2.get_dim() == 3);
  
  unsigned long      i;
  bigrational_vector V1(3);
  bigrational_vector V2(3);
  bigrational_vector p_proto(3);
  unsigned long      r;
  
  for (i = 0; i < 3; i++)
  {
    V1[i] = v1[i] - v0[i];
    V2[i] = v2[i] - v0[i];
  }
  
  p_proto[0] = V1[1] * V2[2] - V1[2] * V2[1];
  p_proto[1] = V1[2] * V2[0] - V1[0] * V2[2];
  p_proto[2] = V1[0] * V2[1] - V1[1] * V2[0];
  
  if (sgn(p_proto[0]) || sgn(p_proto[1]) || sgn(p_proto[2]))
  {
    p = new bigrational [r = 4];
    
    for (i = 0; i < 3; i++)
      p[i] = p_proto[i];
    
    p[3] = - p[0] * v0[0] - p[1] * v0[1] - p[2] * v0[2];
  }
  else  //  if v0, v1 and v2 are collinear
  {
    p = 0;
    r = 0;
  }
  
  return r;
}

static Pascal pascal;

//  K_RATPOLY implicitize(K_RATPOLY& X, K_RATPOLY& Y, K_RATPOLY& Z,
//                        K_RATPOLY& W,
//                        long total_deg)
//     returns the implicit representation of the surface
//               parameterized by (X/W, Y/W, Z/W) and
//               of total degree at most n.

int implicitize(const K_RATPOLY& X, const K_RATPOLY& Y, const K_RATPOLY& Z,
                const K_RATPOLY& W,
                const long total_deg,
                K_RATPOLY*& I)
{
  unsigned long      i, j, k;
  time_t             tm;
  unsigned long      sd;
  long               full_ord;
  unsigned long      num_good_pts;
  bigrational_matrix good_pts;
  bigrational        v[2];
  bigrational        x0, y0, z0, w0;
  unsigned long      deg, ord;
  bigrational        x, y, z, w;
  unsigned long      r, c, piv;
  bigrational_vector r_v, ans, rhs;
  unsigned long      num_non_zero;
  long               d[3];
  long               p[3];
  bigrational        test_v[3];
  int                found;
  
//  cerr << " ksurf: implicitize: X = " << endl << X << endl << flush;
//  cerr << " ksurf: implicitize: Y = " << endl << Y << endl << flush;
//  cerr << " ksurf: implicitize: Z = " << endl << Z << endl << flush;
//  cerr << " ksurf: implicitize: W = " << endl << W << endl << flush;
//  cerr << " ksurf: implicitize: -------------------- " << endl << flush;
  
  time(&tm);
  sd = tm;
  srand(sd);
  
  full_ord     = pascal.get_Pascal(total_deg, 3);
  good_pts     = bigrational_matrix(full_ord, full_ord);
  num_good_pts = 0;
  
  while (num_good_pts < full_ord)
  {
    w0 = 0;
    
    while (!sgn(w0))
    {
      v[0] = bigrational(1, rand());
      v[1] = bigrational(1, rand());
      w0   = W.evaluate(v);
    }
    
    x0 = X.evaluate(v) / w0;
    y0 = Y.evaluate(v) / w0;
    z0 = Z.evaluate(v) / w0;
    
    for (deg = ord = 0; deg <= total_deg; deg++)
    {
      x = 1;
      
      for (i = 0; i <= deg; i++)
      {
        x *= x0;
        y  = 1;
        
        for (j = 0; j <= deg - i; j++)
        {
          y *= y0;
          z  = 1;
          
          for (k = 0; k <= deg - i - j; k++)
            z *= z0;
          
          good_pts(num_good_pts, ord++) = x * y * z;
        }
      }
    }
    
    num_good_pts++;
  }
  
//  cerr << " ksurf: implcitize: num_good_pts = " << num_good_pts << endl << flush;
//  cerr << " ksurf: implcitize: 0: good_pts = " << endl << good_pts << endl << flush;
//  cerr << " ksurf: implicitize: -------------------- " << endl << flush;
  
  //  Gaussian elimination.
  
  r_v = bigrational_vector(full_ord);
  ans = bigrational_vector(full_ord);
  rhs = bigrational_vector(full_ord);
  
  for (r = c = 0; c < full_ord; c++)
  {
    for (piv = r; piv < full_ord && !sgn(good_pts(piv, c)); piv++)
      ;
    
    if (piv < full_ord)
    {
      if (piv != r)
        for (j = 0; j < full_ord; j++)
        {
          r_v[j]           = good_pts(r, j);
          good_pts(r, j)   = good_pts(piv, j);
          good_pts(piv, j) = r_v[j];
        }
      
      for (j = c + 1; j < full_ord; j++)
        good_pts(r, j) /= good_pts(r, c);
      
      good_pts(piv, c) = 1;
      
      for (i = 0; i < full_ord; i++)
        if (i != r)
        {
          for (j = c + 1; j < full_ord; j++)
            good_pts(i, j) -= good_pts(i, c) * good_pts(r, j);
          
          good_pts(i, c) = 0;
        }
      
      r++;
    }
  }
  
//  cerr << " ksurf: implcitize: 1: good_pts = " << endl << good_pts << endl << flush;
//  cerr << " ksurf: implicitize: -------------------- " << endl << flush;
  
  for (i = 0; i < full_ord; i++)
  {
    ans[i] = 0;
    rhs[i] = 0;
  }
  
  r            = 0;
  num_non_zero = 0;
  
//  cerr << " ksurf: implicitize: r = " << r << ", full_ord = "<< full_ord << ", num_non_zero = " << num_non_zero << endl << flush;
//  cerr << " ksurf: implicitize: -------------------- " << endl << flush;
  
  while (r < full_ord && !num_non_zero)
  {
    //  See whether or not there are multiple non-zero elements in row r.
    
    for (j = num_non_zero = 0; j < full_ord; j++)
      if (sgn(good_pts(r, j)))
        num_non_zero++;
    
//    cerr << " ksurf: implcitize: num_non_zero = " << num_non_zero << endl << flush;
    
    if (num_non_zero > 1)
    {
      //  Take care of the first non-zero element in row r.
      
      for (j = 0; !sgn(good_pts(r, j)); j++)
        ;
      
      if (good_pts(r, j) != 1)
        for (k = j + 1; k < full_ord; k++)
          good_pts(r, k) /= good_pts(r, j);
      
      good_pts(r, j) = 0;
      ans[j]         = 1;
      rhs[j]         = - 1;
      
      for (i = r + 1; i < full_ord; i++)
        if (sgn(good_pts(i, j)))
        {
          rhs[i]         -= good_pts(i, j);
          good_pts(i, j)  = 0;
        }
      
      //  Take care of the second non-zero element in row r.
      
      for (++j; !sgn(good_pts(r, j)); j++)
        ;
      
      ans[j] = - 1 / good_pts(r, j);
      
      for (i = 0; i < full_ord; i++)
      {
        rhs[i]         -= ans[j] * good_pts(i, j);
        good_pts(i, j)  = 0;
      }
      
      for (++j; j < full_ord; j++)
        if (sgn(good_pts(r, j)))
          for (i = 0; i < full_ord; i++)
            good_pts(i, j) = 0;
    }
    else  //  if (num_non_zero <= 1)
    {
      if (num_non_zero == 1)
      {
        for (j = 0; sgn(good_pts(r, j)); j++)
          ;
        
        for (i = 0; i < full_ord; i++)
          good_pts(i, j) = 0;
        
        num_non_zero = 0;        
      }
      
      r++;
    }
    
//    cerr << " ksurf: implicitize: r = " << r << ", full_ord = "<< full_ord << ", num_non_zero = " << num_non_zero << endl << flush;
//    cerr << " ksurf: implicitize: -------------------- " << endl << flush;
  }
  
//  cerr << " ksurf: implcitize: 2: good_pts = " << endl << good_pts << endl << flush;
//  cerr << " ksurf: implicitize: -------------------- " << endl << flush;
  
  if (r < full_ord)
  {
    for (++r; r < full_ord; r++)
    {
      for (j = 0; j < full_ord && !sgn(good_pts(r, j)); j++)
        ;
      
      if (j < full_ord)
      {
        ans[j] = rhs[r] / good_pts(r, j);
        
        for (i = r; i < full_ord; i++)
          if (sgn(good_pts(i, j)))
          {
            rhs[i]         -= ans[j] * good_pts(i, j);
            good_pts(i, j)  = 0;
          }
      }
      
      for (++j; j < full_ord; j++)
        if (sgn(good_pts(r, j)))
          for (i = r; i < full_ord; i++)
            good_pts(i, j) = 0;
    }
    
//    cerr << " ksurf: implcitize: 3: good_pts = " << endl << good_pts << endl << flush;
//    cerr << " ksurf: implcitize: 3: ans = " << endl << ans << endl << flush;
//    cerr << " ksurf: implicitize: -------------------- " << endl << flush;
    
    d[0] = d[1] = d[2] = total_deg;
    I    = new K_RATPOLY(3, d);
    
    for (deg = ord = 0; deg <= total_deg; deg++)
      for (i = 0; i <= deg; i++)
      {
        p[0] = i;
        
        for (j = 0; j <= deg - i; j++)
        {
          p[1] = j;
          p[2] = deg - i - j;
          
          I->get_coeff(p) = ans[ord++];
        }
      }
    
    I->reduce_deg();
    
//    cerr << " ksurf: implcitize: *I = " << endl << *I << endl << flush;
//    cerr << " ksurf: implicitize: -------------------- " << endl << flush;
    
    //  Test whether or not I is correct.
    
    if (!(I->deg[0] + I->deg[1] + I->deg[2]) && !sgn(I->coeffs[0]))
    {
//      cerr << " ksurf: implicitize: retry: 0 " << endl << flush;
//      cerr << " ksurf: implicitize: -------------------- " << endl << flush;
      delete I;
      I     = 0;
      found = 0;
    }
    else  //  if *I is not zero
    {
      //  Double check *I vanishes at (X/W, Y/W, Z/W)|v.
      
      w0 = 0;
      
      while (!sgn(w0))
      {
        v[0] = bigrational(1, rand());
        v[1] = bigrational(1, rand());
        w0   = W.evaluate(v);
      }
      
      test_v[0] = X.evaluate(v) / w0;
      test_v[1] = Y.evaluate(v) / w0;
      test_v[2] = Z.evaluate(v) / w0;
//      cerr << " ksurf: implicitize: test_v = (" << test_v[0] << ", " << test_v[1] << ", " << test_v[2] << ") " << endl << flush;
//      cerr << " kusrf: implicitize: *I = " << endl << *I << endl << flush;
//      cerr << " ksurf: implicitize: I->evaluate(test_v) = " << I->evaluate(test_v) << endl << flush;
//      cerr << " ksurf: implicitize: -------------------- " << endl << flush;
      
      if (sgn(I->evaluate(test_v)))
      {
//        cerr << " ksurf: implicitize: retry 1: " << endl << flush;
//        cerr << " ksurf: implicitize: -------------------- " << endl << flush;
        delete I;
        I     = 0;
        found = 0;
      }
      else  //  if (!sgn(I->evaluate(test_v)))
        found = 1;
    }
  }
  else  //  if (num_non_zero)
  {
//    cerr << " ksurf: implicitize: retry 2: " << endl << flush;
//    cerr << " ksurf: implicitize: -------------------- " << endl << flush;
    I     = 0;
    found = 0;
  }
  
  return found;
}

//  K_SURF K_SURF :: split_surf(const bigrational b,
//                              const unsigned long i) const
//     returns a surface that intersects with *this
//       at the curve parametrized by (X_b/W_b, Y_b/W_b, Z_b/W_b)
//       where X_b is obtained by substituting b to the i-th variable of X...

K_SURF K_SURF :: split_surf(const bigrational& b, const unsigned long i) const
{
  assert(i == 0 || i == 1);
  
  unsigned long      j;
  K_RATPOLY          X_b, Y_b, Z_b, W_b;
  K_RATPOLY          QX_b, RX_b, QY_b, RY_b, QZ_b, RZ_b;
  static bigint      cd;
  time_t             tm;
  unsigned long      sd;
  int                found;
  bigrational        W_b_c;
  bigrational        c;
  bigrational_vector x(3);
  bigrational_vector y(3);
  bigrational_vector n(3);
  bigrational_vector z(3);
  bigrational*       plane;
  long               d[3];
  long               p[3];
  K_RATPOLY          Plane;
  unsigned long      t;
  K_RATPOLY          dI_dx, dI_dy, dI_dz;
  K_RATPOLY          T;
  K_RATPOLY          X2_b, Y2_b, Z2_b, W2_b;
  K_RATPOLY          X2, Y2, Z2, W2;
  K_RATPOLY*         I;
  
//  cerr << " ksurf: split_surf: *Impl = " << endl << *Impl << endl << flush;
//  cerr << " ksurf: split_surf: *X = " << endl << *X << endl << flush;
//  cerr << " ksurf: split_surf: *Y = " << endl << *Y << endl << flush;
//  cerr << " ksurf: split_surf: *Z = " << endl << *Z << endl << flush;
//  cerr << " ksurf: split_surf: *W = " << endl << *W << endl << flush;
//  cerr << " ksurf: split_surf: -------------------- " << endl << flush;
  
  cd    = 1;
  time(&tm);
  sd    = tm;
  srand(sd);
  found = 0;
  
  //  1. Compute the curve parametrized by (X0_b/W0_b, Y0_b/W0_b, Z0_b/W0_b)
  //       where X0_b, Y0_b, Z0_b and W0_b are obtained
  //             by substituting b to the i-th variable of X, Y, Z and W, resp.
  
  X_b = X->subst_val(i, b);
  Y_b = Y->subst_val(i, b);
  Z_b = Z->subst_val(i, b);
  W_b = W->subst_val(i, b);
  
//  cerr << " ksurf: split_surf: i = " << i << ", b = " << endl << b << endl << flush;
//  cerr << " ksurf: split_surf: X_b = " << endl << X_b << endl << flush;
//  cerr << " ksurf: split_surf: Y_b = " << endl << Y_b << endl << flush;
//  cerr << " ksurf: split_surf: Z_b = " << endl << Z_b << endl << flush;
//  cerr << " ksurf: split_surf: W_b = " << endl << W_b << endl << flush;
//  cerr << " ksurf: split_surf: -------------------- " << endl << flush;
  
  QX_b = div(X_b, W_b, RX_b);
  QY_b = div(Y_b, W_b, RY_b);
  QZ_b = div(Z_b, W_b, RZ_b);
  
//  cerr << " ksurf: split_surf: QX_b = " << endl << QX_b << endl << flush;
//  cerr << " ksurf: split_surf: RX_b = " << endl << RX_b << endl << flush;
//  cerr << " ksurf: split_surf: QY_b = " << endl << QY_b << endl << flush;
//  cerr << " ksurf: split_surf: RY_b = " << endl << RY_b << endl << flush;
//  cerr << " ksurf: split_surf: QZ_b = " << endl << QZ_b << endl << flush;
//  cerr << " ksurf: split_surf: RZ_b = " << endl << RZ_b << endl << flush;
//  cerr << " ksurf: split_surf: RZ_b.deg[0] = " << RZ_b.deg[0] << ", RZ_b.coeffs[0] = " << RZ_b.coeffs[0] << endl << flush;  
//  cerr << " ksurf: split_surf: -------------------- " << endl << flush;
  
  //  2. See whether or not the curve (X_b/W_b, Y_b/W_b, Z_b/W_b) is a line.
  //     If so then return the plane that intersects with *this at the line.
  
  if (QX_b.deg[0] == 1 && !RX_b.deg[0] && !sgn(RX_b.coeffs[0])
      &&
      QY_b.deg[0] == 1 && !RY_b.deg[0] && !sgn(RY_b.coeffs[0])
      &&
      QZ_b.deg[0] == 1 && !RZ_b.deg[0] && !sgn(RZ_b.coeffs[0])
      ||
      !X_b.deg[0] && Y_b.deg[0] == 1 && Z_b.deg[0] == 1 && !W_b.deg[0]
      ||
      X_b.deg[0] == 1 && !Y_b.deg[0] && Z_b.deg[0] == 1 && !W_b.deg[0]
      ||
      X_b.deg[0] == 1 && Y_b.deg[0] == 1 && !Z_b.deg[0] && !W_b.deg[0]
      ||
      !QX_b.deg[0] && !RX_b.deg[0] && !sgn(RX_b.coeffs[0])
      &&
      !QY_b.deg[0] && !RY_b.deg[0] && !sgn(RY_b.coeffs[0])
      &&
      sgn(QZ_b.coeffs[0])
      ||
      !QX_b.deg[0] && !RX_b.deg[0] && !sgn(RX_b.coeffs[0])
      &&
      sgn(QY_b.coeffs[0])
      &&
      !QZ_b.deg[0] && !RZ_b.deg[0] && !sgn(RZ_b.coeffs[0])
      ||
      sgn(QX_b.coeffs[0])
      &&
      !QY_b.deg[0] && !RY_b.deg[0] && !sgn(RY_b.coeffs[0])
      &&
      !QZ_b.deg[0] && !RZ_b.deg[0] && !sgn(RZ_b.coeffs[0]))
  {
//    cerr << " ksurf: split_surf: the curve is a line. " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  2-1. Choose 2 pts x and y on the curve (X_b/W_b, Y_b/W_b, Z_b/W_b).
    
    W_b_c = 0;
    
    while (!sgn(W_b_c))
    {
      c     = bigrational(1, cd++);
      W_b_c = W_b.evaluate(c);
    }
    
    x[0] = X_b.evaluate(c) / W_b_c;
    x[1] = Y_b.evaluate(c) / W_b_c;
    x[2] = Z_b.evaluate(c) / W_b_c;
    
    W_b_c = 0;
    
    while (!sgn(W_b_c))
    {
      c     = bigrational(1, cd++);
      W_b_c = W_b.evaluate(c);
    }
    
    y[0] = X_b.evaluate(c) / W_b_c;
    y[1] = Y_b.evaluate(c) / W_b_c;
    y[2] = Z_b.evaluate(c) / W_b_c;
    
//    cerr << " ksurf: split_surf: x = (" << x[0] << ", " << x[1] << ", " << x[2] << ") " << endl << flush;
//    cerr << " ksurf: split_surf: y = (" << y[0] << ", " << y[1] << ", " << y[2] << ") " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  2-2. Compute a normal vector n of *this at x, and let z = x + n.
    
    n[0] = Impl->derivative(0).evaluate(x);
    n[1] = Impl->derivative(1).evaluate(x);
    n[2] = Impl->derivative(2).evaluate(x);
    assert(sgn(n[0]) || sgn(n[1]) || sgn(n[2]));
    
//    cerr << " ksurf: split_surf: n = (" << n[0] << ", " << n[1] << ", " << n[2]  << ") " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    for (j = 0; j < 3; j++)
      z[j] = x[j] + n[j];
    
//    cerr << " ksurf: split_surf: z = (" << z[0] << ", " << z[1] << ", " << z[2] << ") " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  2-3. Compute a plane on which 3 pts x, y and z lie.
    
    assert(get_plane_coeffs(x, y, z, plane));
    
    d[0] = d[1] = d[2] = 1;
    I    = new K_RATPOLY(3, d);
    
    p[0] = 1;
    p[1] = 0;
    p[2] = 0;
    I->get_coeff(p) = plane[0];
    
    p[0] = 0;
    p[1] = 1;
    p[2] = 0;
    I->get_coeff(p) = plane[1];
    
    p[0] = 0;
    p[1] = 0;
    p[2] = 1;
    I->get_coeff(p) = plane[2];
    
    p[0] = 0;
    p[1] = 0;
    p[2] = 0;
    I->get_coeff(p) = plane[3];
    
    I->reduce_deg();
    
    found = 1;
  }
  
  //  3. See whether or not
  //       the curve (X_b/W_b, Y_b/W_b, Z_b/W_b) lies on some plane.
  
  if (!found)
  {
//    cerr << " ksurf: split_surf: the curve is on a plane. " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  3-1. Choose 3 pts x, y, and z on the curve (X_b/W_b, Y_b/W_b, Z_b/W_b).
    
    W_b_c = 0;
    
    while (!sgn(W_b_c))
    {
      c     = bigrational(1, cd++);
      W_b_c = W_b.evaluate(c);
    }
    
    x[0] = X_b.evaluate(c) / W_b_c;
    x[1] = Y_b.evaluate(c) / W_b_c;
    x[2] = Z_b.evaluate(c) / W_b_c;
    
    W_b_c = 0;
    
    while (!sgn(W_b_c))
    {
      c     = bigrational(1, cd++);
      W_b_c = W_b.evaluate(c);
    }
    
    y[0] = X_b.evaluate(c) / W_b_c;
    y[1] = Y_b.evaluate(c) / W_b_c;
    y[2] = Z_b.evaluate(c) / W_b_c;
    
    W_b_c = 0;
    
    while (!sgn(W_b_c))
    {
      c     = bigrational(1, cd++);
      W_b_c = W_b.evaluate(c);
    }
    
    z[0] = X_b.evaluate(c) / W_b_c;
    z[1] = Y_b.evaluate(c) / W_b_c;
    z[2] = Z_b.evaluate(c) / W_b_c;
    
//    cerr << " ksurf: split_surf: x = (" << x[0] << ", " << x[1] << ", " << x[2] << ") " << endl << flush;
//    cerr << " ksurf: split_surf: y = (" << y[0] << ", " << y[1] << ", " << y[2] << ") " << endl << flush;
//    cerr << " ksurf: split_surf: z = (" << z[0] << ", " << z[1] << ", " << z[2] << ") " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  3-2. See whether or nor
    //         there exists a plane on which 3 pts x, y, and z lie.
    
    if (get_plane_coeffs(x, y, z, plane))  
    {
      Plane =
        plane[0] * X_b + plane[1] * Y_b + plane[2] * Z_b + plane[3] * W_b;
//      cerr << " ksurf: split_surf: Plane = " << endl << Plane << endl << flush;
      
      if (!Plane.deg[0] && !sgn(Plane.coeffs[0]))
        //  if there exists some plane that intersects with *this at the curve
      {
        d[0] = d[1] = d[2] = 1;
        I    = new K_RATPOLY(3, d);
        
        p[0] = 1;
        p[1] = 0;
        p[2] = 0;
        I->get_coeff(p) = plane[0];
        
        p[0] = 0;
        p[1] = 1;
        p[2] = 0;
        I->get_coeff(p) = plane[1];
        
        p[0] = 0;
        p[1] = 0;
        p[2] = 1;
        I->get_coeff(p) = plane[2];
        
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
        I->get_coeff(p) = plane[3];
        
        I->reduce_deg();
        
        W_b_c = 0;
        
        while (!sgn(W_b_c))
        {
          c     = bigrational(1, cd++);
          W_b_c = W_b.evaluate(c);
        }
        
        x[0] = X_b.evaluate(c) / W_b_c;
        x[1] = Y_b.evaluate(c) / W_b_c;
        x[2] = Z_b.evaluate(c) / W_b_c;
        
        if (!sgn(I->evaluate(x)))
          found = 1;
        else  //  if *I does not vanish at x
          delete I;
      }
    }
  }
  
  t = 1;
  
  while (!found)
  {
//    cerr << " ksurf: split_surf: the curve may not be a line. " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  4. Compute a parametric representation (X2/W2, Y2/W2, Z2/W2) of
    //       the surface that intersects with *this
    //         at the curve (X_b/W_b, Y_b/W_b, Z_b/W_b).
    
    //  4-1. Choose a pt x on the curve (X_b/W_b, Y_b/W_b, Z_b/W_b)
    //         at random.
    
    W_b_c = 0;
    
    while (!sgn(W_b_c))
    {
      c     = bigrational(1, cd++);
      W_b_c = W_b.evaluate(c);
    }
    
    x[0] = X_b.evaluate(c) / W_b_c;
    x[1] = Y_b.evaluate(c) / W_b_c;
    x[2] = Z_b.evaluate(c) / W_b_c;
    
//    cerr << " ksurf: split_surf: x = (" << x[0] << ", " << x[1] << ", " << x[2] << ") " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    X2_b = X_b.add_var(i);
    Y2_b = Y_b.add_var(i);
    Z2_b = Z_b.add_var(i);
    W2_b = W_b.add_var(i);
    
    if (i == 0)
      T = K_RATPOLY(2, 0, 0);
    else  //  if (i == 1)
      T = K_RATPOLY(2, 1, 0);
    
//    cerr << " ksurf: split_surf: X2_b = " << endl << X2_b << endl << flush;
//    cerr << " ksurf: split_surf: Y2_b = " << endl << Y2_b << endl << flush;
//    cerr << " ksurf: split_surf: Z2_b = " << endl << Z2_b << endl << flush;
//    cerr << " ksurf: split_surf: W2_b = " << endl << W2_b << endl << flush;
//    cerr << " ksurf: split_surf: T = " << endl << T << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  4-2. Compute a normal vector n of *this at x (or a random vector).
    
    dI_dx = Impl->derivative(0);
    dI_dy = Impl->derivative(1);
    dI_dz = Impl->derivative(2);
    
//    cerr << " ksurf: split_surf: dI_dx = " << endl << dI_dx << endl << flush;
//    cerr << " ksurf: split_surf: dI_dy = " << endl << dI_dy << endl << flush;
//    cerr << " ksurf: split_surf: dI_dz = " << endl << dI_dz << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    if (dI_dx.deg[0] || dI_dx.deg[1] || dI_dx.deg[2] ||
        dI_dy.deg[0] || dI_dy.deg[1] || dI_dy.deg[2] ||
        dI_dz.deg[0] || dI_dz.deg[1] || dI_dz.deg[2])
    {
      n[0] = dI_dx.evaluate(x);
      n[1] = dI_dy.evaluate(x);
      n[2] = dI_dz.evaluate(x);
    }
    else
    {
      n[0] = rand();
      n[1] = rand();
      n[2] = rand();
    }
    
//    cerr << " ksurf: split_surf: n = (" << n[0] << ", " << n[1] << ", " << n[2]  << ") " << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    //  4-3. Compute a parametric representation (X2/W2, Y2/W2, Z2/W2) of
    //         the surface that intersects with *this at
    //         at the curve (X_b/W_b, Y_b/W_b, Z_b/W_b).
    
    X2 = X2_b + n[0] * T * W2_b;
    Y2 = Y2_b + n[1] * T * W2_b;
    Z2 = Z2_b + n[2] * T * W2_b;
    
//    cerr << " ksurf: split_surf: X2 = " << endl << X2 << endl << flush;
//    cerr << " ksurf: split_surf: Y2 = " << endl << Y2 << endl << flush;
//    cerr << " ksurf: split_surf: Z2 = " << endl << Z2 << endl << flush;
//    cerr << " ksurf: split_surf: X2.deg[0] = " << X2.deg[0] << ", X2.deg[1] = " << X2.deg[1] << endl << flush;
//    cerr << " ksurf: split_surf: Y2.deg[0] = " << Y2.deg[0] << ", Y2.deg[1] = " << Y2.deg[1] << endl << flush;
//    cerr << " ksurf: split_surf: Z2.deg[0] = " << Z2.deg[0] << ", Z2.deg[1] = " << Z2.deg[1] << endl << flush;
//    cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    
    while (!X2.deg[0] && !X2.deg[1] && !sgn(X2.coeffs[0])
           ||
           !Y2.deg[0] && !Y2.deg[1] && !sgn(Y2.coeffs[0])
           ||
           !Z2.deg[0] && !Z2.deg[1] && !sgn(Z2.coeffs[0]))
    {
      //  4-2. Let n be a random vector.
      
      n[0] = rand();
      n[1] = rand();
      n[2] = rand();
      
//      cerr << " ksurf: split_surf: n = (" << n[0] << ", " << n[1] << ", " << n[2]  << ") " << endl << flush;
//      cerr << " ksurf: split_surf: -------------------- " << endl << flush;
      
      X2 = X2_b + n[0] * T * W2_b;
      Y2 = Y2_b + n[1] * T * W2_b;
      Z2 = Z2_b + n[2] * T * W2_b;
      
//      cerr << " ksurf: split_surf: X2 = " << endl << X2 << endl << flush;
//      cerr << " ksurf: split_surf: Y2 = " << endl << Y2 << endl << flush;
//      cerr << " ksurf: split_surf: Z2 = " << endl << Z2 << endl << flush;
//      cerr << " ksurf: split_surf: X2.deg[0] = " << X2.deg[0] << ", X2.deg[1] = " << X2.deg[1] << endl << flush;
//      cerr << " ksurf: split_surf: Y2.deg[0] = " << Y2.deg[0] << ", Y2.deg[1] = " << Y2.deg[1] << endl << flush;
//      cerr << " ksurf: split_surf: Z2.deg[0] = " << Z2.deg[0] << ", Z2.deg[1] = " << Z2.deg[1] << endl << flush;
//      cerr << " ksurf: split_surf: -------------------- " << endl << flush;
    }
    
    W2 = W2_b;
    
    //  4-4. Compute *I by interpolation.
    
    found = implicitize(X2, Y2, Z2, W2, t++, I);
    
//    cerr << " ksurf: split_surf: found = " << found << endl << flush;
//    cerr << " ksurf: split_surf: ==================== " << endl << flush;
  }
  
//  cerr << " ksurf: split_surf: *I = " << endl << *I << endl << flush;
//  cerr << " ksurf: split_surf: -------------------- " << endl << flush;
  
  return K_SURF(*I);
}

//  int K_SURF :: Bezier_output(ostream& out_fs,
//                              const bigrational& l_s,
//                              const bigrational& h_s,
//                              const bigrational& l_t,
//                              const bigrational& h_t) const
//    compute Bezier surface representation of *this.
//    control points are ordered s.t.
//      if the thumb         of our right hand goes l_s -> h_s and
//         the forefinger                           l_t -> h_t then
//         the middle finger                   is the normal to the surface.

int K_SURF :: Bezier_output(ostream& out_fs,
                            const bigrational& l_s,
                            const bigrational& h_s,
                            const bigrational& l_t,
                            const bigrational& h_t) const
{
  //  Let
  //    s_sub = (h_s - l_s) * S + l_s & t_sub = (h_t - l_t) * T + l_t
  //  where S & T are newly introduced variables. Then, substitute
  //    s_sub for s leaving t, S, & T are variables,
  //  and then
  //    t_sub for t leaving S & T are variables.
  
  unsigned long i, j;
  long          d_s_sub[3];
  long          d_t_sub[2];
  long          max_deg_s, max_deg_t;
  long          p[2];
  
  d_s_sub[0] = d_t_sub[0] = 0;
  d_s_sub[1] = d_t_sub[1] = 1;
  d_s_sub[2] = 0;
  
  K_RATPOLY s_sub(3, d_s_sub);
  K_RATPOLY t_sub(2, d_t_sub);
  
  s_sub.coeffs[0] = h_s - l_s;
  s_sub.coeffs[1] = l_s;
  t_sub.coeffs[0] = h_t - l_t;
  t_sub.coeffs[1] = l_t;
  
  K_RATPOLY X_sub, Y_sub, Z_sub, W_sub;
  
  X_sub = *X;
  X_sub = X_sub.add_var(2);
  X_sub = X_sub.add_var(3);
  X_sub = X_sub.subst_expr(0, s_sub);
  X_sub = X_sub.subst_expr(0, t_sub);
  Y_sub = *Y;
  Y_sub = Y_sub.add_var(2);
  Y_sub = Y_sub.add_var(3);
  Y_sub = Y_sub.subst_expr(0, s_sub);
  Y_sub = Y_sub.subst_expr(0, t_sub);
  Z_sub = *Z;
  Z_sub = Z_sub.add_var(2);
  Z_sub = Z_sub.add_var(3);
  Z_sub = Z_sub.subst_expr(0, s_sub);
  Z_sub = Z_sub.subst_expr(0, t_sub);
  W_sub = *W;
  W_sub = W_sub.add_var(2);
  W_sub = W_sub.add_var(3);
  W_sub = W_sub.subst_expr(0, s_sub);
  W_sub = W_sub.subst_expr(0, t_sub);
  
  max_deg_s = X->deg[0];
  max_deg_t = X->deg[1];
  
  if (Y->deg[0] > max_deg_s)
    max_deg_s = Y->deg[0];
  
  if (Y->deg[1] > max_deg_t)
    max_deg_t = Y->deg[1];
  
  if (Z->deg[0] > max_deg_s)
    max_deg_s = Z->deg[0];
  
  if (Z->deg[1] > max_deg_t)
    max_deg_t = Z->deg[1];
  
  if (W->deg[0] > max_deg_s)
    max_deg_s = W->deg[0];
  
  if (W->deg[1] > max_deg_t)
    max_deg_t = W->deg[1];
  
  K_RATPOLY X_Bern, Y_Bern, Z_Bern, W_Bern;
  
  X_Bern = X_sub.conv_to_Bernstein(max_deg_s, max_deg_t);
  Y_Bern = Y_sub.conv_to_Bernstein(max_deg_s, max_deg_t);
  Z_Bern = Z_sub.conv_to_Bernstein(max_deg_s, max_deg_t);
  W_Bern = W_sub.conv_to_Bernstein(max_deg_s, max_deg_t);
  
  //  Dump max. degrees in s and t.
  
  out_fs << max_deg_s << "  " << max_deg_t << endl << flush;
  
  //  Dump control points.
  
  for (i = 0; i <= max_deg_s; i++)
  {
    p[0] = i;
    
    for (j = 0; j <= max_deg_t; j++)
    {
      p[1] = j;
      
      out_fs <<
        X_Bern.get_coeff(p).as_double() << "  " <<
        Y_Bern.get_coeff(p).as_double() << "  " <<
        Z_Bern.get_coeff(p).as_double() << "  " <<
        W_Bern.get_coeff(p).as_double() << "  " <<
        endl << flush;
    }
  }
  
  return 0;
}

unsigned long get_all_int_pts(const K_RATPOLY& P1,
                              const K_SURF& S2, const K_SURF& S3,
                              K_POINT2D**& p)
{
  assert(P1.get_num_vars() == 2);
  
  unsigned long i;
  K_RATPOLY     P2;
  bigrational   l_s, h_s, l_t, h_t;
  unsigned long num_int_pts;
  
//  cerr << " ksurf: get_all_int_pts: ---------- " << endl << flush;
//  cerr << " S2.Impl = " << endl << *S2.Impl << endl << flush;
//  cerr << " S3.Impl = " << endl << *S3.Impl << endl << flush;
//  cerr << " S3.X = " << endl << *S3.X << endl << flush;
//  cerr << " S3.Y = " << endl << *S3.Y << endl << flush;
//  cerr << " S3.Z = " << endl << *S3.Z << endl << flush;
//  cerr << " S3.W = " << endl << *S3.W << endl << flush;
  
  P2 = S2.Impl->subst_param_expr(*S3.X, *S3.Y, *S3.Z, *S3.W);
  
//  cerr << " P1 = " << endl << P1 << endl << flush;
//  cerr << " P2 = " << endl << P2 << endl << flush;
  
  h_s = GoodSylvester(P1, P2, 1).get_Mignotte_bd() + s_epsilon;
  l_s = - h_s - s_epsilon;
  h_t = GoodSylvester(P1, P2, 0).get_Mignotte_bd() + t_epsilon;
  l_t = - h_t - t_epsilon;
  
  num_int_pts = get_pts(l_s, h_s, l_t, h_t, P1, P2, p, 0, 0);
//  for (i = 0; i < num_int_pts; i++)
//    cerr << " p[" << i << "] = " << *p[i] << endl << flush;
//  cerr << " ksurf: get_all_int_pts: ---------- " << endl << flush;
  
//  cerr << " XYZ: S3.X->deg[0] = " << S3.X->deg[0] << ", S3.X->deg[1] = " << S3.X->deg[1] << endl << flush;
//  cerr << " XYZ: S3.Y->deg[0] = " << S3.Y->deg[0] << ", S3.Y->deg[1] = " << S3.Y->deg[1] << endl << flush;
//  cerr << " XYZ: S3.Z->deg[0] = " << S3.Z->deg[0] << ", S3.Z->deg[1] = " << S3.Z->deg[1] << endl << flush;
//  cerr << " XYZ: S3.W->deg[0] = " << S3.W->deg[0] << ", S3.W->deg[1] = " << S3.W->deg[1] << endl << flush;
//  cerr << " ksurf: get_all_int_pts: ---------- " << endl << flush;
  
  return num_int_pts;
}

