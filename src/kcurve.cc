//  file:    kcurve.cc
//  update:  11/13/02

#include <config.h>

#include <kcurve.h>

K_CURVE :: K_CURVE()
{
  poly = 0;
  
  segments     = 0;
  num_segments = 0;
  
  curve_in_other_dom = 0;
  dir_in_other_dom   = 0;
  
  ref_count = 0;
}

K_CURVE :: K_CURVE(const K_RATPOLY& P,
                   K_SEGMENT* const s[], const unsigned long n)
{
  unsigned long i;
  
  poly = new K_RATPOLY(P);
  poly->ref_count++;
  
  if ((num_segments = n) > 0)
  {
    segments = new K_SEGMENT* [num_segments];
    
    for (i = 0; i < num_segments; i++)
    {
      segments[i] = s[i];
      segments[i]->ref_count++;
    }
  }
  else  //  if (num_segments == 0)
    segments = 0;
  
  curve_in_other_dom = 0;
  dir_in_other_dom   = 0;
  
  ref_count = 0;
}

K_CURVE :: K_CURVE(K_RATPOLY* const P,
                   K_SEGMENT* const s[], const unsigned long n)
{
  unsigned long i;
  
  if (poly = P)
    poly->ref_count++;
  
  if ((num_segments = n) > 0)
  {
    segments = new K_SEGMENT* [num_segments];
    
    for (i = 0; i < num_segments; i++)
    {
      segments[i] = s[i];
      segments[i]->ref_count++;
    }
  }
  else  //  if (num_segments == 0)
    segments = 0;
  
  curve_in_other_dom = 0;
  dir_in_other_dom   = 0;
  
  ref_count = 0;
}

K_CURVE :: K_CURVE(const K_CURVE& c,
                   const unsigned long s, const unsigned long e)
{
  unsigned long i;
  
  if (poly = c.poly)
    poly->ref_count++;
  
  if (s <= e)
  {
    segments = new K_SEGMENT* [num_segments = e - s + 1];
    //  s <= e => num_segments >= 1
    
    for (i = s; i <= e; i++)
    {
      segments[i - s] = c.segments[i];
      segments[i - s]->ref_count++;
    }
  }
  else  //  if (s > e)
  {
    assert(c.is_closed());
    
    segments = new K_SEGMENT* [num_segments = c.num_segments - s + e + 1];
    //  num_segments >= 1
    
    for (i = s; i < c.num_segments; i++)
    {
      segments[i - s] = c.segments[i];
      segments[i - s]->ref_count++;
    }
    
    for (i = 0; i <= e; i++)
    {
      segments[i + c.num_segments - s] = c.segments[i];
      segments[i + c.num_segments - s]->ref_count++;
    }
  }
  
  curve_in_other_dom = 0;
  dir_in_other_dom   = 0;
  
  ref_count = 0;
}

K_CURVE :: K_CURVE(const K_CURVE& c)
{
  unsigned long i;
  
  if (poly = c.poly)
    poly->ref_count++;
  
  if ((num_segments = c.num_segments) > 0)
  {
    segments = new K_SEGMENT* [num_segments];
    
    for (i = 0; i < num_segments; i++)
    {
      segments[i] = c.segments[i];
      segments[i]->ref_count++;
    }
  }
  else  //  if (num_segments == 0)
    segments = 0;
  
//  curve_in_other_dom = c.curve_in_other_dom;
//  dir_in_other_dom   = c.dir_in_other_dom;
  curve_in_other_dom = 0;
  dir_in_other_dom   = 0;
  
  ref_count = 0;
}

K_CURVE& K_CURVE :: operator =(const K_CURVE& c)
{
  unsigned long i;
  
  if (this != &c)
  {
    if (poly && !--poly->ref_count)
      delete poly;
    
    if (num_segments > 0)
    {
      for (i = 0; i < num_segments; i++)
        if (!--segments[i]->ref_count)
          delete segments[i];
      
      delete [] segments;
    }
    
    if (curve_in_other_dom)
    {
      curve_in_other_dom->curve_in_other_dom = 0;
      curve_in_other_dom->dir_in_other_dom   = 0;
    }
    
    if (poly = c.poly)
      poly->ref_count++;
    
    if ((num_segments = c.num_segments) > 0)
    {
      unsigned long i;
      
      segments = new K_SEGMENT* [num_segments];
      
      for (i = 0; i < num_segments; i++)
      {
        segments[i] = c.segments[i];
        segments[i]->ref_count++;
      }
    }
    else  //  if (num_segments == 0)
      segments = 0;
    
//    curve_in_other_dom = c.curve_in_other_dom;
//    dir_in_other_dom   = c.dir_in_other_dom;
    curve_in_other_dom = 0;
    dir_in_other_dom   = 0;
  }
  
  return *this;
}

K_CURVE :: ~K_CURVE()
{
  unsigned long i;
  
  if (poly && !--poly->ref_count)
    delete poly;
  
  if (num_segments > 0)
  {
    for (i = 0; i < num_segments; i++)
      if (!--segments[i]->ref_count)
        delete segments[i];
    
    delete [] segments;
  }
//  
//  if (curve_in_other_dom)
//  {
//    curve_in_other_dom->curve_in_other_dom = 0;
//    curve_in_other_dom->dir_in_other_dom   = 0;
//  }
}

ostream& K_CURVE :: output(ostream& o) const
{
  if (num_segments > 0)
  {
    unsigned long i;
    
    o << "<";
    
    for (i = 0; i < num_segments - 1; i++)
//      o << *segments[i] << ", ";
      o << *segments[i] << endl;
    
    o << *segments[num_segments - 1] << ">" << flush;
  }
  else  //  if (num_segments == 0)
    o << " NULL " << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const K_CURVE& x)
{
  return x.output(o);
}

int K_CURVE :: assoc(K_CURVE* const c, const int dir)
{
  assert(c);
  
  curve_in_other_dom    = c;
  c->curve_in_other_dom = this;
  dir_in_other_dom      = c->dir_in_other_dom = dir;
}

K_POINT2D* K_CURVE :: start() const
{
  assert(num_segments > 0);
  
  return segments[0]->start;
}

K_POINT2D* K_CURVE :: end() const
{
  assert(num_segments > 0);
  
  return segments[num_segments - 1]->end;
}

K_BOXCO2 K_CURVE :: bbox() const
{
  assert(num_segments > 0);
  
  unsigned long i;
  K_BOXCO2      b;
  
  b = segments[0]->outer_box();
  
  for (i = 1; i < num_segments; i++)
    b = b.merge(segments[i]->outer_box());
  
  return b;
}

int K_CURVE :: is_closed() const
{
  return start()->equal(*end());
}

int K_CURVE :: rotate_closed_curve(const unsigned long i)
{
  assert(i < num_segments);  //  i >= 0 => num_segments > 0
  assert(is_closed());
  
  unsigned long j;
  K_SEGMENT**   s;
  
  s = new K_SEGMENT* [num_segments];
  
  for (j = i; j < num_segments; j++)
  {
    s[j - i] = segments[j];
    s[j - i]->ref_count++;
  }
  
  for (j = 0; j < i; j++)
  {
    s[j + num_segments - i] = segments[j];
    s[j + num_segments - i]->ref_count++;
  }
  
  for (j = 0; j < num_segments; j++)
    if (!--segments[j]->ref_count)
      delete segments[j];
  
  delete [] segments;  //  num_segments > 0 => segments != 0
  
  segments = s;
  
  return 0;
}

int K_CURVE :: sort_pts(K_POINT2D** const pts,
                        const unsigned long num_pts) const
{
  long           i, j, k;
  int            c, s0, s1;
  K_POINT2D***   pts_on_seg;
  unsigned long* num_pts_on_seg;
  
  pts_on_seg     = new K_POINT2D** [num_segments];
  num_pts_on_seg = new unsigned long [num_segments];
  
  for (i = 0; i < num_segments; i++)
  {
    pts_on_seg[i]     = new K_POINT2D* [num_pts + 2];
    pts_on_seg[i][0]  = segments[i]->start;
    pts_on_seg[i][0]->ref_count++;
    num_pts_on_seg[i] = 1;
  }
  
  for (i = 0; i < num_pts; i++)
  {
    for (c = 0, j = 0; !c && j < num_segments; j++)
      if (c = segments[j]->contains(*pts[i]))
      {
        pts_on_seg[j][num_pts_on_seg[j]] = pts[i];
        pts_on_seg[j][num_pts_on_seg[j]]->ref_count++;
        num_pts_on_seg[j]++;
      }
    
    assert(c);
  }
  
  for (i = 0; i < num_segments; i++)
  {
    if (!sort_s(pts_on_seg[i], num_pts_on_seg[i])
        &&
        !sort_t(pts_on_seg[i], num_pts_on_seg[i]))
    {
      for (j = 0; j < num_pts_on_seg[i] - 1; j++)
        for (k = j + 1; k < num_pts_on_seg[i]; k++)
          if (!pts_on_seg[i][j]->equal(*pts_on_seg[i][k]))
            pts_on_seg[i][j]->separate(*pts_on_seg[i][k]);
      
      assert(sort_s(pts_on_seg[i], num_pts_on_seg[i])
             ||
             sort_t(pts_on_seg[i], num_pts_on_seg[i]));
    }
  }
  
  for (i = k = 0; i < num_segments; i++)
    if (pts_on_seg[i][0]->equiv(*segments[i]->start))
      for (j = 1; j < num_pts_on_seg[i]; j++)
        pts[k++] = pts_on_seg[i][j];
    else
      for (j = num_pts_on_seg[i] - 2; j >= 0; j--)
        pts[k++] = pts_on_seg[i][j];
  
  if (num_segments > 0)
  {
    for (i = 0; i < num_segments; i++)
    {
      for (j = 0; j < num_pts_on_seg[i]; j++)
        if (!--pts_on_seg[i][j]->ref_count)
          delete pts_on_seg[i][j];
      
      delete [] pts_on_seg[i];
    }
    
    delete [] pts_on_seg;
    delete [] num_pts_on_seg;
  }
  
  return 0;
}

unsigned long K_CURVE :: locate_pt_seg_start(K_POINT2D& x)
{
  unsigned long i;
  
  for (i = 0; i < num_segments && !segments[i]->start->equal(x); i++)
    ;
  
  return i;
}

unsigned long K_CURVE :: locate_pt_seg_end(K_POINT2D& x)
{
  unsigned long i;
  
  for (i = 0; i < num_segments && !segments[i]->end->equal(x); i++)
    ;
  
  return i;
}

int K_CURVE :: is_start_or_end(K_POINT2D& x)
{
  return start()->equal(x) || end()->equal(x);
}

//  int K_CURVE :: contains(K_POINT2D& x, const int count_end)
//    tells whether or not x lies on *this, i.e.,
//      there exists segments[i] s.t. segments[i]->contains(x) == 1.
//    x lying at start is considered to be OUT.
//    x lying at end is consider to be IN  if count_at_end == 1, and
//                                     OUT otherwise.

int K_CURVE :: contains(K_POINT2D& x, const int count_end)
{
  unsigned long i;
  int           c;
  
  if (!num_segments)
    c = 0;
  else if (!is_closed() && start()->equal(x))
  //  Unless *this is a loop, start is OUT of *this.
    c = 0;
  else if (end()->equal(x))
  //  end is IN *this if count_at_end > 0 and OUT otherwise.
    c = count_end;
  else
  {
    for (i = 0; i < num_segments && !segments[i]->contains(x); i++)
      ;
    
    if (i < num_segments)
    {
//      cerr << " kcurve: contains: *segments[" << i << "] = " << *segments[i] << endl << flush;
//      cerr << " kcurve: contains: x = " << x << endl << flush;
      c = 1;
    }
    else  //  if (i == num_segments)
      c = 0;
  }
  
  return c;
}

//  K_POINT2D K_CURVE :: pt_on()
//    returns a pt on *this that is neither start nor end.

K_POINT2D K_CURVE :: pt_on()
{
  assert(num_segments > 0);
  
  unsigned long i, j;
  int           go_on;
  K_BOXCO2      b;
  K_RATPOLY     L;
  K_POINT2D**   int_pts;
  unsigned long num_int_pts;
  K_POINT2D     x;
  
  //  If there is some segment[i]->end that is NEITHER start NOR end then
  //    return it.
  //  Otherwise
  //    return an intersection of *this and the line that bisects bbox().
  
  //  See if some segment[i]->end that is NEITHER start NOR end.
  
  for (i = 0;
       i < num_segments - 1
       &&
       (segments[i]->end->overlap(*start())
        ||
        segments[i]->end->overlap(*end()));
       i++)
    ;
  
  if (i < num_segments - 1)
    x = K_POINT2D(*segments[i]->end);
  else  //  if (i == num_segments - 1)
  {
    go_on = 0;
    
    while (!go_on)
    {
      //  Compute intersections of *this and the line that bisects bbox().
      
      b = bbox();
      
      if (b.high[0] - b.low[0] > b.high[1] - b.low[1])
        L = K_RATPOLY(2, 0, (b.low[0] + b.high[0]) / 2);
      else  //  if (b.high[0] - b.low[0] <= b.high[1] - b.low[1])
        L = K_RATPOLY(2, 1, (b.low[1] + b.high[1]) / 2);
      
      if ((num_int_pts = find_intersections(L, int_pts, 1)) > 0)
      {
        //  See if some int_pts[i] that is NEITHER start NOR end.
        
        for (i = 0;
             i < num_int_pts
             &&
             (int_pts[i]->overlap(*start()) || int_pts[i]->overlap(*end()));
             i++)
          ;
        
        if (i < num_int_pts)
        {
          x     = K_POINT2D(*int_pts[i]);
          go_on = 1;
        }
        
        for (i = 0; i < num_int_pts; i++)
          if (!--int_pts[i]->ref_count)
            delete int_pts[i];
        
        delete [] int_pts;
      }
      
      if (!go_on)
      {
        //  start or end or segments[i]->end is too fat. Shrink them.
        
        segments[0]->start->shrink(shrink_step, shrink_step);
        
        for (i = 0; i < num_segments; i++)
          segments[i]->end->shrink(shrink_step, shrink_step);
      }
    }
  }
  
  return x;
}

//  int K_CURVE :: add_pt(K_POINT2D* x)
//    Add x to *this as an endpoint of some segments.

int K_CURVE :: add_pt(K_POINT2D* const x, const int added_at_start)
{
  unsigned long i, j;
  K_SEGMENT**   s;
  int           a;
  
  //  Find the segment[i], if any, on which x lies.
  
  for (i = 0; i < num_segments && !segments[i]->contains(*x); i++)
    ;
  
  if (i < num_segments) //  x \in ( segment[i]->start, segment[i]->end ]
  {
    if (!segments[i]->end->equiv(*x))
    //  Both segments[i]->end & x have already been cached
    //  by segments[i]->contains(*x).
    //  x \in ( segment[i]->start, segment[i]->end )
    //  Chop segment[i] at x into 2 pieces.
    {
      s = new K_SEGMENT* [num_segments + 1];
      
      for (j = 0; j < i; j++)
      {
        s[j] = segments[j];
        s[j]->ref_count++;
      }
      
      s[i]     = new K_SEGMENT(segments[i]->start, x);
      s[i]->ref_count++;
      s[i + 1] = new K_SEGMENT(x, segments[i]->end);
      s[i + 1]->ref_count++;
      
      for (j = i + 1; j < num_segments; j++)
      {
        s[j + 1] = segments[j];
        s[j + 1]->ref_count++;
      }
      
      for (j = 0; j < num_segments; j++)
        if (!--segments[j]->ref_count)
          delete segments[j];
      
      delete [] segments;  //  0 <= i < num_segments => segments != 0
      
      segments = s;
      num_segments++;
    }
    
    a = 1;
  }
  else  //  if (i == num_segments)
    if (segments[0]->start->equiv(*x))
    //  Both segments[0]->start & x have already been cached
    //  by segments[i]->contains(*x).
      a = added_at_start;
    else
      a = 0;
  
  return a;
}

//  unsigned long K_CURVE :: find_intersections(const K_RATPOLY& P,
//                                              K_POINT2D**& intersections,
//                                              const int count_end)
//    returns the length of the array of the intersection of *this and P.

unsigned long K_CURVE :: find_intersections(const K_RATPOLY& P,
                                            K_POINT2D**& intersections,
                                            const int count_at_end)
{
  assert(num_segments > 0);
  
  unsigned long i, j;
  bigrational   l_s, h_s, l_t, h_t, v_s, v_t;
  K_POINT2D**   pts;
  unsigned long num_pts;
  K_POINT2D**   intersections_proto;
  unsigned long num_intersections;
  
  //  Compute the ranges the polynomials of the curves.
  
  l_s = start()->get_low_s();
  h_s = start()->get_high_s();
  l_t = start()->get_low_t();
  h_t = start()->get_high_t();
  
  for (i = 1; i < num_segments; i++)
  {
    if ((v_s = segments[i]->start->get_low_s()) < l_s)
      l_s = v_s;
    
    if ((v_s = segments[i]->start->get_high_s()) > h_s)
      h_s = v_s;
    
    if ((v_t = segments[i]->start->get_low_t()) < l_t)
      l_t = v_t;
    
    if ((v_t = segments[i]->start->get_high_t()) > h_t)
      h_t = v_t;
  }
  
  if (!is_closed())
  {
    if ((v_s = end()->get_low_s()) < l_s)
      l_s = v_s;
    
    if ((v_s = end()->get_high_s()) > h_s)
      h_s = v_s;
    
    if ((v_t = end()->get_low_t()) < l_t)
      l_t = v_t;
    
    if ((v_t = end()->get_high_t()) > h_t)
      h_t = v_t;
  }
  
  //  Compute all the intersections of poly & P in the range.
  
  if (l_s == h_s && l_t == h_t)
  {
    bigrational v[2];
    
    v[0] = l_s;
    v[1] = l_t;
    
    if (!sgn(poly->evaluate(v)) && !sgn(P.evaluate(v)))
    {
      pts    = new K_POINT2D* [num_pts = 1];
      pts[0] = new K_POINT2D(l_s, l_t, *poly, P);
      pts[0]->ref_count++;
    }
    else
    {
      pts     = 0;
      num_pts = 0;
    }
  }
  else if (l_s == h_s && l_t != h_t)
    num_pts = get_pts_proto(l_s,
                            l_t, h_t, P.subst_val(0, l_s),
                            *poly, P,
                            pts,
                            init_tol, 1);
  else if (l_s != h_s && l_t == h_t)
    num_pts = get_pts_proto(l_s, h_s, P.subst_val(1, l_t),
                            l_t,
                            *poly, P,
                            pts,
                            init_tol, 1);
  else  //  if (l_s != h_s && l_t != h_t)
    num_pts = get_pts(l_s, h_s, l_t, h_t, *poly, P, pts, init_tol, 1);
  
  //  Pick up ones lying on *this.
  
  intersections_proto = new K_POINT2D* [num_pts];
  num_intersections   = 0;
  
  for (i = 0; i < num_pts; i++)
    if (contains(*pts[i], count_at_end))
    {
      intersections_proto[num_intersections] = pts[i];
      intersections_proto[num_intersections++]->ref_count++;
      
      if (!--pts[i]->ref_count)
        delete pts[i];
    }
  
  if (num_pts > 0)
    delete [] pts;
  
  intersections = new K_POINT2D* [num_intersections];
  
  for (i = 0; i < num_intersections; i++)
  {
    intersections[i] = intersections_proto[i];
    intersections[i]->ref_count++;
    
    if (!--intersections_proto[i]->ref_count)
      delete intersections_proto[i];
  }
  
  if (num_pts > 0)
    delete [] intersections_proto;
  
  return num_intersections;
}

//  int K_CURVE :: add_on(const K_CURVE& c)
//    returns *this added to c.

int K_CURVE :: add_on(const K_CURVE& c)
{
  assert(*poly == *c.poly);
  assert(num_segments > 0);
  assert(c.num_segments > 0);
  
  unsigned long i, j;
  K_SEGMENT**   s;
  
  s = new K_SEGMENT* [num_segments + c.num_segments];
  
  if (start()->equal(*c.start()))
  {
    for (i = 0; i < num_segments; i++)
    {
      s[i] = new K_SEGMENT(segments[num_segments - 1 - i]->reverse());
      s[i]->ref_count++;
    }
    
    for (i = num_segments, j = 0; j < c.num_segments; i++, j++)
    {
      s[i] = c.segments[j];
      s[i]->ref_count++;
    }
  }
  else if (start()->equal(*c.end()))
  {
    for (i = 0; i < c.num_segments; i++)
    {
      s[i] = c.segments[i];
      s[i]->ref_count++;
    }
    
    for (i = c.num_segments, j = 0; j < num_segments; i++, j++)
    {
      s[i] = segments[j];
      s[i]->ref_count++;
    }
  }
  else if (end()->equal(*c.start()))
  {
    for (i = 0; i < num_segments; i++)
    {
      s[i] = segments[i];
      s[i]->ref_count++;
    }
    
    for (i = num_segments, j = 0; j < c.num_segments; i++, j++)
    {
      s[i] = c.segments[j];
      s[i]->ref_count++;
    }
  }
  else if (end()->equiv(*c.end()))
  //  Both end & c.end have already been cached.
  {
    for (i = 0; i < num_segments; i++)
    {
      s[i] = segments[i];
      s[i]->ref_count++;
    }
    
    for (i = num_segments, j = 0; j < c.num_segments; i++, j++)
    {
      s[i] = new K_SEGMENT(c.segments[c.num_segments - 1 - j]->reverse());
      s[i]->ref_count++;
    }
  }
  else
    assert(1);
  
  for (i = 0; i < num_segments; i++)
    if (!--segments[i]->ref_count)
      delete segments[i];
  
  delete [] segments;  //  num_segments > 0 => segments != 0
  
  num_segments += c.num_segments;
  segments      = s;
  
  curve_in_other_dom = 0;
  dir_in_other_dom   = 0;
  
  return 0;
}

int K_CURVE :: split(const unsigned long s, K_CURVE& c1, K_CURVE&c2)
{
  assert(s < num_segments - 1);  //  num_segments - 1 > s >= 0
  
  unsigned long i, j;
  
  K_SEGMENT** s1;
  K_SEGMENT** s2;
  
  s1 = new K_SEGMENT* [s + 1];
  s2 = new K_SEGMENT* [num_segments - s - 1];
  
  for (i = 0; i < s + 1; i++)
  {
    s1[i] = segments[i];
    s1[i]->ref_count++;
  }
  
  for (i = s + 1; i < num_segments; i++)
  {
    s2[i - s - 1] = segments[i];
    s2[i - s - 1]->ref_count++;
  }
  
  c1 = K_CURVE(poly, s1, s + 1);
  c2 = K_CURVE(poly, s2, num_segments - s - 1);
  
  for (i = 0; i < s + 1; i++)
    if (!--s1[i]->ref_count)
      delete s1[i];
  
  delete [] s1;
  
  for (i = 0; i < num_segments - s - 1; i++)
    if (!--s2[i]->ref_count)
      delete s2[i];
  
  delete [] s2;
  
  return 0;
}

int K_CURVE :: reverse()
{
  unsigned long i;
  K_SEGMENT**   s;

  if (num_segments > 0)
  {
    s = new K_SEGMENT* [num_segments];
    
    for (i = 0; i < num_segments; i++)
    {
      s[i] = new K_SEGMENT(segments[num_segments - 1 - i]->reverse());
      s[i]->ref_count++;
    }
    
    for (i = 0; i < num_segments; i++)
      if (!--segments[i]->ref_count)
        delete segments[i];
    
    delete [] segments;
    
    segments = s;
  }
  
  if (curve_in_other_dom)
    curve_in_other_dom->dir_in_other_dom *= - 1;
  
  dir_in_other_dom *= - 1;
  
  return 0;
}

//  unsigned long gen_curve_topo(const K_RATPOLY& P,
//                               const bigrational& l_s,
//                               const bigrational& h_s,
//                               const bigrational& l_t,
//                               const bigrational& h_t,
//                               K_CURVE**& curves)
//    resolves topology of a bivariate polynomial P
//                      in a region [l_s, h_s] x [l_t, h_t]
//    to obtain monotone curves
//      curves[0], curves[1], ..., curves[num_curves - 1].
//    returns num_curves.

unsigned long gen_curve_topo(const K_RATPOLY& P,
                             const bigrational& l_s, const bigrational& h_s,
                             const bigrational& l_t, const bigrational& h_t,
                             K_CURVE**& curves)
{
  assert(P.num_vars == 2);
  assert(l_s < h_s);
  assert(l_t < h_t);
  
  unsigned long i, j;
  K_RATPOLY     P_l_s, P_h_s, P_l_t, P_h_t, L_s, H_s, L_t, H_t;
  K_POINT2D**   edge_l_s;
  K_POINT2D**   edge_h_s;
  K_POINT2D**   edge_l_t;
  K_POINT2D**   edge_h_t;
  unsigned long num_edge_l_s, num_edge_h_s, num_edge_l_t, num_edge_h_t;
  K_RATPOLY     dP_ds, dP_dt;
  K_POINT2D**   turn_s;
  K_POINT2D**   turn_t;
  unsigned long num_turn_s, num_turn_t;
  unsigned long num_curves;
  
//  cerr << " kcurve: gen_curve_topo: -------------------- " << endl << flush;
//  cerr << "   P = " << endl << P << endl << flush;
//  cerr << "   [" << l_s << ", " << h_s << "] x [" << l_t << ", " << h_t << "]" << endl << flush;
//  cerr << " -------------------------------------------- " << endl << flush;
  
  P_l_s = P.subst_val(0, l_s);
  P_h_s = P.subst_val(0, h_s);
  P_l_t = P.subst_val(1, l_t);
  P_h_t = P.subst_val(1, h_t);
  
  L_s = K_RATPOLY(1, 0, l_s);
  H_s = K_RATPOLY(1, 0, h_s);
  L_t = K_RATPOLY(1, 0, l_t);
  H_t = K_RATPOLY(1, 0, h_t);
  
//  cerr << "   P_l_s = " << endl << P_l_s << endl << flush;
//  cerr << "   P_h_s = " << endl << P_h_s << endl << flush;
//  cerr << "   P_l_t = " << endl << P_l_t << endl << flush;
//  cerr << "   P_h_t = " << endl << P_h_t << endl << flush;
//  cerr << " -------------------------------------------- " << endl << flush;
  
  num_edge_l_s =
    get_pts_proto(l_s, l_t, h_t, P_l_s, P, L_s, edge_l_s, init_tol, 1);
  num_edge_h_s =
    get_pts_proto(h_s, l_t, h_t, P_h_s, P, H_s, edge_h_s, init_tol, 1);
  num_edge_l_t =
    get_pts_proto(l_s, h_s, P_l_t, l_t, P, L_t, edge_l_t, init_tol, 0);
  num_edge_h_t =
    get_pts_proto(l_s, h_s, P_h_t, h_t, P, H_t, edge_h_t, init_tol, 0);
  
  if (num_edge_l_s > 1)
    sort_t(edge_l_s, num_edge_l_s);
  
  if (num_edge_h_s > 1)
    sort_t(edge_h_s, num_edge_h_s);
  
  if (num_edge_l_t > 1)
    sort_s(edge_l_t, num_edge_l_t);
  
  if (num_edge_h_t > 1)
    sort_s(edge_h_t, num_edge_h_t);
  
//  for (i = 0; i < num_edge_l_s; i++)
//    cerr << "   edge_l_s[" << i << "] = " << *edge_l_s[i] << endl << flush;
//  for (i = 0; i < num_edge_h_s; i++)
//    cerr << "   edge_h_s[" << i << "] = " << *edge_h_s[i] << endl << flush;
//  for (i = 0; i < num_edge_l_t; i++)
//    cerr << "   edge_l_t[" << i << "] = " << *edge_l_t[i] << endl << flush;
//  for (i = 0; i < num_edge_h_t; i++)
//    cerr << "   edge_h_t[" << i << "] = " << *edge_h_t[i] << endl << flush;
//  cerr << " -------------------------------------------- " << endl << flush;
  
  dP_ds = P.derivative(0);
  dP_dt = P.derivative(1);
  
//  cerr << "   dP_ds = " << endl << dP_ds << endl << flush;
//  cerr << "   dP_dt = " << endl << dP_dt << endl << flush;
//  cerr << " -------------------------------------------- " << endl << flush;
  
  num_turn_s = get_pts(l_s, h_s, l_t, h_t, P, dP_ds, turn_s, init_tol, 0);
  num_turn_t = get_pts(l_s, h_s, l_t, h_t, P, dP_dt, turn_t, init_tol, 0);
  
//  for (i = 0; i < num_turn_s; i++)
//    cerr << "   turn_s[" << i << "] = " << *turn_s[i] << endl << flush;
//  for (i = 0; i < num_turn_t; i++)
//    cerr << "   turn_t[" << i << "] = " << *turn_t[i] << endl << flush;
//  cerr << " -------------------------------------------- " << endl << flush;
  
  num_curves = gen_curve_topo_proto(P, l_s, h_s, l_t, h_t,
                                    edge_l_s, num_edge_l_s,
                                    edge_h_s, num_edge_h_s,
                                    edge_l_t, num_edge_l_t,
                                    edge_h_t, num_edge_h_t,
                                    turn_s, num_turn_s,
                                    turn_t, num_turn_t,
                                    curves);
  
  if (num_edge_l_s > 0)
  {
    for (i = 0; i < num_edge_l_s; i++)
      if (!--edge_l_s[i]->ref_count)
        delete edge_l_s[i];
    
    delete [] edge_l_s;
  }
  
  if (num_edge_h_s > 0)
  {
    for (i = 0; i < num_edge_h_s; i++)
      if (!--edge_h_s[i]->ref_count)
        delete edge_h_s[i];
    
    delete [] edge_h_s;
  }
  
  if (num_edge_l_t > 0)
  {
    for (i = 0; i < num_edge_l_t; i++)
      if (!--edge_l_t[i]->ref_count)
        delete edge_l_t[i];
    
    delete [] edge_l_t;
  }
  
  if (num_edge_h_t > 0)
  {
    for (i = 0; i < num_edge_h_t; i++)
      if (!--edge_h_t[i]->ref_count)
        delete edge_h_t[i];
    
    delete [] edge_h_t;
  }
  
  if (num_turn_s > 0)
  {
    for (i = 0; i < num_turn_s; i++)
      if (!--turn_s[i]->ref_count)
        delete turn_s[i];
    
    delete [] turn_s;
  }
  
  if (num_turn_t > 0)
  {
    for (i = 0; i < num_turn_t; i++)
      if (!--turn_t[i]->ref_count)
        delete turn_t[i];
    
    delete [] turn_t;
  }
  
//  cerr << "   num_curves = " << num_curves << endl << flush;
//  for (i = 0; i < num_curves; i++)
//    cerr << "   curves[" << i << "]->num_segments = " << curves[i]->num_segments << endl << flush;
//  cerr << " kcurve: gen_curve_topo: -------------------- " << endl << flush;
  
  return num_curves;
}

unsigned long gen_curve_topo_proto(const K_RATPOLY& P,
                                   const bigrational& l_s,
                                   const bigrational& h_s,
                                   const bigrational& l_t,
                                   const bigrational& h_t,
                                   K_POINT2D** const edge_l_s,
                                   const unsigned long num_edge_l_s,
                                   K_POINT2D** const edge_h_s,
                                   const unsigned long num_edge_h_s,
                                   K_POINT2D** const edge_l_t,
                                   const unsigned long num_edge_l_t,
                                   K_POINT2D** const edge_h_t,
                                   const unsigned long num_edge_h_t,
                                   K_POINT2D** const turn_s,
                                   const unsigned long num_turn_s,
                                   K_POINT2D** const turn_t,
                                   const unsigned long num_turn_t,
                                   K_CURVE**& curves)
{
  assert(P.num_vars == 2);
//  assert(l_s < h_s);
//  assert(l_t < h_t);
  
  static unsigned int cut_num = 0;
  
  unsigned long i, j, k, l;
  unsigned long num_turn;
  unsigned long num_edge;
  unsigned long num_curves;
  
//  cerr << " kcurve: gen_curve_topo_proto: -------------------- " << endl << flush;
//  cerr << "   P = " << endl << P << endl << flush;
//  cerr << "   [" << l_s << ", " << h_s << "] x [" << l_t << ", " << h_t << "]" << endl << flush;
//  cerr << "   ( " << num_edge_l_s << ", " << num_edge_h_s << ", " << num_edge_l_t << ", " << num_edge_h_t << ", " << num_turn_s << ", " << num_turn_t << " )" << endl << flush;
//  cerr << " -------------------------------------------------- " << endl << flush;
  
  num_turn = num_turn_s + num_turn_t;
  num_edge = num_edge_l_s + num_edge_h_s + num_edge_l_t + num_edge_h_t;
  
//  cerr << "   num_turn = " << num_turn << ", num_edge = " << num_edge << endl << flush;
//  cerr << " -------------------------------------------------- " << endl << flush;
  
  if (num_turn > 1)
  //  If the region contains more than 1 turn-points then split the region s.t.
  //    each subregion will contain less turn-points.
  {
    K_POINT2D**   turn;
    unsigned long d_split;
    bigrational   v_split;
    int           splittable;
    
    turn = new K_POINT2D* [num_turn];  //  num_turn > 1 => turn != 0
    
    for (i = 0; i < num_turn_s; i++)
    {
      turn[i] = turn_s[i];
      turn[i]->ref_count++;
    }
    
    for (i = num_turn_s, j = 0; j < num_turn_t; i++, j++)
    {
      turn[i] = turn_t[j];
      turn[i]->ref_count++;
    }
    
    if (sort_s(turn, num_turn))
    //  Turn-points are sortable in s. Split the region by a vertical line.
    {
      d_split    = 0;
      v_split    = (turn[num_turn / 2 - 1]->get_high_s() +
                    turn[num_turn / 2]->get_low_s()) / 2;
                     //  num_turn >= 2 => num_turn / 2 - 1 > 0
      splittable = 1;
    }
    else if (sort_t(turn, num_turn))
    //  Turn-points are sortable in t. Split the region by a horizontal line.
    {
      d_split    = 1;
      v_split    = (turn[num_turn / 2 - 1]->get_high_t() +
                    turn[num_turn / 2]->get_low_t()) / 2;
                  //  num_turn >= 2 => num_turn / 2 - 1 > 0
      splittable = 1;
    }
    else  //  if (!sort_s(turn, num_turn) && !sort_t(turn, num_turn))
    {
      //  See if trun-points are separable by some vertical line.
      //  Turn-points have already been sorted in t.
      
      for (i = 1; i < num_turn && !turn[i - 1]->compare_t(*turn[i]); i++)
        ;
      
      if (i < num_turn)
      {
        d_split    = 1;
        v_split    = (turn[i - 1]->get_high_t() + turn[i]->get_low_t()) / 2;
                       //  i >= 1 => i - 1 >= 0
        splittable = 1;
      }
      else  //  if (i == num_turn)
      {
        //  See if turn-points are separable by some horizontal line.
        //  Turn-points have already been sorted in t. Re-sort them in s.
        
        sort_s(turn, num_turn);
        
        for (j = 1; j < num_turn && !turn[j - 1]->compare_s(*turn[j]); j++)
          ;
        
        if (j < num_turn)
        {
          d_split    = 0;
          v_split    = (turn[j - 1]->get_high_s() + turn[j]->get_low_s()) / 2;
                         //  j >= 1 => j - 1 >= 0
          splittable = 1;
        }
        else  //  if (j == num_turn)
        {
          //  See if turn-points are separable.
          
          for (k = 1; k < num_turn && turn[k - 1]->equal(*turn[k]); k++)
            ;
          
          if (k < num_turn)
          //  Turn-points are separable. Shrink them and recurse.
          {
            for (l = 0; l < num_turn; l++)
              turn[l]->shrink(shrink_step, shrink_step);
            
//            cerr << " kcurve: gen_curve_topo_proto: recursive call after shrink: 0 " << endl << flush;
            num_curves = gen_curve_topo_proto(P, l_s, h_s, l_t, h_t,
                                              edge_l_s, num_edge_l_s,
                                              edge_h_s, num_edge_h_s,
                                              edge_l_t, num_edge_l_t,
                                              edge_h_t, num_edge_h_t,
                                              turn_s, num_turn_s,
                                              turn_t, num_turn_t,
                                              curves);
          }
          else  //  if (k == num_turn)
          //  All turn-points are equal, i.e.,
          //  P, dP / ds and dP / dt all vanish at *turn[0].
          {
//            cerr << "   All turn-points are equal. " << endl << flush;
//            cerr << "   num_turn = " << num_turn << endl << flush;
//            for (k = 0; k < num_turn; k++)
//              cerr << "   *turn[" << k << "] = " << *turn[k] << endl << flush;
            
            if (num_edge > 0)
            //  P is singular at *turn[0].
            {
//              cerr << "   P = " << endl << P << endl << flush;
//              cerr << "   [" << l_s << ", " << h_s << "] x [" << l_t << ", " << h_t << "]" << endl << flush;
//              cerr << " -------------------------------------------------- " << endl << flush;
              cerr << "   Singularity detected. " << endl << flush;
              abort();
            }
            else  //  if (num_edge == 0)
            //  P = 0 is just a point.
            {
              cerr << "   Degeneracy detected. " << endl << flush;
              
              num_curves = 0;
              curves     = 0;
            }
          }
          
          splittable = 0;
        }
      }
    }
    
    for (i = 0; i < num_turn; i++)
      if (!--turn[i]->ref_count)
        delete turn[i];
    
    delete [] turn;  //  num_turn > 1 => turn != 0
    
    if (splittable && d_split == 0)
    //  Split the region by some horizontal line.
    //  Thus, each subregion will contain less turn-points and edge-t-points.
    {
      K_POINT2D**   split_s;
      unsigned long num_split_s;
      K_POINT2D**   edge_l_t_b;
      K_POINT2D**   edge_h_t_b;
      K_POINT2D**   edge_l_t_a;
      K_POINT2D**   edge_h_t_a;
      unsigned long num_edge_l_t_b, num_edge_h_t_b;
      unsigned long num_edge_h_t_a, num_edge_l_t_a;
      K_POINT2D**   turn_s_b;
      K_POINT2D**   turn_t_b;
      K_POINT2D**   turn_s_a;
      K_POINT2D**   turn_t_a;
      unsigned long num_turn_s_b, num_turn_t_b;
      unsigned long num_turn_s_a, num_turn_t_a;
      K_CURVE**     curves_b;
      K_CURVE**     curves_a;
      unsigned long num_curves_b, num_curves_a, num_curves_ba, num_curves_bb;
      
      if (num_turn_s > 1)
        sort_s(turn_s, num_turn_s);
      
      if (num_turn_t > 1)
        sort_s(turn_t, num_turn_t);
      
      num_split_s = get_pts_proto(v_split, l_t, h_t, P.subst_val(0, v_split),
                                  P, K_RATPOLY(1, 0, v_split),
                                  split_s,
                                  init_tol, 0);
      
      if (num_split_s > 1)
        sort_t(split_s, num_split_s);
      
      //  Divide edge_l_t & edge_h_t.
      
      for (num_edge_l_t_b = 0;
           num_edge_l_t_b < num_edge_l_t
           &&
           edge_l_t[num_edge_l_t_b]->compare_s(v_split) < 0;
           num_edge_l_t_b++)
        ;
      
      if (num_edge_l_t_b > 0)
        edge_l_t_b = new K_POINT2D* [num_edge_l_t_b];
      else  //  if (num_edge_l_t_b == 0)
        edge_l_t_b = 0;
      
      for (i = 0; i < num_edge_l_t_b; i++)
      {
        edge_l_t_b[i] = edge_l_t[i];
        edge_l_t_b[i]->ref_count++;
      }
      
      if ((num_edge_l_t_a = num_edge_l_t - num_edge_l_t_b) > 0)
        edge_l_t_a = new K_POINT2D* [num_edge_l_t_a];
      else  //  if (num_edge_l_t_a == 0)
        edge_l_t_a = 0;
      
      for (i = 0, j = num_edge_l_t_b; i < num_edge_l_t_a; i++, j++)
      {
        edge_l_t_a[i] = edge_l_t[j];
        edge_l_t_a[i]->ref_count++;
      }
      
      assert(j == num_edge_l_t);
      
      for (num_edge_h_t_b = 0;
           num_edge_h_t_b < num_edge_h_t
           &&
           edge_h_t[num_edge_h_t_b]->compare_s(v_split) < 0;
           num_edge_h_t_b++)
        ;
      
      if (num_edge_h_t_b > 0)
        edge_h_t_b = new K_POINT2D* [num_edge_h_t_b];
      else  //  if (num_edge_h_t_b == 0)
        edge_h_t_b = 0;
      
      for (i = 0; i < num_edge_h_t_b; i++)
      {
        edge_h_t_b[i] = edge_h_t[i];
        edge_h_t_b[i]->ref_count++;
      }
      
      if ((num_edge_h_t_a = num_edge_h_t - num_edge_h_t_b) > 0)
        edge_h_t_a = new K_POINT2D* [num_edge_h_t_a];
      else  //  if (num_edge_h_t_a == 0)
        edge_h_t_a = 0;
      
      for (i = 0, j = num_edge_h_t_b; i < num_edge_h_t_a; i++, j++)
      {
        edge_h_t_a[i] = edge_h_t[j];
        edge_h_t_a[i]->ref_count++;
      }
      
      assert(j == num_edge_h_t);
      
      //  Divide turn_s & turn_t.
      
      for (num_turn_s_b = 0;
           num_turn_s_b < num_turn_s
           &&
           turn_s[num_turn_s_b]->compare_s(v_split) < 0;
           num_turn_s_b++)
        ;
      
      if (num_turn_s_b > 0)
        turn_s_b = new K_POINT2D* [num_turn_s_b];
      else  //  if (num_turn_s_b == 0)
        turn_s_b = 0;
      
      for (i = 0; i < num_turn_s_b; i++)
      {
        turn_s_b[i] = turn_s[i];
        turn_s_b[i]->ref_count++;
      }
      
      if ((num_turn_s_a = num_turn_s - num_turn_s_b) > 0)
        turn_s_a = new K_POINT2D* [num_turn_s_a];
      else  //  if (num_turn_s_a == 0)
        turn_s_a = 0;
      
      for (i = 0, j = num_turn_s_b; i < num_turn_s_a; i++, j++)
      {
        turn_s_a[i] = turn_s[j];
        turn_s_a[i]->ref_count++;
      }
      
      assert(j == num_turn_s);
      
      for (num_turn_t_b = 0;
           num_turn_t_b < num_turn_t
           &&
           turn_t[num_turn_t_b]->compare_s(v_split) < 0;
           num_turn_t_b++)
        ;
      
      if (num_turn_t_b > 0)
        turn_t_b = new K_POINT2D* [num_turn_t_b];
      else  //  if (num_turn_t_b == 0)
        turn_t_b = 0;
      
      for (i = 0; i < num_turn_t_b; i++)
      {
        turn_t_b[i] = turn_t[i];
        turn_t_b[i]->ref_count++;
      }
      
      if ((num_turn_t_a = num_turn_t - num_turn_t_b) > 0)
        turn_t_a = new K_POINT2D* [num_turn_t_a];
      else  //  if (num_turn_t_a == 0)
        turn_t_a = 0;
      
      for (i = 0, j = num_turn_t_b; i < num_turn_t_a; i++, j++)
      {
        turn_t_a[i] = turn_t[j];
        turn_t_a[i]->ref_count++;
      }
      
      assert(j == num_turn_t);
      
      //  Make recursive calls to subregions.
      
      if (l_s < v_split)
        num_curves_b = gen_curve_topo_proto(P, l_s, v_split, l_t, h_t,
                                            edge_l_s, num_edge_l_s,
                                            split_s, num_split_s,
                                            edge_l_t_b, num_edge_l_t_b,
                                            edge_h_t_b, num_edge_h_t_b,
                                            turn_s_b, num_turn_s_b,
                                            turn_t_b, num_turn_t_b,
                                            curves_b);
      else  //  if (l_s == v_split)
      {
        num_curves_b = 0;
        curves_b     = 0;
      }
      
      if (h_s > v_split)
        num_curves_a = gen_curve_topo_proto(P, v_split, h_s, l_t, h_t,
                                            split_s, num_split_s,
                                            edge_h_s, num_edge_h_s,
                                            edge_l_t_a, num_edge_l_t_a,
                                            edge_h_t_a, num_edge_h_t_a,
                                            turn_s_a, num_turn_s_a,
                                            turn_t_a, num_turn_t_a,
                                            curves_a);
      else  //  if (h_s == v_split)
      {
        num_curves_a = 0;
        curves_a     = 0;
      }
      
      //  Merge curves.
      
      num_curves_ba = num_curves_bb = 0;
      
      for (i = 0; i < num_split_s; i++)
      {
        for (j = 0;
             j < num_curves_b && !curves_b[j]->is_start_or_end(*split_s[i]);
             j++)
          ;
        
        if (j < num_curves_b)
        {
          for (k = 0;
               k < num_curves_a && !curves_a[k]->is_start_or_end(*split_s[i]);
               k++)
            ;
          
          if (k < num_curves_a)
          {
            curves_b[j]->add_on(*curves_a[k]);
            
            if (k < num_curves_a - 1)
            {
              curves_a[k] = curves_a[num_curves_a - 1];
              curves_a[k]->ref_count++;
            }
            
            num_curves_a--;
            num_curves_ba++;
          }
          else if (j < num_curves_b - 1 && !curves_b[j]->is_closed())
          {
            for (l = j + 1;
                 l < num_curves_b
                 &&
                 !curves_b[l]->is_start_or_end(*split_s[i]);
                 l++)
              ;
            
            if (l < num_curves_b)
            {
              curves_b[j]->add_on(*curves_b[l]);
              
              if (l < num_curves_b - 1)
              {
                curves_b[l] = curves_b[num_curves_b - 1];
                curves_b[l]->ref_count++;
              }
              
              num_curves_b--;
              num_curves_bb++;
            }
          }
        }
      }
      
      if ((num_curves = num_curves_b + num_curves_a) > 0)
        curves = new K_CURVE* [num_curves];
      else  //  if (num_curves == 0)
        curves = 0;
      
      for (i = 0; i < num_curves_b; i++)
      {
        curves[i] = curves_b[i];
        curves[i]->ref_count++;
      }
      
      for (i = num_curves_b, j = 0; j < num_curves_a; i++, j++)
      {
        curves[i] = curves_a[j];
        curves[i]->ref_count++;
      }
      
      assert(i == num_curves);
      
      if (num_curves_b + num_curves_bb > 0)
      {
        for (i = 0; i < num_curves_b + num_curves_bb; i++)
          if (!--curves_b[i]->ref_count)
            delete curves_b[i];
        
        delete [] curves_b;
      }
      
      if (num_curves_a + num_curves_ba > 0)
      {
        for (i = 0; i < num_curves_a + num_curves_ba; i++)
          if (!--curves_a[i]->ref_count)
            delete curves_a[i];
        
        delete [] curves_a;
      }
      
      if (num_split_s > 0)
      {
        for (i = 0; i < num_split_s; i++)
          if (!--split_s[i]->ref_count)
            delete split_s[i];
        
        delete [] split_s;
      }
      
      if (num_edge_l_t_b > 0)
      {
        for (i = 0; i < num_edge_l_t_b; i++)
          if (!--edge_l_t_b[i]->ref_count)
            delete edge_l_t_b[i];
        
        delete [] edge_l_t_b;
      }
      
      if (num_edge_l_t_a > 0)
      {
        for (i = 0; i < num_edge_l_t_a; i++)
          if (!--edge_l_t_a[i]->ref_count)
            delete edge_l_t_a[i];
        
        delete [] edge_l_t_a;
      }
      
      if (num_edge_h_t_b > 0)
      {
        for (i = 0; i < num_edge_h_t_b; i++)
          if (!--edge_h_t_b[i]->ref_count)
            delete edge_h_t_b[i];
        
        delete [] edge_h_t_b;
      }
      
      if (num_edge_h_t_a > 0)
      {
        for (i = 0; i < num_edge_h_t_a; i++)
          if (!--edge_h_t_a[i]->ref_count)
            delete edge_h_t_a[i];
        
        delete [] edge_h_t_a;
      }
      
      if (num_turn_s_b > 0)
      {
        for (i = 0; i < num_turn_s_b; i++)
          if (!--turn_s_b[i]->ref_count)
            delete turn_s_b[i];
        
        delete [] turn_s_b;
      }
      
      if (num_turn_t_b > 0)
      {
        for (i = 0; i < num_turn_t_b; i++)
          if (!--turn_t_b[i]->ref_count)
            delete turn_t_b[i];
        
        delete [] turn_t_b;
      }
      
      if (num_turn_s_a > 0)
      {
        for (i = 0; i < num_turn_s_a; i++)
          if (!--turn_s_a[i]->ref_count)
            delete turn_s_a[i];
        
        delete [] turn_s_a;
      }
      
      if (num_turn_t_a > 0)
      {
        for (i = 0; i < num_turn_t_a; i++)
          if (!--turn_t_a[i]->ref_count)
            delete turn_t_a[i];
        
        delete [] turn_t_a;
      }
    }
    else if (splittable && d_split == 1)
    //  Split the region by some vertical line.
    //  Thus, each subregion will contain less turn-points and edge-s-points.
    {
      K_POINT2D**   split_t;
      unsigned long num_split_t;
      K_POINT2D**   edge_l_s_b;
      K_POINT2D**   edge_h_s_b;
      K_POINT2D**   edge_l_s_a;
      K_POINT2D**   edge_h_s_a;
      unsigned long num_edge_l_s_b, num_edge_h_s_b;
      unsigned long num_edge_h_s_a, num_edge_l_s_a;
      K_POINT2D**   turn_s_b;
      K_POINT2D**   turn_t_b;
      K_POINT2D**   turn_s_a;
      K_POINT2D**   turn_t_a;
      unsigned long num_turn_s_b, num_turn_t_b;
      unsigned long num_turn_s_a, num_turn_t_a;
      K_CURVE**     curves_b;
      K_CURVE**     curves_a;
      unsigned long num_curves_b, num_curves_a, num_curves_ba, num_curves_bb;
      
      if (num_turn_s > 1)
        sort_t(turn_s, num_turn_s);
      
      if (num_turn_t > 1)
        sort_t(turn_t, num_turn_t);
      
      num_split_t = get_pts_proto(l_s, h_s, P.subst_val(1, v_split), v_split,
                                  P, K_RATPOLY(1, 0, v_split),
                                  split_t,
                                  init_tol, 0);
      
      if (num_split_t > 1)
        sort_s(split_t, num_split_t);
      
      //  Divide edge_l_s & edge_h_s.
      
      for (num_edge_l_s_b = 0;
           num_edge_l_s_b < num_edge_l_s
           &&
           edge_l_s[num_edge_l_s_b]->compare_t(v_split) < 0;
           num_edge_l_s_b++)
        ;
      
      if (num_edge_l_s_b > 0)
        edge_l_s_b = new K_POINT2D* [num_edge_l_s_b];
      else  //  if (num_edge_l_s_b == 0)
        edge_l_s_b = 0;
      
      for (i = 0; i < num_edge_l_s_b; i++)
      {
        edge_l_s_b[i] = edge_l_s[i];
        edge_l_s_b[i]->ref_count++;
      }
      
      if ((num_edge_l_s_a = num_edge_l_s - num_edge_l_s_b) > 0)
        edge_l_s_a = new K_POINT2D* [num_edge_l_s_a];
      else  //  if (num_edge_l_s_a == 0)
        edge_l_s_a = 0;
      
      for (i = 0, j = num_edge_l_s_b; i < num_edge_l_s_a; i++, j++)
      {
        edge_l_s_a[i] = edge_l_s[j];
        edge_l_s_a[i]->ref_count++;
      }
      
      assert(j == num_edge_l_s);
      
      for (num_edge_h_s_b = 0;
           num_edge_h_s_b < num_edge_h_s
           &&
           edge_h_s[num_edge_h_s_b]->compare_t(v_split) < 0;
           num_edge_h_s_b++)
        ;
      
      if (num_edge_h_s_b > 0)
        edge_h_s_b = new K_POINT2D* [num_edge_h_s_b];
      else  //  if (num_edge_h_s_b == 0)
        edge_h_s_b = 0;
      
      for (i = 0; i < num_edge_h_s_b; i++)
      {
        edge_h_s_b[i] = edge_h_s[i];
        edge_h_s_b[i]->ref_count++;
      }
      
      if ((num_edge_h_s_a = num_edge_h_s - num_edge_h_s_b) > 0)
        edge_h_s_a = new K_POINT2D* [num_edge_h_s_a];
      else  //  if (num_edge_h_s_a == 0)
        edge_h_s_a = 0;
      
      for (i = 0, j = num_edge_h_s_b; i < num_edge_h_s_a; i++, j++)
      {
        edge_h_s_a[i] = edge_h_s[j];
        edge_h_s_a[i]->ref_count++;
      }
      
      assert(j == num_edge_h_s);
      
      //  Divide turn_s & turn_t.
      
      for (num_turn_s_b = 0;
           num_turn_s_b < num_turn_s
           &&
           turn_s[num_turn_s_b]->compare_t(v_split) < 0;
           num_turn_s_b++)
        ;
      
      if (num_turn_s_b > 0)
        turn_s_b = new K_POINT2D* [num_turn_s_b];
      else  //  if (num_turn_s_b == 0)
        turn_s_b = 0;
      
      for (i = 0; i < num_turn_s_b; i++)
      {
        turn_s_b[i] = turn_s[i];
        turn_s_b[i]->ref_count++;
      }
      
      if ((num_turn_s_a = num_turn_s - num_turn_s_b) > 0)
        turn_s_a = new K_POINT2D* [num_turn_s_a];
      else  //  if (num_turn_s_a == 0)
        turn_s_a = 0;
      
      for (i = 0, j = num_turn_s_b; i < num_turn_s_a; i++, j++)
      {
        turn_s_a[i] = turn_s[j];
        turn_s_a[i]->ref_count++;
      }
      
      assert(j == num_turn_s);
      
      for (num_turn_t_b = 0;
           num_turn_t_b < num_turn_t
           &&
           turn_t[num_turn_t_b]->compare_t(v_split) < 0;
           num_turn_t_b++)
        ;
      
      if (num_turn_t_b > 0)
        turn_t_b = new K_POINT2D* [num_turn_t_b];
      else  //  if (num_turn_t_b == 0)
        turn_t_b = 0;
      
      for (i = 0; i < num_turn_t_b; i++)
      {
        turn_t_b[i] = turn_t[i];
        turn_t_b[i]->ref_count++;
      }
      
      if ((num_turn_t_a = num_turn_t - num_turn_t_b) > 0)
        turn_t_a = new K_POINT2D* [num_turn_t_a];
      else  //  if (num_turn_t_a == 0)
        turn_t_a = 0;
      
      for (i = 0, j = num_turn_t_b; i < num_turn_t_a; i++, j++)
      {
        turn_t_a[i] = turn_t[j];
        turn_t_a[i]->ref_count++;
      }
      
      assert(j == num_turn_t);
      
      //  Make recursive calls to subregions.
      
      if (l_t < v_split)
        num_curves_b = gen_curve_topo_proto(P, l_s, h_s, l_t, v_split,
                                            edge_l_s_b, num_edge_l_s_b,
                                            edge_h_s_b, num_edge_h_s_b,
                                            edge_l_t, num_edge_l_t,
                                            split_t, num_split_t,
                                            turn_s_b, num_turn_s_b,
                                            turn_t_b, num_turn_t_b,
                                            curves_b);
      else  //  if (l_t == v_split)
      {
        num_curves_b = 0;
        curves_b     = 0;
      }
      
      if (h_t > v_split)
        num_curves_a = gen_curve_topo_proto(P, l_s, h_s, v_split, h_t,
                                            edge_l_s_a, num_edge_l_s_a,
                                            edge_h_s_a, num_edge_h_s_a,
                                            split_t, num_split_t,
                                            edge_h_t, num_edge_h_t,
                                            turn_s_a, num_turn_s_a,
                                            turn_t_a, num_turn_t_a,
                                            curves_a);
      else  //  if (h_t == v_split)
      {
        num_curves_a = 0;
        curves_a     = 0;
      }
      
      //  Merge curves.
      
      num_curves_ba = num_curves_bb = 0;
      
      for (i = 0; i < num_split_t; i++)
      {
        for (j = 0;
             j < num_curves_b && !curves_b[j]->is_start_or_end(*split_t[i]);
             j++)
          ;
        
        if (j < num_curves_b)
        {
          for (k = 0;
               k < num_curves_a && !curves_a[k]->is_start_or_end(*split_t[i]);
               k++)
            ;
          
          if (k < num_curves_a)
          {
            curves_b[j]->add_on(*curves_a[k]);
            
            if (k < num_curves_a - 1)
            {
              curves_a[k] = curves_a[num_curves_a - 1];
              curves_a[k]->ref_count++;
            }
            
            num_curves_a--;
            num_curves_ba++;
          }
          else if (j < num_curves_b - 1 && !curves_b[j]->is_closed())
          {
            for (l = j + 1;
                 l < num_curves_b
                 &&
                 !curves_b[l]->is_start_or_end(*split_t[i]);
                 l++)
              ;
            
            if (l < num_curves_b)
            {
              curves_b[j]->add_on(*curves_b[l]);
              
              if (l < num_curves_b - 1)
              {
                curves_b[l] = curves_b[num_curves_b - 1];
                curves_b[l]->ref_count++;
              }
              
              num_curves_b--;
              num_curves_bb++;
            }
          }
        }
      }
      
      if ((num_curves = num_curves_b + num_curves_a) > 0)
        curves = new K_CURVE* [num_curves];
      else  //  if (num_curves == 0)
        curves = 0;
      
      for (i = 0; i < num_curves_b; i++)
      {
        curves[i] = curves_b[i];
        curves[i]->ref_count++;
      }
      
      for (i = num_curves_b, j = 0; j < num_curves_a; i++, j++)
      {
        curves[i] = curves_a[j];
        curves[i]->ref_count++;
      }
      
      assert(i == num_curves);
      
      if (num_curves_b + num_curves_bb > 0)
      {
        for (i = 0; i < num_curves_b + num_curves_bb; i++)
          if (!--curves_b[i]->ref_count)
            delete curves_b[i];
        
        delete [] curves_b;
      }
      
      if (num_curves_a + num_curves_ba > 0)
      {
        for (i = 0; i < num_curves_a + num_curves_ba; i++)
          if (!--curves_a[i]->ref_count)
            delete curves_a[i];
        
        delete [] curves_a;
      }
      
      if (num_split_t > 0)
      {
        for (i = 0; i < num_split_t; i++)
          if (!--split_t[i]->ref_count)
            delete split_t[i];
        
        delete [] split_t;
      }
      
      if (num_edge_l_s_b > 0)
      {
        for (i = 0; i < num_edge_l_s_b; i++)
          if (!--edge_l_s_b[i]->ref_count)
            delete edge_l_s_b[i];
        
        delete [] edge_l_s_b;
      }
      
      if (num_edge_l_s_a > 0)
      {
        for (i = 0; i < num_edge_l_s_a; i++)
          if (!--edge_l_s_a[i]->ref_count)
            delete edge_l_s_a[i];
        
        delete [] edge_l_s_a;
      }
      
      if (num_edge_h_s_b > 0)
      {
        for (i = 0; i < num_edge_h_s_b; i++)
          if (!--edge_h_s_b[i]->ref_count)
            delete edge_h_s_b[i];
        
        delete [] edge_h_s_b;
      }
      
      if (num_edge_h_s_a > 0)
      {
        for (i = 0; i < num_edge_h_s_a; i++)
          if (!--edge_h_s_a[i]->ref_count)
            delete edge_h_s_a[i];
        
        delete [] edge_h_s_a;
      }
      
      if (num_turn_s_b > 0)
      {
        for (i = 0; i < num_turn_s_b; i++)
          if (!--turn_s_b[i]->ref_count)
            delete turn_s_b[i];
        
        delete [] turn_s_b;
      }
      
      if (num_turn_t_b > 0)
      {
        for (i = 0; i < num_turn_t_b; i++)
          if (!--turn_t_b[i]->ref_count)
            delete turn_t_b[i];
        
        delete [] turn_t_b;
      }
      
      if (num_turn_s_a > 0)
      {
        for (i = 0; i < num_turn_s_a; i++)
          if (!--turn_s_a[i]->ref_count)
            delete turn_s_a[i];
        
        delete [] turn_s_a;
      }
      
      if (num_turn_t_a > 0)
      {
        for (i = 0; i < num_turn_t_a; i++)
          if (!--turn_t_a[i]->ref_count)
            delete turn_t_a[i];
        
        delete [] turn_t_a;
      }
    }
  }
  else if (num_turn == 1)  //  if (num_turn_s == 1 || num_turn_t == 1)
    if (num_edge == 2)
    //  If the region contains 1 turn-point and exactly 2 edge-points then
    //  the curve must go from one edge-point
    //                    through the turn-point
    //                    to the other edge-point.
    {
      cut_num = 0;
      
      K_POINT2D* turn;
      K_POINT2D* edge[2];
      K_SEGMENT* s[2];
      
      if (num_turn_s == 1)  //  if (num_turn_s == 1 && num_turn_t == 0)
      {
        turn = turn_s[0];
        turn->ref_count++;
      }
      else  //  if (num_turn_s == 0 && num_turn_t == 1)
      {
        turn = turn_t[0];
        turn->ref_count++;
      }
      
      for (i = 0; i < num_edge_l_s; i++)
      {
        edge[i] = edge_l_s[i];
        edge[i]->ref_count++;
      }
      
      for (j = 0; j < num_edge_h_s; i++, j++)
      {
        edge[i] = edge_h_s[j];
        edge[i]->ref_count++;
      }
      
      for (j = 0; j < num_edge_l_t; i++, j++)
      {
        edge[i] = edge_l_t[j];
        edge[i]->ref_count++;
      }
      
      for (j = 0; j < num_edge_h_t; i++, j++)
      {
        edge[i] = edge_h_t[j];
        edge[i]->ref_count++;
      }
      
      assert(i == 2);
      
      s[0] = new K_SEGMENT(edge[0], turn);
      s[0]->ref_count++;
      s[1] = new K_SEGMENT(turn, edge[1]);
      s[1]->ref_count++;
      
      curves    = new K_CURVE* [num_curves = 1];
      curves[0] = new K_CURVE(P, s, 2);
      curves[0]->ref_count++;
      
      if (!--turn->ref_count)
        delete turn;
      
      for (i = 0; i < 2; i++)
        if (!--edge[i]->ref_count)
          delete edge[i];
    }
    else if (num_edge > 2)
    //  If the region contains 1 turn-point and more than 2 edge-points then
    //  split the region along the boundaries of the box for the turn-point.
    {
      if (cut_num == 0 || cut_num == 2)
      //  Split the region by some vertical line.
      {
        bigrational v_split;
        
        if (cut_num == 0)
          if (num_turn_s == 1)
            v_split = turn_s[0]->get_high_s();
          else  //  if (num_turn_t == 1)
            v_split = turn_t[0]->get_high_s();
        else  //  if (cut_num == 2)
          if (num_turn_s == 1)
            v_split = turn_s[0]->get_low_s();
          else  //  if (num_turn_t == 1)
            v_split = turn_t[0]->get_low_s();
        
        K_POINT2D**   split_s;
        unsigned long num_split_s;
        K_POINT2D**   edge_l_t_b;
        K_POINT2D**   edge_h_t_b;
        K_POINT2D**   edge_l_t_a;
        K_POINT2D**   edge_h_t_a;
        unsigned long num_edge_l_t_b, num_edge_h_t_b;
        unsigned long num_edge_h_t_a, num_edge_l_t_a;
        K_POINT2D**   turn_s_b;
        K_POINT2D**   turn_t_b;
        K_POINT2D**   turn_s_a;
        K_POINT2D**   turn_t_a;
        unsigned long num_turn_s_b, num_turn_t_b;
        unsigned long num_turn_s_a, num_turn_t_a;
        K_CURVE**     curves_b;
        K_CURVE**     curves_a;
        unsigned long num_curves_b, num_curves_a, num_curves_ba, num_curves_bb;
        
        num_split_s = get_pts_proto(v_split, l_t, h_t, P.subst_val(0, v_split),
                                    P, K_RATPOLY(1, 0, v_split),
                                    split_s,
                                    init_tol, 0);
        
        if (num_split_s > 1)
          sort_t(split_s, num_split_s);
        
        //  Divide edge_l_t & edge_h_t.
        
        for (num_edge_l_t_b = 0;
             num_edge_l_t_b < num_edge_l_t
             &&
             edge_l_t[num_edge_l_t_b]->compare_s(v_split) < 0;
             num_edge_l_t_b++)
          ;
        
        if (num_edge_l_t_b > 0)
          edge_l_t_b = new K_POINT2D* [num_edge_l_t_b];
        else  //  if (num_edge_l_t_b == 0)
          edge_l_t_b = 0;
        
        for (i = 0; i < num_edge_l_t_b; i++)
        {
          edge_l_t_b[i] = edge_l_t[i];
          edge_l_t_b[i]->ref_count++;
        }
        
        if ((num_edge_l_t_a = num_edge_l_t - num_edge_l_t_b) > 0)
          edge_l_t_a = new K_POINT2D* [num_edge_l_t_a];
        else  //  if (num_edge_l_t_a == 0)
          edge_l_t_a = 0;
        
        for (i = 0, j = num_edge_l_t_b; i < num_edge_l_t_a; i++, j++)
        {
          edge_l_t_a[i] = edge_l_t[j];
          edge_l_t_a[i]->ref_count++;
        }
        
        assert(j == num_edge_l_t);
        
        for (num_edge_h_t_b = 0;
             num_edge_h_t_b < num_edge_h_t
             &&
             edge_h_t[num_edge_h_t_b]->compare_s(v_split) < 0;
             num_edge_h_t_b++)
          ;
        
        if (num_edge_h_t_b > 0)
          edge_h_t_b = new K_POINT2D* [num_edge_h_t_b];
        else  //  if (num_edge_h_t_b == 0)
          edge_h_t_b = 0;
        
        for (i = 0; i < num_edge_h_t_b; i++)
        {
          edge_h_t_b[i] = edge_h_t[i];
          edge_h_t_b[i]->ref_count++;
        }
        
        if ((num_edge_h_t_a = num_edge_h_t - num_edge_h_t_b) > 0)
          edge_h_t_a = new K_POINT2D* [num_edge_h_t_a];
        else  //  if (num_edge_h_t_a == 0)
          edge_h_t_a = 0;
        
        for (i = 0, j = num_edge_h_t_b; i < num_edge_h_t_a; i++, j++)
        {
          edge_h_t_a[i] = edge_h_t[j];
          edge_h_t_a[i]->ref_count++;
        }
        
        assert(j == num_edge_h_t);
        
        //  Divide turn_s & turn_t.
        
        if (cut_num == 0)
          if (num_turn_s == 1)
          {
            num_turn_s_b = 1;
            num_turn_s_a = num_turn_t_b = num_turn_t_a = 0;
          }
          else  //  if (num_turn_t == 1)
          {
            num_turn_s_b = num_turn_s_a = num_turn_t_a = 0;
            num_turn_t_b = 1;
          }
        else  //  if (cut_num == 2)
          if (num_turn_s == 1)
          {
            num_turn_s_b = num_turn_t_b = num_turn_t_a = 0;
            num_turn_s_a = 1;
          }
          else  //  if (num_turn_t == 1)
          {
            num_turn_s_b = num_turn_s_a = num_turn_t_b = 0;
            num_turn_t_a = 1;
          }
        
        if (num_turn_s_b > 0)
          turn_s_b = new K_POINT2D* [num_turn_s_b];
        else  //  if (num_turn_s_b == 0)
          turn_s_b = 0;
        
        for (i = 0; i < num_turn_s_b; i++)
        {
          turn_s_b[i] = turn_s[i];
          turn_s_b[i]->ref_count++;
        }
        
        if (num_turn_s_a > 0)
          turn_s_a = new K_POINT2D* [num_turn_s_a];
        else  //  if (num_turn_s_a == 0)
          turn_s_a = 0;
        
        for (i = 0, j = num_turn_s_b; i < num_turn_s_a; i++, j++)
        {
          turn_s_a[i] = turn_s[j];
          turn_s_a[i]->ref_count++;
        }
        
        assert(j == num_turn_s);
        
        if (num_turn_t_b > 0)
          turn_t_b = new K_POINT2D* [num_turn_t_b];
        else  //  if (num_turn_t_b == 0)
          turn_t_b = 0;
        
        for (i = 0; i < num_turn_t_b; i++)
        {
          turn_t_b[i] = turn_t[i];
          turn_t_b[i]->ref_count++;
        }
        
        if (num_turn_t_a > 0)
          turn_t_a = new K_POINT2D* [num_turn_t_a];
        else  //  if (num_turn_t_a == 0)
          turn_t_a = 0;
        
        for (i = 0, j = num_turn_t_b; i < num_turn_t_a; i++, j++)
        {
          turn_t_a[i] = turn_t[j];
          turn_t_a[i]->ref_count++;
        }
        
        assert(j == num_turn_t);
        
        //  Make recursive calls to subregions.
        
        cut_num++;
        
        if (l_s < v_split)
          num_curves_b = gen_curve_topo_proto(P, l_s, v_split, l_t, h_t,
                                              edge_l_s, num_edge_l_s,
                                              split_s, num_split_s,
                                              edge_l_t_b, num_edge_l_t_b,
                                              edge_h_t_b, num_edge_h_t_b,
                                              turn_s_b, num_turn_s_b,
                                              turn_t_b, num_turn_t_b,
                                              curves_b);
        else  //  if (l_s == v_split)
        {
          num_curves_b = 0;
          curves_b     = 0;
        }
        
        if (h_s > v_split)
          num_curves_a = gen_curve_topo_proto(P, v_split, h_s, l_t, h_t,
                                              split_s, num_split_s,
                                              edge_h_s, num_edge_h_s,
                                              edge_l_t_a, num_edge_l_t_a,
                                              edge_h_t_a, num_edge_h_t_a,
                                              turn_s_a, num_turn_s_a,
                                              turn_t_a, num_turn_t_a,
                                              curves_a);
        else  //  if (h_s == v_split)
        {
          num_curves_a = 0;
          curves_a     = 0;
        }
        
        //  Merge curves.
        
        num_curves_ba = num_curves_bb = 0;
        
        for (i = 0; i < num_split_s; i++)
        {
          for (j = 0;
               j < num_curves_b && !curves_b[j]->is_start_or_end(*split_s[i]);
               j++)
            ;
          
          if (j < num_curves_b)
          {
            for (k = 0;
                 k < num_curves_a
                 &&
                 !curves_a[k]->is_start_or_end(*split_s[i]);
                 k++)
              ;
            
            if (k < num_curves_a)
            {
              curves_b[j]->add_on(*curves_a[k]);
              
              if (k < num_curves_a - 1)
              {
                curves_a[k] = curves_a[num_curves_a - 1];
                curves_a[k]->ref_count++;
              }
              
              num_curves_a--;
              num_curves_ba++;
            }
            else if (j < num_curves_b - 1 && !curves_b[j]->is_closed())
            {
              for (l = j + 1;
                   l < num_curves_b
                   &&
                   !curves_b[l]->is_start_or_end(*split_s[i]);
                   l++)
                ;
              
              if (l < num_curves_b)
              {
                curves_b[j]->add_on(*curves_b[l]);
                
                if (l < num_curves_b - 1)
                {
                  curves_b[l] = curves_b[num_curves_b - 1];
                  curves_b[l]->ref_count++;
                }
                
                num_curves_b--;
                num_curves_bb++;
              }
            }
          }
        }
        
        if ((num_curves = num_curves_b + num_curves_a) > 0)
          curves = new K_CURVE* [num_curves];
        else  //  if (num_curves == 0)
          curves = 0;
        
        for (i = 0; i < num_curves_b; i++)
        {
          curves[i] = curves_b[i];
          curves[i]->ref_count++;
        }
        
        for (i = num_curves_b, j = 0; j < num_curves_a; i++, j++)
        {
          curves[i] = curves_a[j];
          curves[i]->ref_count++;
        }
        
        assert(i == num_curves);
        
        if (num_curves_b + num_curves_bb > 0)
        {
          for (i = 0; i < num_curves_b + num_curves_bb; i++)
            if (!--curves_b[i]->ref_count)
              delete curves_b[i];
          
          delete [] curves_b;
        }
        
        if (num_curves_a + num_curves_ba > 0)
        {
          for (i = 0; i < num_curves_a + num_curves_ba; i++)
            if (!--curves_a[i]->ref_count)
              delete curves_a[i];
          
          delete [] curves_a;
        }
        
        if (num_split_s > 0)
        {
          for (i = 0; i < num_split_s; i++)
            if (!--split_s[i]->ref_count)
              delete split_s[i];
          
          delete [] split_s;
        }
        
        if (num_edge_l_t_b > 0)
        {
          for (i = 0; i < num_edge_l_t_b; i++)
            if (!--edge_l_t_b[i]->ref_count)
              delete edge_l_t_b[i];
          
          delete [] edge_l_t_b;
        }
        
        if (num_edge_l_t_a > 0)
        {
          for (i = 0; i < num_edge_l_t_a; i++)
            if (!--edge_l_t_a[i]->ref_count)
              delete edge_l_t_a[i];
          
          delete [] edge_l_t_a;
        }
        
        if (num_edge_h_t_b > 0)
        {
          for (i = 0; i < num_edge_h_t_b; i++)
            if (!--edge_h_t_b[i]->ref_count)
              delete edge_h_t_b[i];
          
          delete [] edge_h_t_b;
        }
        
        if (num_edge_h_t_a > 0)
        {
          for (i = 0; i < num_edge_h_t_a; i++)
            if (!--edge_h_t_a[i]->ref_count)
              delete edge_h_t_a[i];
          
          delete [] edge_h_t_a;
        }
        
        if (num_turn_s_b > 0)
        {
          for (i = 0; i < num_turn_s_b; i++)
            if (!--turn_s_b[i]->ref_count)
              delete turn_s_b[i];
          
          delete [] turn_s_b;
        }
        
        if (num_turn_t_b > 0)
        {
          for (i = 0; i < num_turn_t_b; i++)
            if (!--turn_t_b[i]->ref_count)
              delete turn_t_b[i];
          
          delete [] turn_t_b;
        }
        
        if (num_turn_s_a > 0)
        {
          for (i = 0; i < num_turn_s_a; i++)
            if (!--turn_s_a[i]->ref_count)
              delete turn_s_a[i];
          
          delete [] turn_s_a;
        }
        
        if (num_turn_t_a > 0)
        {
          for (i = 0; i < num_turn_t_a; i++)
            if (!--turn_t_a[i]->ref_count)
              delete turn_t_a[i];
          
          delete [] turn_t_a;
        }
      }
      else if (cut_num == 1 || cut_num == 3)
      //  Split the region by some horizontal line.
      {
        bigrational v_split;
        
        if (cut_num == 1)
          if (num_turn_s == 1)
            v_split = turn_s[0]->get_high_t();
          else  //  if (num_turn_t == 1)
            v_split = turn_t[0]->get_high_t();
        else  //  if (cut_num == 3)
          if (num_turn_s == 1)
            v_split = turn_s[0]->get_low_t();
          else  //  if (num_turn_t == 1)
            v_split = turn_t[0]->get_low_t();
        
        K_POINT2D**   split_t;
        unsigned long num_split_t;
        K_POINT2D**   edge_l_s_b;
        K_POINT2D**   edge_h_s_b;
        K_POINT2D**   edge_l_s_a;
        K_POINT2D**   edge_h_s_a;
        unsigned long num_edge_l_s_b, num_edge_h_s_b;
        unsigned long num_edge_h_s_a, num_edge_l_s_a;
        K_POINT2D**   turn_s_b;
        K_POINT2D**   turn_t_b;
        K_POINT2D**   turn_s_a;
        K_POINT2D**   turn_t_a;
        unsigned long num_turn_s_b, num_turn_t_b;
        unsigned long num_turn_s_a, num_turn_t_a;
        K_CURVE**     curves_b;
        K_CURVE**     curves_a;
        unsigned long num_curves_b, num_curves_a, num_curves_ba, num_curves_bb;
        
        num_split_t = get_pts_proto(l_s, h_s, P.subst_val(1, v_split), v_split,
                                    P, K_RATPOLY(1, 0, v_split),
                                    split_t,
                                    init_tol, 0);
        
        if (num_split_t > 1)
          sort_s(split_t, num_split_t);
        
        //  Divide edge_l_s & edge_h_s.
        
        for (num_edge_l_s_b = 0;
             num_edge_l_s_b < num_edge_l_s
             &&
             edge_l_s[num_edge_l_s_b]->compare_t(v_split) < 0;
             num_edge_l_s_b++)
          ;
        
        if (num_edge_l_s_b > 0)
          edge_l_s_b = new K_POINT2D* [num_edge_l_s_b];
        else  //  if (num_edge_l_s_b == 0)
          edge_l_s_b = 0;
        
        for (i = 0; i < num_edge_l_s_b; i++)
        {
          edge_l_s_b[i] = edge_l_s[i];
          edge_l_s_b[i]->ref_count++;
        }
        
        if ((num_edge_l_s_a = num_edge_l_s - num_edge_l_s_b) > 0)
          edge_l_s_a = new K_POINT2D* [num_edge_l_s_a];
        else  //  if (num_edge_l_s_a == 0)
          edge_l_s_a = 0;
        
        for (i = 0, j = num_edge_l_s_b; i < num_edge_l_s_a; i++, j++)
        {
          edge_l_s_a[i] = edge_l_s[j];
          edge_l_s_a[i]->ref_count++;
        }
        
        assert(j == num_edge_l_s);
        
        for (num_edge_h_s_b = 0;
             num_edge_h_s_b < num_edge_h_s
             &&
             edge_h_s[num_edge_h_s_b]->compare_t(v_split) < 0;
             num_edge_h_s_b++)
          ;
        
        if (num_edge_h_s_b > 0)
          edge_h_s_b = new K_POINT2D* [num_edge_h_s_b];
        else  //  if (num_edge_h_s_b == 0)
          edge_h_s_b = 0;
        
        for (i = 0; i < num_edge_h_s_b; i++)
        {
          edge_h_s_b[i] = edge_h_s[i];
          edge_h_s_b[i]->ref_count++;
        }
        
        if ((num_edge_h_s_a = num_edge_h_s - num_edge_h_s_b) > 0)
          edge_h_s_a = new K_POINT2D* [num_edge_h_s_a];
        else  //  if (num_edge_h_s_a == 0)
          edge_h_s_a = 0;
        
        for (i = 0, j = num_edge_h_s_b; i < num_edge_h_s_a; i++, j++)
        {
          edge_h_s_a[i] = edge_h_s[j];
          edge_h_s_a[i]->ref_count++;
        }
        
        assert(j == num_edge_h_s);
        
        //  Divide turn_s & turn_t.
        
        if (cut_num == 1)
          if (num_turn_s == 1)
          {
            num_turn_s_b = 1;
            num_turn_s_a = num_turn_t_b = num_turn_t_a = 0;
          }
          else  //  if (num_turn_t == 1)
          {
            num_turn_s_b = num_turn_s_a = num_turn_t_a = 0;
            num_turn_t_b = 1;
          }
        else  //  if (cut_num == 3)
          if (num_turn_s == 1)
          {
            num_turn_s_b = num_turn_t_b = num_turn_t_a = 0;
            num_turn_s_a = 1;
          }
          else  //  if (num_turn_t == 1)
          {
            num_turn_s_b = num_turn_s_a = num_turn_t_b = 0;
            num_turn_t_a = 1;
          }
        
        if (num_turn_s_b > 0)
          turn_s_b = new K_POINT2D* [num_turn_s_b];
        else  //  if (num_turn_s_b == 0)
          turn_s_b = 0;
        
        for (i = 0; i < num_turn_s_b; i++)
        {
          turn_s_b[i] = turn_s[i];
          turn_s_b[i]->ref_count++;
        }
        
        if (num_turn_s_a > 0)
          turn_s_a = new K_POINT2D* [num_turn_s_a];
        else  //  if (num_turn_s_a == 0)
          turn_s_a = 0;
        
        for (i = 0, j = num_turn_s_b; i < num_turn_s_a; i++, j++)
        {
          turn_s_a[i] = turn_s[j];
          turn_s_a[i]->ref_count++;
        }
        
        assert(j == num_turn_s);
        
        if (num_turn_t_b > 0)
          turn_t_b = new K_POINT2D* [num_turn_t_b];
        else  //  if (num_turn_t_b == 0)
          turn_t_b = 0;
        
        for (i = 0; i < num_turn_t_b; i++)
        {
          turn_t_b[i] = turn_t[i];
          turn_t_b[i]->ref_count++;
        }
        
        if (num_turn_t_a > 0)
          turn_t_a = new K_POINT2D* [num_turn_t_a];
        else  //  if (num_turn_t_a == 0)
          turn_t_a = 0;
        
        for (i = 0, j = num_turn_t_b; i < num_turn_t_a; i++, j++)
        {
          turn_t_a[i] = turn_t[j];
          turn_t_a[i]->ref_count++;
        }
        
        assert(j == num_turn_t);
        
        //  Make recursive calls to subregions.
        
        cut_num++;
        
        if (l_t < v_split)
          num_curves_b = gen_curve_topo_proto(P, l_s, h_s, l_t, v_split,
                                              edge_l_s_b, num_edge_l_s_b,
                                              edge_h_s_b, num_edge_h_s_b,
                                              edge_l_t, num_edge_l_t,
                                              split_t, num_split_t,
                                              turn_s_b, num_turn_s_b,
                                              turn_t_b, num_turn_t_b,
                                              curves_b);
        else  //  if (l_t == v_split)
        {
          num_curves_b = 0;
          curves_b     = 0;
        }
        
        if (h_t > v_split)
          num_curves_a = gen_curve_topo_proto(P, l_s, h_s, v_split, h_t,
                                              edge_l_s_a, num_edge_l_s_a,
                                              edge_h_s_a, num_edge_h_s_a,
                                              split_t, num_split_t,
                                              edge_h_t, num_edge_h_t,
                                              turn_s_a, num_turn_s_a,
                                              turn_t_a, num_turn_t_a,
                                              curves_a);
        else  //  if (h_t == v_split)
        {
          num_curves_a = 0;
          curves_a     = 0;
        }
        
        //  Merge curves.
        
        num_curves_ba = num_curves_bb = 0;
        
        for (i = 0; i < num_split_t; i++)
        {
          for (j = 0;
               j < num_curves_b && !curves_b[j]->is_start_or_end(*split_t[i]);
               j++)
            ;
          
          if (j < num_curves_b)
          {
            for (k = 0;
                 k < num_curves_a
                 &&
                 !curves_a[k]->is_start_or_end(*split_t[i]);
                 k++)
              ;
            
            if (k < num_curves_a)
            {
              curves_b[j]->add_on(*curves_a[k]);
              
              if (k < num_curves_a - 1)
              {
                curves_a[k] = curves_a[num_curves_a - 1];
                curves_a[k]->ref_count++;
              }
              
              num_curves_a--;
              num_curves_ba++;
            }
            else if (j < num_curves_b - 1 && !curves_b[j]->is_closed())
            {
              for (l = j + 1;
                   l < num_curves_b
                   &&
                   !curves_b[l]->is_start_or_end(*split_t[i]);
                   l++)
                ;
              
              if (l < num_curves_b)
              {
                curves_b[j]->add_on(*curves_b[l]);
                
                if (l < num_curves_b - 1)
                {
                  curves_b[l] = curves_b[num_curves_b - 1];
                  curves_b[l]->ref_count++;
                }
                
                num_curves_b--;
                num_curves_bb++;
              }
            }
          }
        }
        
        if ((num_curves = num_curves_b + num_curves_a) > 0)
          curves = new K_CURVE* [num_curves];
        else  //  if (num_curves == 0)
          curves = 0;
        
        for (i = 0; i < num_curves_b; i++)
        {
          curves[i] = curves_b[i];
          curves[i]->ref_count++;
        }
        
        for (i = num_curves_b, j = 0; j < num_curves_a; i++, j++)
        {
          curves[i] = curves_a[j];
          curves[i]->ref_count++;
        }
        
        assert(i == num_curves);
        
        if (num_curves_b + num_curves_bb > 0)
        {
          for (i = 0; i < num_curves_b + num_curves_bb; i++)
            if (!--curves_b[i]->ref_count)
              delete curves_b[i];
          
          delete [] curves_b;
        }
        
        if (num_curves_a + num_curves_ba > 0)
        {
          for (i = 0; i < num_curves_a + num_curves_ba; i++)
            if (!--curves_a[i]->ref_count)
              delete curves_a[i];
          
          delete [] curves_a;
        }
        
        if (num_split_t > 0)
        {
          for (i = 0; i < num_split_t; i++)
            if (!--split_t[i]->ref_count)
              delete split_t[i];
          
          delete [] split_t;
        }
        
        if (num_edge_l_s_b > 0)
        {
          for (i = 0; i < num_edge_l_s_b; i++)
            if (!--edge_l_s_b[i]->ref_count)
              delete edge_l_s_b[i];
          
          delete [] edge_l_s_b;
        }
        
        if (num_edge_l_s_a > 0)
        {
          for (i = 0; i < num_edge_l_s_a; i++)
            if (!--edge_l_s_a[i]->ref_count)
              delete edge_l_s_a[i];
          
          delete [] edge_l_s_a;
        }
        
        if (num_edge_h_s_b > 0)
        {
          for (i = 0; i < num_edge_h_s_b; i++)
            if (!--edge_h_s_b[i]->ref_count)
              delete edge_h_s_b[i];
          
          delete [] edge_h_s_b;
        }
        
        if (num_edge_h_s_a > 0)
        {
          for (i = 0; i < num_edge_h_s_a; i++)
            if (!--edge_h_s_a[i]->ref_count)
              delete edge_h_s_a[i];
          
          delete [] edge_h_s_a;
        }
        
        if (num_turn_s_b > 0)
        {
          for (i = 0; i < num_turn_s_b; i++)
            if (!--turn_s_b[i]->ref_count)
              delete turn_s_b[i];
          
          delete [] turn_s_b;
        }
        
        if (num_turn_t_b > 0)
        {
          for (i = 0; i < num_turn_t_b; i++)
            if (!--turn_t_b[i]->ref_count)
              delete turn_t_b[i];
          
          delete [] turn_t_b;
        }
        
        if (num_turn_s_a > 0)
        {
          for (i = 0; i < num_turn_s_a; i++)
            if (!--turn_s_a[i]->ref_count)
              delete turn_s_a[i];
          
          delete [] turn_s_a;
        }
        
        if (num_turn_t_a > 0)
        {
          for (i = 0; i < num_turn_t_a; i++)
            if (!--turn_t_a[i]->ref_count)
              delete turn_t_a[i];
          
          delete [] turn_t_a;
        }
      }
      else
      //  if (cut_num != 0 || cut_num != 1 || cut_num != 2 || cut_num != 3)
      //  All 4 boundaries of the box for the turn-point have been used.
      //  Shrink the turn-point and try again.
      {
        cut_num = 0;
        
        if (num_turn_s == 1)
          turn_s[0]->shrink(shrink_step, shrink_step);
        else  //  if (num_turn_t == 1)
          turn_t[0]->shrink(shrink_step, shrink_step);
        
//        cerr << " kcurve: gen_curve_topo_proto: recursive call after shrink: 1 " << endl << flush;
        num_curves = gen_curve_topo_proto(P, l_s, h_s, l_t, h_t,
                                          edge_l_s, num_edge_l_s,
                                          edge_h_s, num_edge_h_s,
                                          edge_l_t, num_edge_l_t,
                                          edge_h_t, num_edge_h_t,
                                          turn_s, num_turn_s,
                                          turn_t, num_turn_t,
                                          curves);
      }
    }
    else  //  if (num_edge < 2)
    //  If the region contains 1 turn-point then
    //    it must contain at least 2 edge-points
    //  since P is assumed to be regular.
    {
//      cerr << "   P = " << endl << P << endl << flush;
//      cerr << "   [" << l_s << ", " << h_s << "] x [" << l_t << ", " << h_t << "]" << endl << flush;
//      cerr << " -------------------------------------------------- " << endl << flush;
      cerr << "   Singularity detected. " << endl << flush;
      abort();
    }
  else  //  if (!num_turn_s && !num_turn_t)
    if (num_edge == 0)
    {
      curves     = 0;
      num_curves = 0;
    }
    else if (num_edge == 1)
    //  P = 0 is just a point.
    {
      cerr << "   Degeneracy detected. " << endl << flush;
      
      curves     = 0;
      num_curves = 0;
    }
    //  If the region contains no turn-point and exactly 2 edge-points then
    //  the curve must go from one edge-point to the other.
    else if (num_edge_l_s == 1 && num_edge_h_s == 1 &&
             !num_edge_l_t && !num_edge_h_t)
    {
      K_SEGMENT* s[1];
      
      curves    = new K_CURVE* [num_curves = 1];
      s[0]      = new K_SEGMENT(edge_l_s[0], edge_h_s[0]);
      s[0]->ref_count++;
      curves[0] = new K_CURVE(P, s, 1);
      curves[0]->ref_count++;
    }
    else if (num_edge_l_s == 1 && !num_edge_h_s &&
             num_edge_l_t == 1 && !num_edge_h_t)
    {
      K_SEGMENT* s[1];
      
      curves    = new K_CURVE* [num_curves = 1];
      s[0]      = new K_SEGMENT(edge_l_t[0], edge_l_s[0]);
      s[0]->ref_count++;
      curves[0] = new K_CURVE(P, s, 1);
      curves[0]->ref_count++;
    }
    else if (num_edge_l_s == 1 && !num_edge_h_s &&
             !num_edge_l_t && num_edge_h_t == 1)
    {
      K_SEGMENT* s[1];
      
      curves    = new K_CURVE* [num_curves = 1];
      s[0]      = new K_SEGMENT(edge_l_s[0], edge_h_t[0]);
      s[0]->ref_count++;
      curves[0] = new K_CURVE(P, s, 1);
      curves[0]->ref_count++;
    }
    else if (!num_edge_l_s && num_edge_h_s == 1 &&
             num_edge_l_t == 1 && !num_edge_h_t)
    {
      K_SEGMENT* s[1];
      
      curves    = new K_CURVE* [num_curves = 1];
      s[0]      = new K_SEGMENT(edge_l_t[0], edge_h_s[0]);
      s[0]->ref_count++;
      curves[0] = new K_CURVE(P, s, 1);
      curves[0]->ref_count++;
    }
    else if (!num_edge_l_s && num_edge_h_s == 1 &&
             !num_edge_l_t && num_edge_h_t == 1)
    {
      K_SEGMENT* s[1];
      
      curves    = new K_CURVE* [num_curves = 1];
      s[0]      = new K_SEGMENT(edge_h_s[0], edge_h_t[0]);
      s[0]->ref_count++;
      curves[0] = new K_CURVE(P, s, 1);
      curves[0]->ref_count++;
    }
    else if (!num_edge_l_s && !num_edge_h_s &&
             num_edge_l_t == 1 && num_edge_h_t == 1)
    {
      K_SEGMENT* s[1];
      
      curves    = new K_CURVE* [num_curves = 1];
      s[0]      = new K_SEGMENT(edge_l_t[0], edge_h_t[0]);
      s[0]->ref_count++;
      curves[0] = new K_CURVE(P, s, 1);
      curves[0]->ref_count++;
    }
    else if (num_edge % 2 == 0
             &&
             (num_edge_l_t > num_edge_h_t
              &&
              num_edge_l_s + num_edge_h_s == num_edge_l_t - num_edge_h_t
              ||
              num_edge_l_t <= num_edge_h_t
              &&
              num_edge_l_s + num_edge_h_s == num_edge_h_t - num_edge_l_t
              ||
              num_edge_l_s + num_edge_h_s > num_edge_l_t + num_edge_h_t
              &&
              (num_edge_l_s > num_edge_h_s
               &&
               num_edge_l_s - num_edge_h_s != num_edge_l_t + num_edge_h_t
               ||
               num_edge_l_s <= num_edge_h_s
               &&
               num_edge_h_s - num_edge_l_s != num_edge_l_t + num_edge_h_t)))
    //  Split the region by some vertical line.
    //  Thus, each subregion will contain less t-edge-points.
    {
      bigrational v_split;
      
      if (num_edge_l_t > 1)
        v_split = (edge_l_t[num_edge_l_t / 2 - 1]->get_high_s() +
                   edge_l_t[num_edge_l_t / 2]->get_low_s()) / 2;
      else if (num_edge_h_t > 1)
        v_split = (edge_h_t[num_edge_h_t / 2 - 1]->get_high_s() +
                   edge_h_t[num_edge_h_t / 2]->get_low_s()) / 2;
      else
      {
//        cerr << " kcurve: gen_curve_topo_proto: num_edge_l_t = " << num_edge_l_t << ", num_edge_h_t = " << num_edge_h_t << endl << flush;
        assert(num_edge_l_t >= 1 && num_edge_h_t >= 1);
        
        v_split = (edge_l_t[0]->get_low_s() + edge_l_t[0]->get_high_s() +
                   edge_h_t[0]->get_low_s() + edge_h_t[0]->get_high_s()) / 4;
      }
      
      K_POINT2D**   split_s;
      unsigned long num_split_s;
      K_POINT2D**   edge_l_t_b;
      K_POINT2D**   edge_h_t_b;
      K_POINT2D**   edge_l_t_a;
      K_POINT2D**   edge_h_t_a;
      unsigned long num_edge_l_t_b, num_edge_h_t_b;
      unsigned long num_edge_h_t_a, num_edge_l_t_a;
      K_CURVE**     curves_b;
      K_CURVE**     curves_a;
      unsigned long num_curves_b, num_curves_a, num_curves_ba, num_curves_bb;
      
      num_split_s = get_pts_proto(v_split, l_t, h_t, P.subst_val(0, v_split),
                                  P, K_RATPOLY(1, 0, v_split),
                                  split_s,
                                  init_tol, 0);
      
      if (num_split_s > 1)
        sort_t(split_s, num_split_s);
      
      for (num_edge_l_t_b = 0;
           num_edge_l_t_b < num_edge_l_t
           &&
           edge_l_t[num_edge_l_t_b]->compare_s(v_split) < 0;
           num_edge_l_t_b++)
        ;
      
      if (num_edge_l_t_b > 0)
        edge_l_t_b = new K_POINT2D* [num_edge_l_t_b];
      else  //  if (num_edge_l_t_b == 0)
        edge_l_t_b = 0;
      
      for (i = 0; i < num_edge_l_t_b; i++)
      {
        edge_l_t_b[i] = edge_l_t[i];
        edge_l_t_b[i]->ref_count++;
      }
      
      if ((num_edge_l_t_a = num_edge_l_t - num_edge_l_t_b) > 0)
        edge_l_t_a = new K_POINT2D* [num_edge_l_t_a];
      else  //  if (num_edge_l_t_a == 0)
        edge_l_t_a = 0;
      
      for (i = 0, j = num_edge_l_t_b; i < num_edge_l_t_a; i++, j++)
      {
        edge_l_t_a[i] = edge_l_t[j];
        edge_l_t_a[i]->ref_count++;
      }
      
      assert(j == num_edge_l_t);
      
      for (num_edge_h_t_b = 0;
           num_edge_h_t_b < num_edge_h_t
           &&
           edge_h_t[num_edge_h_t_b]->compare_s(v_split) < 0;
           num_edge_h_t_b++)
        ;
      
      if (num_edge_h_t_b > 0)
        edge_h_t_b = new K_POINT2D* [num_edge_h_t_b];
      else  //  if (num_edge_h_t_b == 0)
        edge_h_t_b = 0;
      
      for (i = 0; i < num_edge_h_t_b; i++)
      {
        edge_h_t_b[i] = edge_h_t[i];
        edge_h_t_b[i]->ref_count++;
      }
      
      if ((num_edge_h_t_a = num_edge_h_t - num_edge_h_t_b) > 0)
        edge_h_t_a = new K_POINT2D* [num_edge_h_t_a];
      else  //  if (num_edge_h_t_a == 0)
        edge_h_t_a = 0;
      
      for (i = 0, j = num_edge_h_t_b; i < num_edge_h_t_a; i++, j++)
      {
        edge_h_t_a[i] = edge_h_t[j];
        edge_h_t_a[i]->ref_count++;
      }
      
      assert(j == num_edge_h_t);
      
      if (l_s < v_split)
        num_curves_b = gen_curve_topo_proto(P, l_s, v_split, l_t, h_t,
                                            edge_l_s, num_edge_l_s,
                                            split_s, num_split_s,
                                            edge_l_t_b, num_edge_l_t_b,
                                            edge_h_t_b, num_edge_h_t_b,
                                            turn_s, num_turn_s,
                                            turn_t, num_turn_t,
                                            curves_b);
      else  //  if (l_s == v_split)
      {
        num_curves_b = 0;
        curves_b     = 0;
      }
      
      if (h_s > v_split)
        num_curves_a = gen_curve_topo_proto(P, v_split, h_s, l_t, h_t,
                                            split_s, num_split_s,
                                            edge_h_s, num_edge_h_s,
                                            edge_l_t_a, num_edge_l_t_a,
                                            edge_h_t_a, num_edge_h_t_a,
                                            turn_s, num_turn_s,
                                            turn_t, num_turn_t,
                                            curves_a);
      else  //  if (h_s == v_split)
      {
        num_curves_a = 0;
        curves_a     = 0;
      }
      
      //  Merge curves.
      
      num_curves_ba = num_curves_bb = 0;
      
      for (i = 0; i < num_split_s; i++)
      {
        for (j = 0;
             j < num_curves_b && !curves_b[j]->is_start_or_end(*split_s[i]);
             j++)
          ;
        
        if (j < num_curves_b)
        {
          for (k = 0;
               k < num_curves_a && !curves_a[k]->is_start_or_end(*split_s[i]);
               k++)
            ;
          
          if (k < num_curves_a)
          {
            curves_b[j]->add_on(*curves_a[k]);
            
            if (k < num_curves_a - 1)
            {
              curves_a[k] = curves_a[num_curves_a - 1];
              curves_a[k]->ref_count++;
            }
            
            num_curves_a--;
            num_curves_ba++;
          }
          else if (j < num_curves_b - 1 && !curves_b[j]->is_closed())
          {
            for (l = j + 1;
                 l < num_curves_b
                 &&
                 !curves_b[l]->is_start_or_end(*split_s[i]);
                 l++)
              ;
            
            if (l < num_curves_b)
            {
              curves_b[j]->add_on(*curves_b[l]);
              
              if (l < num_curves_b - 1)
              {
                curves_b[l] = curves_b[num_curves_b - 1];
                curves_b[l]->ref_count++;
              }
              
              num_curves_b--;
              num_curves_bb++;
            }
          }
        }
      }
      
      if ((num_curves = num_curves_b + num_curves_a) > 0)
        curves = new K_CURVE* [num_curves];
      else  //  if (num_curves == 0)
        curves = 0;
      
      for (i = 0; i < num_curves_b; i++)
      {
        curves[i] = curves_b[i];
        curves[i]->ref_count++;
      }
      
      for (i = num_curves_b, j = 0; j < num_curves_a; i++, j++)
      {
        curves[i] = curves_a[j];
        curves[i]->ref_count++;
      }
      
      assert(i == num_curves);
      
      if (num_curves_b + num_curves_bb > 0)
      {
        for (i = 0; i < num_curves_b + num_curves_bb; i++)
          if (!--curves_b[i]->ref_count)
            delete curves_b[i];
        
        delete [] curves_b;
      }
      
      if (num_curves_a + num_curves_ba > 0)
      {
        for (i = 0; i < num_curves_a + num_curves_ba; i++)
          if (!--curves_a[i]->ref_count)
            delete curves_a[i];
        
        delete [] curves_a;
      }
      
      if (num_split_s > 0)
      {
        for (i = 0; i < num_split_s; i++)
          if (!--split_s[i]->ref_count)
            delete split_s[i];
        
        delete [] split_s;
      }
      
      if (num_edge_l_t_b > 0)
      {
        for (i = 0; i < num_edge_l_t_b; i++)
          if (!--edge_l_t_b[i]->ref_count)
            delete edge_l_t_b[i];
        
        delete [] edge_l_t_b;
      }
      
      if (num_edge_l_t_a > 0)
      {
        for (i = 0; i < num_edge_l_t_a; i++)
          if (!--edge_l_t_a[i]->ref_count)
            delete edge_l_t_a[i];
        
        delete [] edge_l_t_a;
      }
      
      if (num_edge_h_t_b > 0)
      {
        for (i = 0; i < num_edge_h_t_b; i++)
          if (!--edge_h_t_b[i]->ref_count)
            delete edge_h_t_b[i];
        
        delete [] edge_h_t_b;
      }
      
      if (num_edge_h_t_a > 0)
      {
        for (i = 0; i < num_edge_h_t_a; i++)
          if (!--edge_h_t_a[i]->ref_count)
            delete edge_h_t_a[i];
        
        delete [] edge_h_t_a;
      }
    }
    else if (num_edge % 2 == 0)
    //  Split the region by some horozontal line.
    //  Thus, each subregion will contain less s-edge-points.
    {
      bigrational v_split;
      
      if (num_edge_l_s > 1)
        v_split = (edge_l_s[num_edge_l_s / 2 - 1]->get_high_t() +
                   edge_l_s[num_edge_l_s / 2]->get_low_t()) / 2;
      else if (num_edge_h_s > 1)
        v_split = (edge_h_s[num_edge_h_s / 2 - 1]->get_high_t() +
                   edge_h_s[num_edge_h_s / 2]->get_low_t()) / 2;
      else
      {
//        cerr << " kcurve: gen_curve_topo_proto: num_edge_l_s = " << num_edge_l_s << ", num_edge_h_s = " << num_edge_h_s << endl << flush;
        assert(num_edge_l_s >= 1 && num_edge_h_s >= 1);
        
        v_split = (edge_l_s[0]->get_low_t() + edge_l_s[0]->get_high_t() +
                   edge_h_s[0]->get_low_t() + edge_h_s[0]->get_high_t()) / 4;
      }
      
      K_POINT2D**   split_t;
      unsigned long num_split_t;
      K_POINT2D**   edge_l_s_b;
      K_POINT2D**   edge_h_s_b;
      K_POINT2D**   edge_l_s_a;
      K_POINT2D**   edge_h_s_a;
      unsigned long num_edge_l_s_b, num_edge_h_s_b;
      unsigned long num_edge_h_s_a, num_edge_l_s_a;
      K_CURVE**     curves_b;
      K_CURVE**     curves_a;
      unsigned long num_curves_b, num_curves_a, num_curves_ba, num_curves_bb;
      
      num_split_t = get_pts_proto(l_s, h_s, P.subst_val(1, v_split), v_split,
                                  P, K_RATPOLY(1, 0, v_split),
                                  split_t,
                                  init_tol, 0);
      
      if (num_split_t > 1)
        sort_s(split_t, num_split_t);
      
      for (num_edge_l_s_b = 0;
           num_edge_l_s_b < num_edge_l_s
           &&
           edge_l_s[num_edge_l_s_b]->compare_t(v_split) < 0;
           num_edge_l_s_b++)
        ;
      
      if (num_edge_l_s_b > 0)
        edge_l_s_b = new K_POINT2D* [num_edge_l_s_b];
      else  //  if (num_edge_l_s_b == 0)
        edge_l_s_b = 0;
      
      for (i = 0; i < num_edge_l_s_b; i++)
      {
        edge_l_s_b[i] = edge_l_s[i];
        edge_l_s_b[i]->ref_count++;
      }
      
      if ((num_edge_l_s_a = num_edge_l_s - num_edge_l_s_b) > 0)
        edge_l_s_a = new K_POINT2D* [num_edge_l_s_a];
      else  //  if (num_edge_l_s_a == 0)
        edge_l_s_a = 0;
      
      for (i = 0, j = num_edge_l_s_b; i < num_edge_l_s_a; i++, j++)
      {
        edge_l_s_a[i] = edge_l_s[j];
        edge_l_s_a[i]->ref_count++;
      }
      
      assert(j == num_edge_l_s);
      
      for (num_edge_h_s_b = 0;
           num_edge_h_s_b < num_edge_h_s
           &&
           edge_h_s[num_edge_h_s_b]->compare_t(v_split) < 0;
           num_edge_h_s_b++)
        ;
      
      if (num_edge_h_s_b > 0)
        edge_h_s_b = new K_POINT2D* [num_edge_h_s_b];
      else  //  if (num_edge_h_s_b == 0)
        edge_h_s_b = 0;
      
      for (i = 0; i < num_edge_h_s_b; i++)
      {
        edge_h_s_b[i] = edge_h_s[i];
        edge_h_s_b[i]->ref_count++;
      }
      
      if ((num_edge_h_s_a = num_edge_h_s - num_edge_h_s_b) > 0)
        edge_h_s_a = new K_POINT2D* [num_edge_h_s_a];
      else  //  if (num_edge_h_s_a == 0)
        edge_h_s_a = 0;
      
      for (i = 0, j = num_edge_h_s_b; i < num_edge_h_s_a; i++, j++)
      {
        edge_h_s_a[i] = edge_h_s[j];
        edge_h_s_a[i]->ref_count++;
      }
      
      assert(j == num_edge_h_s);
      
      if (l_t < v_split)
        num_curves_b = gen_curve_topo_proto(P, l_s, h_s, l_t, v_split,
                                            edge_l_s_b, num_edge_l_s_b,
                                            edge_h_s_b, num_edge_h_s_b,
                                            edge_l_t, num_edge_l_t,
                                            split_t, num_split_t,
                                            turn_s, num_turn_s,
                                            turn_t, num_turn_t,
                                            curves_b);
      else  //  if (l_t == v_split)
      {
        num_curves_b = 0;
        curves_b     = 0;
      }
      
      if (h_t > v_split)
        num_curves_a = gen_curve_topo_proto(P, l_s, h_s, v_split, h_t,
                                            edge_l_s_a, num_edge_l_s_a,
                                            edge_h_s_a, num_edge_h_s_a,
                                            split_t, num_split_t,
                                            edge_h_t, num_edge_h_t,
                                            turn_s, num_turn_s,
                                            turn_t, num_turn_t,
                                            curves_a);
      else  //  if (h_t == v_split)
      {
        num_curves_a = 0;
        curves_a     = 0;
      }
      
      //  Merge curves.
      
      num_curves_ba = num_curves_bb = 0;
      
      for (i = 0; i < num_split_t; i++)
      {
        for (j = 0;
             j < num_curves_b && !curves_b[j]->is_start_or_end(*split_t[i]);
             j++)
          ;
        
        if (j < num_curves_b)
        {
          for (k = 0;
               k < num_curves_a && !curves_a[k]->is_start_or_end(*split_t[i]);
               k++)
            ;
          
          if (k < num_curves_a)
          {
            curves_b[j]->add_on(*curves_a[k]);
            
            if (k < num_curves_a - 1)
            {
              curves_a[k] = curves_a[num_curves_a - 1];
              curves_a[k]->ref_count++;
            }
            
            num_curves_a--;
            num_curves_ba++;
          }
          else if (j < num_curves_b - 1 && !curves_b[j]->is_closed())
          {
            for (l = j + 1;
                 l < num_curves_b
                 &&
                 !curves_b[l]->is_start_or_end(*split_t[i]);
                 l++)
              ;
            
            if (l < num_curves_b)
            {
              curves_b[j]->add_on(*curves_b[l]);
              
              if (l < num_curves_b - 1)
              {
                curves_b[l] = curves_b[num_curves_b - 1];
                curves_b[l]->ref_count++;
              }
              
              num_curves_b--;
              num_curves_bb++;
            }
          }
        }
      }
      
      if ((num_curves = num_curves_b + num_curves_a) > 0)
        curves = new K_CURVE* [num_curves];
      else  //  if (num_curves == 0)
        curves = 0;
      
      for (i = 0; i < num_curves_b; i++)
      {
        curves[i] = curves_b[i];
        curves[i]->ref_count++;
      }
      
      for (i = num_curves_b, j = 0; j < num_curves_a; i++, j++)
      {
        curves[i] = curves_a[j];
        curves[i]->ref_count++;
      }
      
      assert(i == num_curves);
      
      if (num_curves_b + num_curves_bb > 0)
      {
        for (i = 0; i < num_curves_b + num_curves_bb; i++)
          if (!--curves_b[i]->ref_count)
            delete curves_b[i];
        
        delete [] curves_b;
      }
      
      if (num_curves_a + num_curves_ba > 0)
      {
        for (i = 0; i < num_curves_a + num_curves_ba; i++)
          if (!--curves_a[i]->ref_count)
            delete curves_a[i];
        
        delete [] curves_a;
      }
      
      if (num_split_t > 0)
      {
        for (i = 0; i < num_split_t; i++)
          if (!--split_t[i]->ref_count)
            delete split_t[i];
        
        delete [] split_t;
      }
      
      if (num_edge_l_s_b > 0)
      {
        for (i = 0; i < num_edge_l_s_b; i++)
          if (!--edge_l_s_b[i]->ref_count)
            delete edge_l_s_b[i];
        
        delete [] edge_l_s_b;
      }
      
      if (num_edge_l_s_a > 0)
      {
        for (i = 0; i < num_edge_l_s_a; i++)
          if (!--edge_l_s_a[i]->ref_count)
            delete edge_l_s_a[i];
        
        delete [] edge_l_s_a;
      }
      
      if (num_edge_h_s_b > 0)
      {
        for (i = 0; i < num_edge_h_s_b; i++)
          if (!--edge_h_s_b[i]->ref_count)
            delete edge_h_s_b[i];
        
        delete [] edge_h_s_b;
      }
      
      if (num_edge_h_s_a > 0)
      {
        for (i = 0; i < num_edge_h_s_a; i++)
          if (!--edge_h_s_a[i]->ref_count)
            delete edge_h_s_a[i];
        
        delete [] edge_h_s_a;
      }
    }
    else  //  if num_edge is more than 1 and odd
    {
      cerr << "   Cannot connect edge-poiunts. " << endl << flush;
      cerr << "   Degeneracy detected. " << endl << flush;
      abort();
    }
  
//  for (i = 0; i < num_curves; i++)
//    cerr << "   *curves[" << i << "] = " << endl << *curves[i] << endl << flush;
//  cerr << " kcurve: gen_curve_topo_proto: -------------------- " << endl << flush;
//  cerr << " kcurve: gen_curve_topo_proto: ==================== " << endl << flush;
  
  return num_curves;
}

//  int K_CURVE :: pt_inside_trim_curves(K_POINT2D& x,
//                                       K_CURVE** const trim_curves,
//                                       const unsigned long num_trim_curves)
//    returns 1 if x lies inside trim_curves
//            0 if x lies outside trim_curves.
//    POSSIBLY DOES NOT TERMINATE!

int pt_inside_trim_curves(K_POINT2D& x,
                          K_CURVE** const trim_curves,
                          const unsigned long num_trim_curves)
{
  unsigned long i, j;
  unsigned long c;
  bigrational   x_s, x_t, x_low_s, x_high_s, x_low_t, x_high_t;
  unsigned long num_int_pts, num_int_pts_alt;
  K_BOXCO2      b;
  K_RATPOLY     split_poly;
  K_RATPOLY     dummy;
  K_POINT2D**   split_pts;
  unsigned long num_split_pts;
  int           overlapping;
  
//  cerr << " kcurve: pt_inside_trim_curves: 0: -------------------- " << endl << flush;
  
  //  Make sure that
  //    x is contained at most one MONOTONIC segments of some trim_curves.
  
  c = 0;
//  c = 2;
  
//  cerr << " kcurve: pt_inside_trim_curves: num_trim_curves = " << num_trim_curves << endl << flush;
//  while (c > 1)
//  {
  for (c = 0, i = 0; i < num_trim_curves; i++)
  {
//    cerr << " kcurve: pt_inside_trim_curves: i = " << i << endl << flush;
    if (trim_curves[i]->contains(x))
    {
//      cerr << " kcurve: pt_inside_trim_curves: *trim_curves[" << i << "] contains x. " << endl << flush;
//      cerr << " kcurve: pt_inside_trim_curves: *trim_curves[" << i << "] = " << endl << *trim_curves[i] << endl << flush;
//      cerr << " kcurve: pt_inside)trim_curves: x = " << x << endl << flush;
      c++;
    }
//    cerr << " kcurve: pt_inside_trim_curves: c = " << c << endl << flush;
  }
//  if (c > 1)
//    x.shrink(shrink_step, shrink_step);
//  }
  
//  cerr << " kcurve: pt_inside_trim_curves: 1: ------------------------- " << endl << flush;
  
  assert(!c || c == 1);
  
  //  ray shooting
  
  overlapping = 0;
  
  if (!c)
  //  No segment of trim_curves contains x.
  //  Draw a vertical line through x and
  //    count the number of intersections of the line and trim_curves
  //      above x.
  {
    x_s         = (x.get_low_s() + x.get_high_s()) / 2;
    x_t         = (x.get_low_t() + x.get_high_t()) / 2;
    num_int_pts = 0;
    
    for (i = 0; i < num_trim_curves; i++)
    {
      b = trim_curves[i]->bbox();
      
      if (b.low[0] <= x_s && b.high[0] >= x_s)
        if ((num_split_pts =
               get_pts_proto(x_s,
                             x_t, b.high[1],
                             trim_curves[i]->poly->subst_val(0, x_s),
                             dummy, dummy,
                             split_pts,
                             init_tol, 1)) > 0)
        {
          for (j = 0; j < num_split_pts; j++)
            if (trim_curves[i]->contains(*split_pts[j]))
              num_int_pts++;
          
          for (j = 0; j < num_split_pts; j++)
            if (!--split_pts[j]->ref_count)
              delete split_pts[j];
          
          delete [] split_pts;
        }
    }
  }
  else  //  if (c == 1)
  //  There exists a MONOTONIC segment of some trim_curves which contains x.
    if (x.type == 1)
    //  For each of 4 corners of x,
    //    draw a vertical line through x and
    //      count the number of intersections of the line and trim_curves
    //        above it.
    {
      x_low_s     = x.get_low_s();
      x_high_s    = x.get_high_s();
      x_low_t     = x.get_low_t();
      x_high_t    = x.get_high_t();
      num_int_pts = num_int_pts_alt = 0;
      
      //  left side
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_low_s && b.high[0] >= x_low_s)
        {
          split_poly = trim_curves[i]->poly->subst_val(0, x_low_s);
          
          if ((num_split_pts = get_pts_proto(x_low_s,
                                             x_low_t, b.high[1],
                                             split_poly,
                                             dummy, dummy,
                                             split_pts,
                                             init_tol, 1)) > 0)
          {            
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_low_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts++;
                  
                  if (split_pts[j]->get_low_t() < x_low_t
                      &&
                      split_pts[j]->get_high_t() > x_low_t)
                    split_pts[j]->cut_t(x_low_t);
                  
                  if (split_pts[j]->get_low_t() < x_high_t
                      &&
                      split_pts[j]->get_high_t() > x_high_t)
                    split_pts[j]->cut_t(x_high_t);
                  
                  if (split_pts[j]->get_low_t() >= x_low_t
                      &&
                      split_pts[j]->get_high_t() <= x_high_t)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
      
      //  right side
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_high_s && b.high[0] >= x_high_s)
        {
          split_poly = trim_curves[i]->poly->subst_val(0, x_high_s);
          
          if ((num_split_pts = get_pts_proto(x_high_s,
                                             x_low_t, b.high[1],
                                             split_poly,
                                             dummy, dummy,
                                             split_pts,
                                             init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_low_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts_alt++;
                  
                  if (split_pts[j]->get_low_t() < x_low_t
                      &&
                      split_pts[j]->get_high_t() > x_low_t)
                    split_pts[j]->cut_t(x_low_t);
                  
                  if (split_pts[j]->get_low_t() < x_high_t
                      &&
                      split_pts[j]->get_high_t() > x_high_t)
                    split_pts[j]->cut_t(x_high_t);
                  
                  if (split_pts[j]->get_low_t() >= x_low_t
                      &&
                      split_pts[j]->get_high_t() <= x_high_t)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
      
      if (num_int_pts != num_int_pts_alt)
        overlapping = 1;
    }
    else if (x.type == 2)
    //  For each of 2 corners of x,
    //    draw a horizontal line through x and
    //      count the number of intersections of the line and trim_curves
    //        right-side of it.
    {
      x_low_s     = x.get_low_s();
      x_high_s    = x.get_high_s();
      x_t         = x.PtBt;
      num_int_pts = 0;
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[1] <= x_t && b.high[1] >= x_t)
        {
          split_poly = trim_curves[i]->poly->subst_val(1, x_t);
          
          if ((num_split_pts = get_pts_proto(x_low_s, b.high[0],
                                             split_poly,
                                             x_t,
                                             dummy, dummy,
                                             split_pts,
                                             init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBs != x_low_s)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts++;
                  
                  if (split_pts[j]->get_low_s() < x_low_s
                      &&
                      split_pts[j]->get_high_s() > x_low_s)
                    split_pts[j]->cut_s(x_low_s);
                  
                  if (split_pts[j]->get_low_s() < x_high_s
                      &&
                      split_pts[j]->get_high_s() > x_high_s)
                    split_pts[j]->cut_s(x_high_s);
                  
                  if (split_pts[j]->get_low_s() >= x_low_s
                      &&
                      split_pts[j]->get_high_s() <= x_high_s)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
    }
    else if (x.type == 3)
    //  For each of 2 corners of x,
    //    draw a vertical line through x and
    //      count the number of intersections of the line and trim_curves
    //        above it.
    {
      x_s         = x.PtBs;
      x_low_t     = x.get_low_t();
      x_high_t    = x.get_high_t();
      num_int_pts = 0;
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_s && b.high[0] >= x_s)
        {
          split_poly = trim_curves[i]->poly->subst_val(0, x_s);
          
          if ((num_split_pts = get_pts_proto(x_s,
                                             x_low_t, b.high[1],
                                             split_poly,
                                             dummy, dummy,
                                             split_pts,
                                             init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_low_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts++;
                  
                  if (split_pts[j]->get_low_t() < x_low_t
                      &&
                      split_pts[j]->get_high_t() > x_low_t)
                    split_pts[j]->cut_t(x_low_t);
                  
                  if (split_pts[j]->get_low_t() < x_high_t
                      &&
                      split_pts[j]->get_high_t() > x_high_t)
                    split_pts[j]->cut_t(x_high_t);
                  
                  if (split_pts[j]->get_low_t() >= x_low_t
                      &&
                      split_pts[j]->get_high_t() <= x_high_t)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
    }
    else  //  if (x.type == 4)
    //  Draw a vertical line through x and
    //    count the number of intersections of the line and trim_curves
    //      above it.
    {
      x_s         = x.PtBs;
      x_t         = x.PtBt;
      num_int_pts = 0;
      
      for (i = 0; i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_s && b.high[0] >= x_s)
        {
          if ((num_split_pts =
                 get_pts_proto(x_s,
                               x_t, b.high[1],
                               trim_curves[i]->poly->subst_val(0, x_s),
                               dummy, dummy,
                               split_pts,
                               init_tol, 1)) > 0)
          {
            for (j = 0; j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                  num_int_pts++;
              
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
    }
  
  if (!overlapping)
    if (!(num_int_pts % 2))
      return 0;
    else  //  if (num_int_pts % 2)
      return 1;
  else  //  if (overlapping)
  {
    x.shrink(shrink_step, shrink_step);
    
    return pt_inside_trim_curves(x, trim_curves, num_trim_curves);
  }
}

int pt_in_on_out_trim_curves(K_POINT2D& x,
                             K_CURVE** const trim_curves,
                             const unsigned long num_trim_curves)
{
  unsigned long i, j;
  unsigned long c;
  bigrational   x_s, x_t, x_low_s, x_high_s, x_low_t, x_high_t;
  unsigned long num_int_pts, num_int_pts_alt;
  K_BOXCO2      b;
  K_RATPOLY     split_poly;
  K_RATPOLY     dummy;
  K_POINT2D**   split_pts;
  unsigned long num_split_pts;
  int           overlapping;
  
  //  Make sure that
  //    x is contained at most one MONOTONIC segments of some trim_curves.
  
  c = 0;
  
  for (i = 0; i < num_trim_curves; i++)
    if (trim_curves[i]->contains(x))
      c++;
  
  assert(!c || c == 1);
  
  //  ray shooting
  
  overlapping = 0;
  
  if (!c)
  //  No segment of trim_curves contains x.
  //  Draw a vertical line through x and
  //    count the number of intersections of the line and trim_curves
  //      above x.
  {
    x_s         = (x.get_low_s() + x.get_high_s()) / 2;
    x_t         = (x.get_low_t() + x.get_high_t()) / 2;
    num_int_pts = 0;
    
    for (i = 0; i < num_trim_curves; i++)
    {
      b = trim_curves[i]->bbox();
      
      if (b.low[0] <= x_s && b.high[0] >= x_s)
        if ((num_split_pts =
               get_pts_proto(x_s,
                             x_t, b.high[1],
                             trim_curves[i]->poly->subst_val(0, x_s),
                             dummy, dummy,
                             split_pts,
                             init_tol, 1)) > 0)
        {
          for (j = 0; j < num_split_pts; j++)
            if (trim_curves[i]->contains(*split_pts[j]))
              num_int_pts++;
          
          for (j = 0; j < num_split_pts; j++)
            if (!--split_pts[j]->ref_count)
              delete split_pts[j];
          
          delete [] split_pts;
        }
    }
  }
  else  //  if (c == 1)
  //  There exists a MONOTONIC segment of some trim_curves which contains x.
    if (x.type == 1)
    //  For each of 4 corners of x,
    //    draw a vertical line through x and
    //      count the number of intersections of the line and trim_curves
    //        above it.
    {
      x_low_s     = x.get_low_s();
      x_high_s    = x.get_high_s();
      x_low_t     = x.get_low_t();
      x_high_t    = x.get_high_t();
      num_int_pts = num_int_pts_alt = 0;
      
      //  left side
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_low_s && b.high[0] >= x_low_s)
        {
          split_poly = trim_curves[i]->poly->subst_val(0, x_low_s);
          
          if (split_poly.is_zero())
            overlapping = 1;
          else if ((num_split_pts = get_pts_proto(x_low_s,
                                                  x_low_t, b.high[1],
                                                  split_poly,
                                                  dummy, dummy,
                                                  split_pts,
                                                  init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_low_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts++;
                  
                  if (split_pts[j]->get_low_t() < x_low_t
                      &&
                      split_pts[j]->get_high_t() > x_low_t)
                    split_pts[j]->cut_t(x_low_t);
                  
                  if (split_pts[j]->get_low_t() < x_high_t
                      &&
                      split_pts[j]->get_high_t() > x_high_t)
                    split_pts[j]->cut_t(x_high_t);
                  
                  if (split_pts[j]->get_low_t() >= x_low_t
                      &&
                      split_pts[j]->get_high_t() <= x_high_t)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
      
      //  right side
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_high_s && b.high[0] >= x_high_s)
        {
          split_poly = trim_curves[i]->poly->subst_val(0, x_high_s);
          
          if (split_poly.is_zero())
            overlapping = 1;
          else if ((num_split_pts = get_pts_proto(x_high_s,
                                                  x_low_t, b.high[1],
                                                  split_poly,
                                                  dummy, dummy,
                                                  split_pts,
                                                  init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_low_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts_alt++;
                  
                  if (split_pts[j]->get_low_t() < x_low_t
                      &&
                      split_pts[j]->get_high_t() > x_low_t)
                    split_pts[j]->cut_t(x_low_t);
                  
                  if (split_pts[j]->get_low_t() < x_high_t
                      &&
                      split_pts[j]->get_high_t() > x_high_t)
                    split_pts[j]->cut_t(x_high_t);
                  
                  if (split_pts[j]->get_low_t() >= x_low_t
                      &&
                      split_pts[j]->get_high_t() <= x_high_t)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
      
      if (num_int_pts != num_int_pts_alt)
        overlapping = 1;
    }
    else if (x.type == 2)
    //  For each of 2 corners of x,
    //    draw a horizontal line through x and
    //      count the number of intersections of the line and trim_curves
    //        right-side of it.
    {
      x_low_s     = x.get_low_s();
      x_high_s    = x.get_high_s();
      x_t         = x.PtBt;
      num_int_pts = 0;
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[1] <= x_t && b.high[1] >= x_t)
        {
          split_poly = trim_curves[i]->poly->subst_val(1, x_t);
          
          if (split_poly.is_zero())
            overlapping = 1;
          else if ((num_split_pts = get_pts_proto(x_low_s, b.high[0],
                                                  split_poly,
                                                  x_t,
                                                  dummy, dummy,
                                                  split_pts,
                                                  init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBs != x_low_s)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts++;
                  
                  if (split_pts[j]->get_low_s() < x_low_s
                      &&
                      split_pts[j]->get_high_s() > x_low_s)
                    split_pts[j]->cut_s(x_low_s);
                  
                  if (split_pts[j]->get_low_s() < x_high_s
                      &&
                      split_pts[j]->get_high_s() > x_high_s)
                    split_pts[j]->cut_s(x_high_s);
                  
                  if (split_pts[j]->get_low_s() >= x_low_s
                      &&
                      split_pts[j]->get_high_s() <= x_high_s)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
    }
    else if (x.type == 3)
    //  For each of 2 corners of x,
    //    draw a vertical line through x and
    //      count the number of intersections of the line and trim_curves
    //        above it.
    {
      x_s         = x.PtBs;
      x_low_t     = x.get_low_t();
      x_high_t    = x.get_high_t();
      num_int_pts = 0;
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_s && b.high[0] >= x_s)
        {
          split_poly = trim_curves[i]->poly->subst_val(0, x_s);
          
          if (split_poly.is_zero())
            overlapping = 1;
          else if ((num_split_pts = get_pts_proto(x_s,
                                                  x_low_t, b.high[1],
                                                  split_poly,
                                                  dummy, dummy,
                                                  split_pts,
                                                  init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_low_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                {
                  num_int_pts++;
                  
                  if (split_pts[j]->get_low_t() < x_low_t
                      &&
                      split_pts[j]->get_high_t() > x_low_t)
                    split_pts[j]->cut_t(x_low_t);
                  
                  if (split_pts[j]->get_low_t() < x_high_t
                      &&
                      split_pts[j]->get_high_t() > x_high_t)
                    split_pts[j]->cut_t(x_high_t);
                  
                  if (split_pts[j]->get_low_t() >= x_low_t
                      &&
                      split_pts[j]->get_high_t() <= x_high_t)
                    overlapping = 1;
                }
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
        }
      }
    }
    else  //  if (x.type == 4)
    //  Draw a vertical line through x and
    //    count the number of intersections of the line and trim_curves
    //      above it.
    {
      x_s         = x.PtBs;
      x_t         = x.PtBt;
      num_int_pts = 0;
      
      for (i = 0; !overlapping && i < num_trim_curves; i++)
      {
        b = trim_curves[i]->bbox();
        
        if (b.low[0] <= x_s && b.high[0] >= x_s)
          if ((num_split_pts =
               get_pts_proto(x_s,
                             x_t, b.high[1],
                             trim_curves[i]->poly->subst_val(0, x_s),
                             dummy, dummy,
                             split_pts,
                             init_tol, 1)) > 0)
          {
            for (j = 0; !overlapping && j < num_split_pts; j++)
              if (split_pts[j]->type != 4 || split_pts[j]->PtBt != x_t)
              {
                if (trim_curves[i]->contains(*split_pts[j]))
                  num_int_pts++;
              }
              else  //  if (split_pts[j]->type == 4
                    //      &&
                    //      split_pts[j]->PtBt != x_t)
                if (trim_curves[i]->contains(*split_pts[j]))
                  overlapping = 1;
            
            for (j = 0; j < num_split_pts; j++)
              if (!--split_pts[j]->ref_count)
                delete split_pts[j];
            
            delete [] split_pts;
          }
      }
    }
  
  if (!overlapping)
    if (!(num_int_pts % 2))
      return - 1;
    else  //  if (num_int_pts % 2)
      return 1;
  else  //  if (overlapping)
    return 0;
}

int K_CURVE :: mk_seg_fw(int* const marks, const int good_bad,
                         const long top_seg,
                         K_POINT2D** const tail_pts1,
                         const unsigned long num_tail_pts1,
                         K_POINT2D** const tail_pts2,
                         const unsigned long num_tail_pts2) const
{
  assert(!good_bad || good_bad == 1);
  
  long          i;
  unsigned long j, k;
  int           c;
  
  i = top_seg;
  
  if (i == num_segments && is_closed())
    i = 0;
  else if (i == - 1 && is_closed())
    i = num_segments - 1;
  
  assert(i >= 0 && i < num_segments);
  c = 1;
  
  while (c)
  {
    marks[i] = good_bad;
    i++;
    
    if (i >= num_segments && is_closed())
      i = 0;
    
    if (i >= 0 && i < num_segments)
    {
      for (j = 0;
           j < num_tail_pts1 && !segments[i]->start->equiv(*tail_pts1[j]);
           j++)
        ;
      
      if (j < num_tail_pts1)
        c = 0;
      else  //  if (j == num_tail_pts1)
      {
        for (k = 0;
             k < num_tail_pts2 && !segments[i]->start->equiv(*tail_pts2[k]);
             k++)
          ;
        
        if (k < num_tail_pts2)
          c = 0;
      }
    }
    else
      c = 0;
  }
  
  return 0;
}

int K_CURVE :: mk_seg_bw(int* const marks, const int good_bad,
                         const long top_seg,
                         K_POINT2D** const tail_pts1,
                         const unsigned long num_tail_pts1,
                         K_POINT2D** const tail_pts2,
                         const unsigned long num_tail_pts2) const
{
  assert(!good_bad || good_bad == 1);
  
  long          i;
  unsigned long j, k;
  int           c;
  
  i = top_seg;
  
  if (i == num_segments && is_closed())
    i = 0;
  
  if (i == - 1 && is_closed())
    i = num_segments - 1;
  
  assert(i >= 0 && i < num_segments);
  c = 1;
  
  while (c)
  {
    marks[i] = good_bad;
    i--;
    
    if (i < 0 && is_closed())
      i = num_segments - 1;
    
    if (i >= 0 && i < num_segments)
    {
      for (j = 0;
           j < num_tail_pts1 && !segments[i]->end->equiv(*tail_pts1[j]);
           j++)
        ;
      
      if (j < num_tail_pts1)
        c = 0;
      else
      {
        for (k = 0;
             k < num_tail_pts2 && !segments[i]->end->equiv(*tail_pts2[k]);
             k++)
          ;
        
        if (k < num_tail_pts2)
          c = 0;
      }
    }
    else
      c = 0;
  }
  
  return 0;
}

int K_CURVE :: subdivide(const unsigned long num_cut)
{
  assert(num_cut > 0);
  
  if (num_cut > 1)
  {
    unsigned long i, j;
    K_BOXCO2      b;
    bigrational   b_width[2];
    unsigned long cut_dir;
    bigrational   val_step, val_cut;
    K_POINT2D**   cut_pts;
    unsigned long num_cut_pts;
    
    b = bbox();
    
    for (i = 0; i < 2; i++)
      b_width[i] = b.high[i] - b.low[i];
    
    if (b_width[0] < b_width[1])
      cut_dir = 1;
    else  //  if (b_width[0] >= b_width[1])
      cut_dir = 0;
    
    if ((val_step  = b_width[cut_dir] / (num_cut + 1)) > 0)
    {
      val_cut = b.low[cut_dir] + val_step;
      
      for (i = 0; i < num_cut; i++)
      {
        if ((num_cut_pts = find_intersections(K_RATPOLY(2, cut_dir, val_cut),
                                              cut_pts,
                                              1)) > 0)
        {
          for (j = 0; j < num_cut_pts; j++)
            add_pt(cut_pts[j]);
          
          for (j = 0; j < num_cut_pts; j++)
            if (!--cut_pts[j]->ref_count)
              delete cut_pts[j];
          
          delete [] cut_pts;
        }
        
        val_cut += val_step;
      }
    }
  }
  
  return 0;
}

