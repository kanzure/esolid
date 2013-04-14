//  file:    kpatch.cc
//  update:  11/18/02

#include <config.h>

#include <kpatch.h>
#include <pt_surf_assoc.h>

static unsigned long PATCH_ID_COUNT = 0;
static PT_SURF_ASSOC pt_surf_assoc;

K_PATCH :: K_PATCH()
  : ID(PATCH_ID_COUNT++),
    surf(0), low_s(0), high_s(0), low_t(0), high_t(0),
    num_trim_curves(0), trim_curves(0), adj_surfs(0), adj_patches(0),
    num_int_curves(0), int_curves(0), adj_int_surfs(0), adj_int_patches(0),
    num_merged(0), ic(0), len_ic(0), ic_closed(0),
    ref_count(0)
{ }

K_PATCH :: K_PATCH(K_SURF* const s,
                   K_CURVE* const t[], const unsigned long n)
  : ID(PATCH_ID_COUNT++)
{
  unsigned long i, j;
  
  surf = s;
  surf->ref_count++;
  
  if ((num_trim_curves = n) > 0)
  {
    trim_curves = new K_CURVE* [num_trim_curves];
    adj_surfs   = new K_SURF* [num_trim_curves];
    adj_patches = new K_PATCH* [num_trim_curves];
    
    for (i = 0; i < num_trim_curves; i++)
    {
      trim_curves[i] = t[i];
      trim_curves[i]->ref_count++;
      adj_surfs[i]   = 0;
      adj_patches[i] = 0;
    }
  }
  else  //  if (!num_trim_curves)
  {
    trim_curves = 0;
    adj_surfs   = 0;
    adj_patches = 0;
  }
  
  num_int_curves  = 0;
  int_curves      = 0;
  adj_int_surfs   = 0;
  adj_int_patches = 0;
  
  num_merged = 0;
  ic        = 0;
  len_ic    = 0;
  ic_closed = 0;
  
  set_range();
  
  ref_count = 0;
}

K_PATCH :: K_PATCH(const K_PARTITION& p)
  : ID(PATCH_ID_COUNT++)
{
  unsigned long i;
  
  surf = p.from->surf;
  surf->ref_count++;
  
  if ((num_trim_curves = p.num_trim_curves) > 0)
  {
    trim_curves = new K_CURVE* [num_trim_curves];
    adj_surfs   = new K_SURF* [num_trim_curves];
    adj_patches = new K_PATCH* [num_trim_curves];
    
    for (i = 0; i < num_trim_curves; i++)
    {
      trim_curves[i] = p.trim_curves[i];
      trim_curves[i]->ref_count++;
      trim_curves[i]->assoc(p.trim_curves[i]->curve_in_other_dom,
                            p.trim_curves[i]->dir_in_other_dom);
      
      if (p.rev_curves[i])
        trim_curves[i]->reverse();
      
      adj_surfs[i]   = 0;
      adj_patches[i] = 0;
    }
  }
  else  //  if (!num_trim_curves)
  {
    trim_curves = 0;
    adj_surfs   = 0;
    adj_patches = 0;
  }
  
  num_int_curves  = 0;
  int_curves      = 0;
  adj_int_surfs   = 0;
  adj_int_patches = 0;
  
  num_merged = 0;
  ic        = 0;
  len_ic    = 0;
  ic_closed = 0;
  
  set_range();
  
  ref_count = 0;
}

K_PATCH :: K_PATCH(const K_PATCH& p)
  : ID(PATCH_ID_COUNT++)
{
  unsigned long i, j;
  
  surf   = p.surf;
  surf->ref_count++;
  low_s  = p.low_s;
  high_s = p.high_s;
  low_t  = p.low_t;
  high_t = p.high_t;
  
  if ((num_trim_curves = p.num_trim_curves) > 0)
  {
    trim_curves = new K_CURVE* [num_trim_curves];
    adj_surfs   = new K_SURF* [num_trim_curves];
    adj_patches = new K_PATCH* [num_trim_curves];
    
    for (i = 0; i < num_trim_curves; i++)
    {
      trim_curves[i] = p.trim_curves[i];
      trim_curves[i]->ref_count++;
      
      if (adj_surfs[i] = p.adj_surfs[i])
        adj_surfs[i]->ref_count++;
      
      if (adj_patches[i] = p.adj_patches[i])
        adj_patches[i]->ref_count++;
    }
  }
  else  //  if (!num_trim_curves)
  {
    trim_curves = 0;
    adj_surfs   = 0;
    adj_patches = 0;
  }
  
  if ((num_int_curves = p.num_int_curves) > 0)
  {
    int_curves      = new K_CURVE* [num_int_curves];
    adj_int_surfs   = new K_SURF* [num_int_curves];
    adj_int_patches = new K_PATCH* [num_int_curves];
    
    for (i = 0; i < num_int_curves; i++)
    {
      int_curves[i] = p.int_curves[i];
      int_curves[i]->ref_count++;
      
      if (adj_int_surfs[i] = p.adj_int_surfs[i])
        adj_int_surfs[i]->ref_count++;
      
      if (adj_int_patches[i] = p.adj_int_patches[i])
        adj_int_patches[i]->ref_count++;
    }
  }
  else  //  if (!num_int_curves)
  {
    int_curves      = 0;
    adj_int_surfs   = 0;
    adj_int_patches = 0;
  }
  
  if ((num_merged = p.num_merged) > 0)
  {
    ic        = new K_CURVE** [num_merged];
    len_ic    = new unsigned long [num_merged];
    ic_closed = new int [num_merged];
    
    for (i = 0; i < num_merged; i++)
    {
      if ((len_ic[i] = p.len_ic[i]) > 0)
        ic[i] = new K_CURVE* [len_ic[i]];
      else
        ic[i] = 0;
      
      ic_closed[i] = p.ic_closed[i];
      
      for (j = 0; j < len_ic[i]; j++)
      {
        ic[i][j] = p.ic[i][j];
        ic[i][j]->ref_count++;
      }
    }
  }
  else  //  if (!num_merged)
  {
    ic        = 0;
    len_ic    = 0;
    ic_closed = 0;
  }
  
  ref_count = 0;
}

K_PATCH& K_PATCH :: operator =(const K_PATCH& p)
{
  if (this != &p)
  {
    unsigned long i, j;
    
    if (surf && !--surf->ref_count)
      delete surf;
    
    if (num_trim_curves > 0)
    {
      for (i = 0; i < num_trim_curves; i++)
      {
        if (!--trim_curves[i]->ref_count)
          delete trim_curves[i];
        
        if (adj_surfs[i] && !--adj_surfs[i]->ref_count)
          delete adj_surfs[i];
        
        if (adj_patches[i] && !--adj_patches[i]->ref_count)
          delete adj_patches[i];
      }
      
      delete [] trim_curves;
      delete [] adj_surfs;
      delete [] adj_patches;
    }
    
    if (num_int_curves > 0)
    {
      for (i = 0; i < num_int_curves; i++)
      {
        if (!--int_curves[i]->ref_count)
          delete int_curves[i];
        
        if (adj_int_surfs[i] && !--adj_int_surfs[i]->ref_count)
          delete adj_int_surfs[i];
        
        if (adj_int_patches[i] && !--adj_int_patches[i]->ref_count)
          delete adj_int_patches[i];
      }
      
      delete [] int_curves;
      delete [] adj_int_surfs;
      delete [] adj_int_patches;
    }
    
    if (num_merged > 0)
    {
      for (i = 0; i < num_merged; i++)
      {
        for (j = 0; j < len_ic[i]; j++)
          if (!--ic[i][j]->ref_count)
            delete ic[i][j];
        
        if (len_ic[i] > 0)
          delete [] ic[i];
      }
      
      delete [] ic;
      delete [] len_ic;
      delete [] ic_closed;
    }
    
    surf   = p.surf;
    surf->ref_count++;
    low_s  = p.low_s;
    high_s = p.high_s;
    low_t  = p.low_t;
    high_t = p.high_t;
    
    if ((num_trim_curves = p.num_trim_curves) > 0)
    {
      trim_curves = new K_CURVE* [num_trim_curves];
      adj_surfs   = new K_SURF* [num_trim_curves];
      adj_patches = new K_PATCH* [num_trim_curves];
      
      for (i = 0; i < num_trim_curves; i++)
      {
        trim_curves[i] = p.trim_curves[i];
        trim_curves[i]->ref_count++;
        
        if (adj_surfs[i] = p.adj_surfs[i])
          adj_surfs[i]->ref_count++;
        
        if (adj_patches[i] = p.adj_patches[i])
          adj_patches[i]->ref_count++;
      }
    }
    else  //  if (!num_trim_curves)
    {
      trim_curves = 0;
      adj_surfs   = 0;
      adj_patches = 0;
    }
    
    if ((num_int_curves = p.num_int_curves) > 0)
    {
      int_curves      = new K_CURVE* [num_int_curves];
      adj_int_surfs   = new K_SURF* [num_int_curves];
      adj_int_patches = new K_PATCH* [num_int_curves];
      
      for (i = 0; i < num_int_curves; i++)
      {
        int_curves[i] = p.int_curves[i];
        int_curves[i]->ref_count++;
        
        if (adj_int_surfs[i] = p.adj_int_surfs[i])
          adj_int_surfs[i]->ref_count++;
        
        if (adj_int_patches[i] = p.adj_int_patches[i])
          adj_int_patches[i]->ref_count++;
      }
    }
    else  //  if (!num_int_curves)
    {
      int_curves      = 0;
      adj_int_surfs   = 0;
      adj_int_patches = 0;
    }
    
    if ((num_merged = p.num_merged) > 0)
    {
      ic        = new K_CURVE** [num_merged];
      len_ic    = new unsigned long [num_merged];
      ic_closed = new int [num_merged];
      
      for (i = 0; i < num_merged; i++)
      {
        if (len_ic[i] = p.len_ic[i])
          ic[i] = new K_CURVE* [len_ic[i]];
        else
          ic[i] = 0;
        
        ic_closed[i] = p.ic_closed[i];
        
        for (j = 0; j < len_ic[i]; j++)
        {
          ic[i][j] = p.ic[i][j];
          ic[i][j]->ref_count++;
        }
      }
    }
    else  //  if (!num_merged)
    {
      ic        = 0;
      len_ic    = 0;
      ic_closed = 0;
    }
  }
  
  return *this;
}

K_PATCH :: ~K_PATCH()
{
  unsigned long i, j;
  
  if (surf && !--surf->ref_count)
    delete surf;
  
  if (num_trim_curves > 0)
  {
    for (i = 0; i < num_trim_curves; i++)
    {
      if (!--trim_curves[i]->ref_count)
        delete trim_curves[i];
      
      if (adj_surfs[i] && !--adj_surfs[i]->ref_count)
        delete adj_surfs[i];
      
      if (adj_patches[i] && !--adj_patches[i]->ref_count)
        delete adj_patches[i];
    }
    
    delete [] trim_curves;
    delete [] adj_surfs;
    delete [] adj_patches;
  }
  
  if (num_int_curves > 0)
  {
    for (i = 0; i < num_int_curves; i++)
    {
      if (!--int_curves[i]->ref_count)
        delete int_curves[i];
      
      if (adj_int_surfs[i] && !--adj_int_surfs[i]->ref_count)
        delete adj_int_surfs[i];
      
      if (adj_int_patches[i] && !--adj_int_patches[i]->ref_count)
        delete adj_int_patches[i];
    }
    
    delete [] int_curves;
    delete [] adj_int_surfs;
    delete [] adj_int_patches;
  }
  
  if (num_merged > 0)
  {
    for (i = 0; i < num_merged; i++)
    {
      for (j = 0; j < len_ic[i]; j++)
        if (!--ic[i][j]->ref_count)
          delete ic[i][j];
      
      if (len_ic[i] > 0)
        delete [] ic[i];
    }
    
    delete [] ic;
    delete [] len_ic;
    delete [] ic_closed;
  }
}

ostream& K_PATCH :: output(ostream& o) const
{
  unsigned long i;
  
  o << " patch ID: " << ID << flush;
//  o << endl << "   [ (" << low_s << ", " << low_t;
//  o << "), (" << high_s << ", " << high_t << ") ]" << flush;
//  
//  if (surf)
//    o << endl << "   surf: " << *surf << flush;
//  
  if (num_trim_curves > 0)
    for (i = 0; i < num_trim_curves; i++)
    {
      o << "   trim_curves[" << i << "] = " << endl << flush;
      o << *trim_curves[i] << endl << flush;
    }
  
  return o;
}

ostream& operator <<(ostream& o, const K_PATCH& p)
{
  return p.output(o);
}

int K_PATCH :: set_range()
{
  unsigned long i;
  K_BOXCO2      b;
  bigrational   diff_s, diff_t;
  
  if (num_trim_curves > 0)
  {
    b      = trim_curves[0]->bbox();
    low_s  = b.low[0];
    high_s = b.high[0];
    low_t  = b.low[1];
    high_t = b.high[1];
    
    for (i = 1; i < num_trim_curves; i++)
    {
      b = trim_curves[i]->bbox();
      
      if (b.low[0] < low_s)
        low_s = b.low[0];
      
      if (b.high[0] > high_s)
        high_s = b.high[0];
      
      if (b.low[1] < low_t)
        low_t = b.low[1];
      
      if (b.high[1] > high_t)
        high_t = b.high[1];
    }
  }
  
  assert(sgn(diff_s = high_s - low_s) > 0);
  assert(sgn(diff_t = high_t - low_t) > 0);
  
  diff_s /= 10;
  diff_t /= 10;
  low_s  -= diff_s;
  high_s += diff_s;
  low_t  -= diff_t;
  high_t += diff_t;
  
  return 0;
}

bigrational K_PATCH :: get_low_s() const
{
  return low_s;
}

bigrational K_PATCH :: get_high_s() const
{
  return high_s;
}

bigrational K_PATCH :: get_low_t() const
{
  return low_t;
}

bigrational K_PATCH :: get_high_t() const
{
  return high_t;
}

int K_PATCH :: in_dom(K_POINT2D& x)
{
  assert(x.type > 0);
  
  bigrational x_l_s, x_h_s, x_l_t, x_h_t;
  int         c;
  
  x.cut_s(low_s);
  x.cut_s(high_s);
  x.cut_t(low_t);
  x.cut_t(high_t);
  
  x_l_s = x.get_low_s();
  x_h_s = x.get_high_s();
  x_l_t = x.get_low_t();
  x_h_t = x.get_high_t();
  
  if (x.type == 1)
    if (x_l_s >= low_s && x_h_s <= high_s && x_l_t >= low_t && x_h_t <= high_t)
      c = 1;
    else
      c = 0;
  else if (x.type == 2)
    if (x_l_s >= low_s && x_h_s <= high_s && x_l_t > low_t && x_h_t < high_t)
      c = 1;
    else
      c = 0;
  else if (x.type == 3)
    if (x_l_s > low_s && x_h_s < high_s && x_l_t >= low_t && x_h_t <= high_t)
      c = 1;
    else
      c = 0;
  else  //  if (x.type == 4)
    if (x_l_s > low_s && x_h_s < high_s && x_l_t > low_t && x_h_t < high_t)
      c = 1;
    else
      c = 0;
  
  return c;
}

int K_PATCH :: contains(K_POINT2D& x)
{
  int c;
  
//  cerr << " kpatch: contains: -------------------- " << endl << flush;
  
  if (in_dom(x))
  {
//    cerr << " kpatch: contains: in_dom " << endl << flush;
    c = pt_inside_trim_curves(x, trim_curves, num_trim_curves);
  }
  else  //  if (!in_dom(x))
  {
//    cerr << " kpatch: contains: not in_dom " << endl << flush;
    c = 0;
  }
  
//  cerr << " kpatch: contains: c = " << c << endl << flush;
//  cerr << " kpatch: contains: ==================== " << endl << flush;
  
  return c;
}

int K_PATCH :: in_on_out(K_POINT2D& x)
{
  if (in_dom(x))
    return pt_in_on_out_trim_curves(x, trim_curves, num_trim_curves);
  else
    return - 1;
}

#define MAX_NUM_NEW_INT_PTS 64

int K_PATCH :: intersect(K_PATCH& p)
{
  unsigned long i, j, k, l;
  
  //  1. Compute intersection curves;
  //       int_curve1 on this domain and int_curve2 on p's domain.
  
  K_RATPOLY   int_curve1, int_curve2;
  bigrational v_l[2];
  bigrational v_h[2];
  bigrational v_min, v_max;
  
//  cerr << " Debug0: kpatch: intersect: ---------- " << endl << flush;
//  cerr << "   surf->Impl = " << endl << *surf->Impl << endl << flush;
//  cerr << "   p.surf->Impl = " << endl << *p.surf->Impl << endl << flush;
//  cerr << " ------------------------------------ " << endl << flush;
//  cerr << "   surf->X = " << endl << *surf->X << endl << flush;
//  cerr << "   surf->Y = " << endl << *surf->Y << endl << flush;
//  cerr << "   surf->Z = " << endl << *surf->Z << endl << flush;
//  cerr << "   surf->W = " << endl << *surf->W << endl << flush;
//  cerr << "   p.surf->X = " << endl << *p.surf->X << endl << flush;
//  cerr << "   p.surf->Y = " << endl << *p.surf->Y << endl << flush;
//  cerr << "   p.surf->Z = " << endl << *p.surf->Z << endl << flush;
//  cerr << "   p.surf->W = " << endl << *p.surf->W << endl << flush;
//  cerr << " ------------------------------------- " << endl << flush;
  
  int_curve1 =
    p.surf->Impl->subst_param_expr(*surf->X, *surf->Y, *surf->Z, *surf->W);
//  cerr << "   int_curve1 = " << endl << int_curve1 << endl << flush;
  
  if (!int_curve1.deg[0] && !int_curve1.deg[1] && !sgn(int_curve1.coeffs[0]))
  //  if int_curve1 is 0
  {
    cerr << " Error: kpatch: intersect: Patches overlap and coplanar. " << endl << flush;
    abort();
  }
  
  v_l[0] = low_s;
  v_l[1] = low_t;
  v_h[0] = high_s;
  v_h[1] = high_t;
  int_curve1.eval_range(v_l, v_h, v_min, v_max);
  
//  cerr << "   int_curve1 at [(" << v_l[0] << ", " << v_l[1] << ") x (" << v_h[0] << ", " << v_h[1] << ")] ranges (" << v_min << ", " << v_max << ") " << endl << flush;
  
  if (sgn(v_min) > 0 || sgn(v_max) < 0)
  {
    cerr << " Error: kpatch: intersect: Patches do not intersect. " << endl << flush;
    return 0;
  }
  
  int_curve2 =
    surf->Impl->subst_param_expr(*p.surf->X, *p.surf->Y, *p.surf->Z, *p.surf->W);
//  cerr << "   int_curve2 = " << endl << int_curve2 << endl << flush;
  
  if (!int_curve2.deg[0] && !int_curve2.deg[1] && !sgn(int_curve2.coeffs[0]))
  //  if int_curve2 is 0
  {
    cerr << " Error: kpatch: intersect: Patches overlap and coplanar. " << endl << flush;
    abort();
  }
  
  v_l[0] = p.low_s;
  v_l[1] = p.low_t;
  v_h[0] = p.high_s;
  v_h[1] = p.high_t;
  int_curve2.eval_range(v_l, v_h, v_min, v_max);
  
//  cerr << "   int_curve2 at [(" << v_l[0] << ", " << v_l[1] << ") x (" << v_h[0] << ", " << v_h[1] << ")] ranges (" << v_min << ", " << v_max << ") " << endl << flush;
  
  if (sgn(v_min) > 0 || sgn(v_max) < 0)
  {
    cerr << " Error: kpatch: intersect: Patches do not intersect. " << endl << flush;
    return 0;
  }
  
//  cerr << " ------------------------------------- " << endl << flush;
  
  //  2. Resolve curve topologies.
  
  K_CURVE**     int_curves1;
  K_CURVE**     int_curves2;
  unsigned long num_int_curves1, num_int_curves2;
  
//  cerr << " Debug1: kpatch: intersect: ---------- " << endl << flush;
  
  num_int_curves1 = gen_curve_topo(int_curve1,
                                   low_s, high_s, low_t, high_t,
                                   int_curves1);
//  cerr << "   num_int_curves1 = " << num_int_curves1 << endl << flush;
  
//  if (!(num_int_curves1 = gen_curve_topo(int_curve1,
//                                         low_s, high_s, low_t, high_t,
//                                         int_curves1)))
  if (!num_int_curves1)
  {
    cerr << " Error: kpatch: intersect: num_inter_curves1 = 0 " << endl << flush;
    return 0;
  }
  
  num_int_curves2 = gen_curve_topo(int_curve2,
                                   p.low_s, p.high_s, p.low_t, p.high_t,
                                   int_curves2);
//  cerr << "   num_int_curves2 = " << num_int_curves2 << endl << flush;
  
//  if (!(num_int_curves2 = gen_curve_topo(int_curve2,
//                                         p.low_s, p.high_s, p.low_t, p.high_t,
//                                         int_curves2)))
  if (!num_int_curves2)
  {
    cerr << " Error: kpatch: intersect: num_inter_curves2 = 0 " << endl << flush;
    return 0;
  }
  
//  for (i = 0; i < num_int_curves1; i++)
//  {
//    cerr << "   *int_curves1[" << i << "]->poly = " << endl << *int_curves1[i]->poly << endl << flush;
//    cerr << "   *int_curves1[" << i << "] = " << endl << *int_curves1[i] << endl << flush;
//  }
//  cerr << " ------------------------------------- " << endl << flush;
//  for (i = 0; i < num_int_curves2; i++)
//  {
//    cerr << "   *int_curves2[" << i << "]->poly = " << endl << *int_curves2[i]->poly << endl << flush;
//    cerr << "   *int_curves2[" << i << "] = " << endl << *int_curves2[i] << endl << flush;
//  }
//  cerr << " ------------------------------------- " << endl << flush;
  
  //  3. Clip curves to trim_curves and invert points.
  
  K_POINT2D**   int_pts1;
  K_POINT2D**   int_pts2;
  unsigned long num_int_pts1, num_int_pts2;
  
  unsigned long* num_new_int_pts1;
  unsigned long* num_new_trim_pts1;
  unsigned long* num_new_int_pts2;
  unsigned long* num_new_trim_pts2;
  
  K_POINT2D**    new_int_pts1_this;
  K_POINT2D**    new_int_pts1_other;
  unsigned long* ind_int_curves1_this;
  unsigned long* ind_int_curves1_other;
  unsigned long  num_new_int_pts1_this, num_new_int_pts1_other;
  K_POINT2D**    new_int_pts2_this;
  K_POINT2D**    new_int_pts2_other;
  unsigned long* ind_int_curves2_this;
  unsigned long* ind_int_curves2_other;
  unsigned long  num_new_int_pts2_this, num_new_int_pts2_other;
  
  num_new_int_pts1  = new unsigned long [num_int_curves1];
  num_new_trim_pts1 = new unsigned long [num_int_curves1];
  
  for (i = 0; i < num_int_curves1; i++)
    num_new_int_pts1[i] = num_new_trim_pts1[i] = 0;
  
  num_new_int_pts2  = new unsigned long [num_int_curves2];
  num_new_trim_pts2 = new unsigned long [num_int_curves2];
  
  for (i = 0; i < num_int_curves2; i++)
    num_new_int_pts2[i] = num_new_trim_pts2[i] = 0;
  
  new_int_pts1_this     = new K_POINT2D* [MAX_NUM_NEW_INT_PTS];
  new_int_pts1_other    = new K_POINT2D* [MAX_NUM_NEW_INT_PTS];
  ind_int_curves1_this  = new unsigned long [MAX_NUM_NEW_INT_PTS];
  ind_int_curves1_other = new unsigned long [MAX_NUM_NEW_INT_PTS];
  
  new_int_pts2_this     = new K_POINT2D* [MAX_NUM_NEW_INT_PTS];
  new_int_pts2_other    = new K_POINT2D* [MAX_NUM_NEW_INT_PTS];
  ind_int_curves2_this  = new unsigned long [MAX_NUM_NEW_INT_PTS];
  ind_int_curves2_other = new unsigned long [MAX_NUM_NEW_INT_PTS];
  
  num_new_int_pts1_this = num_new_int_pts1_other = 0;
  num_new_int_pts2_this = num_new_int_pts2_other = 0;
  
//  cerr << " Debug2: kpatch: intersect: ---------- " << endl << flush;
  
  for (i = 0; i < num_trim_curves; i++)
  {
    //  Find intersections of trim_curves and int_curve1 on this domain.
    
    if (trim_curves[(i + 1) % num_trim_curves]->poly->
        eq_upto_const(int_curve1))
    //  The next trim_curve is identical to int_curve1.
    //  Don't count an intersection at the end.
      num_int_pts1 =
        trim_curves[i]->find_intersections(int_curve1, int_pts1, 0);
    else
    //  The next trim_curve is different from int_curve1.
    //  Count an intersection at the end.
      num_int_pts1 =
        trim_curves[i]->find_intersections(int_curve1, int_pts1, 1);
    
//    cerr << endl << "   int_curve1 = " << endl << int_curve1 << endl << flush;
//    cerr << "   *trim_curves[" << i << "]->poly = " << endl << *trim_curves[i]->poly << flush;
//    cerr << "   *trim_curves[" << i << "] = " << *trim_curves[i] << endl << flush;
//    cerr << " ------------------------------------- " << endl << flush;
//    cerr << "   int_curve1 has " << num_int_pts1 << " intersections with trim_curves[" << i << "]. " << endl << flush;
//    cerr << " ------------------------------------- " << endl << flush;
    
    if (num_int_pts1 > 0)
    {
      //  Find intersections on p's domain.
      
//      num_int_pts2 =
//        get_all_pts(int_curve2,
//                    adj_surfs[i]->Impl->subst_param_expr(*p.surf->X,
//                                                         *p.surf->Y,
//                                                         *p.surf->Z,
//                                                         *p.surf->W),
//                    int_pts2,
//                    0);
      num_int_pts2 =
        get_all_int_pts(int_curve2, *adj_surfs[i], *p.surf, int_pts2);
//      cerr << "   num_int_pts1 = " << num_int_pts1 << ", num_int_pts2 = " << num_int_pts2 << endl << flush;
      assert(num_int_pts1 <= num_int_pts2);
      
      //  Match intersections on this domain and p's domain.
      
      match_pts(int_pts1, num_int_pts1, surf, int_pts2, num_int_pts2, p.surf);
      
//      for (j = 0; j < num_int_pts1; j++)
//        cerr << "   int_pts1[" << j << "] = " << *int_pts1[j] << ", int_pts2[" << j << "] = " << *int_pts2[j] << endl << flush;
//      cerr << " ------------------------------------- " << endl << flush;
      
      //  Add intersections to trim_curves, int_curves1 and int_curves2.
      
      for (j = 0; j < num_int_pts1; j++)
      {
        //  Associate intersections on this domain and p's domain.
        
        pt_surf_assoc.record_as_assoc(int_pts1[j], surf, int_pts2[j], p.surf);
        
        //  Add intersections on this domain to trim_curves.
        
        assert(trim_curves[i]->add_pt(int_pts1[j]));
        
        //  Add intersections on this domain to int_curves1.
        
        for (k = 0;
             k < num_int_curves1 && !int_curves1[k]->add_pt(int_pts1[j]);
             k++)
          ;
        
        assert(k < num_int_curves1);
        num_new_int_pts1[k]++;
        num_new_trim_pts1[k]++;
        assert(num_new_int_pts1_this < MAX_NUM_NEW_INT_PTS);
        new_int_pts1_this[num_new_int_pts1_this]    = int_pts1[j];
        new_int_pts1_this[num_new_int_pts1_this]->ref_count++;
        ind_int_curves1_this[num_new_int_pts1_this] = k;
        num_new_int_pts1_this++;
        
        //  Add intersections on p's domain to int_curves2.
        
        if (p.in_dom(*int_pts2[j]))
        {
          for (k = 0;
               k < num_int_curves2 && !int_curves2[k]->add_pt(int_pts2[j], 1);
               k++)
            ;
          
          assert(k < num_int_curves2);
          num_new_int_pts2[k]++;
          assert(num_new_int_pts2_other < MAX_NUM_NEW_INT_PTS);
          new_int_pts2_other[num_new_int_pts2_other]    = int_pts2[j];
          new_int_pts2_other[num_new_int_pts2_other]->ref_count++;
          ind_int_curves2_other[num_new_int_pts2_other] = k;
          num_new_int_pts2_other++;
        }
      }
      
//      for (j = 0; j < num_int_curves1; j++)
//      {
//        cerr << "   *int_curves1[" << j << "]->poly = " << endl << *int_curves1[j]->poly << endl << flush;
//        cerr << "   *int_curves1[" << j << "] = " << endl << *int_curves1[j] << endl << flush;
//      }
//      cerr << " ------------------------------------- " << endl << flush;
//      for (j = 0; j < num_int_curves2; j++)
//      {
//        cerr << "   *int_curves2[" << j << "]->poly = " << endl << *int_curves2[j]->poly << endl << flush;
//        cerr << "   *int_curves2[" << j << "] = " << endl << *int_curves2[j] << endl << flush;
//      }
//      cerr << " ------------------------------------- " << endl << flush;
    }
    else  //  if (!num_int_pts1)
    {
      num_int_pts2 = 0;
      int_pts2     = 0;
    }
    
    if (num_int_pts1 > 0)
    {
      for (j = 0; j < num_int_pts1; j++)
        if (!--int_pts1[j]->ref_count)
          delete int_pts1[j];
      
      delete [] int_pts1;
    }
    
    if (num_int_pts2 > 0)
    {
      for (j = 0; j < num_int_pts2; j++)
        if (!--int_pts2[j]->ref_count)
          delete int_pts2[j];
      
      delete [] int_pts2;
    }
  }
  
  for (i = 0; i < p.num_trim_curves; i++)
  {
    //  Find intersections of p.trim_curves and int_curve2 on p's domain.
    
    if (p.trim_curves[(i + 1) % p.num_trim_curves]->poly->
        eq_upto_const(int_curve2))
    //  The next trim_curve is identical to int_curve2.
    //  Don't count an intersection at the end.
      num_int_pts2 =
        p.trim_curves[i]->find_intersections(int_curve2, int_pts2, 0);
    else
    //  The next trim_curves is different from int_curve2.
    //  Count an intersection at the end.
      num_int_pts2 =
        p.trim_curves[i]->find_intersections(int_curve2, int_pts2, 1);
    
//    cerr << endl << "   int_curve2 = " << endl << int_curve2 << endl << flush;
//    cerr << "   *p.trim_curves[" << i << "]->poly = " << endl << *p.trim_curves[i]->poly << flush;
//    cerr << "   *p.trim_curves[" << i << "] = " << endl << *p.trim_curves[i] << flush;
//    cerr << " ------------------------------------- " << endl << flush;
//    cerr << "   int_curve2 has " << num_int_pts2 << " intersections with p.trim_curve[" << i << "]. " << endl << flush;
//    cerr << " ------------------------------------- " << endl << flush;
    
    if (num_int_pts2 > 0)
    {
      //  Find intersections on this domain.
      
//      num_int_pts1 =
//        get_all_pts(int_curve1,
//                    p.adj_surfs[i]->Impl->subst_param_expr(*surf->X,
//                                                           *surf->Y,
//                                                           *surf->Z,
//                                                           *surf->W),
//                    int_pts1,
//                    0);
      num_int_pts1 =
        get_all_int_pts(int_curve1, *p.adj_surfs[i], *surf, int_pts1);
//      cerr << "   num_int_pts2 = " << num_int_pts2 << ", num_int_pts1 = " << num_int_pts1 << endl << flush;
      assert(num_int_pts2 <= num_int_pts1);
      
      //  Match intersections on p's domain and this domain.
      
      match_pts(int_pts2, num_int_pts2, p.surf, int_pts1, num_int_pts1, surf);
      
//      for (j = 0; j < num_int_pts2; j++)
//        cerr << "   int_pts2[" << j << "] = " << *int_pts2[j] << ", int_pts1[" << j << "] = " << *int_pts1[j] << endl << flush;
//      cerr << " ------------------------------------- " << endl << flush;
      
      //  Add intersections to p.trim_curves, int_curves1 and int_curves2.
      
      for (j = 0; j < num_int_pts2; j++)
      {
        //  Associate intersections on this domain and p's domain.
        
        pt_surf_assoc.record_as_assoc(int_pts1[j], surf, int_pts2[j], p.surf);
        
        //  Add intersections on p's domain to p.trim_curves.
        
        assert(p.trim_curves[i]->add_pt(int_pts2[j]));
        
        //  Add intersections on this domain to int_curves1.
        
        if (in_dom(*int_pts1[j]))
        {
          for (k = 0;
               k < num_int_curves1 && !int_curves1[k]->add_pt(int_pts1[j]);
               k++)
            ;
          
          assert(k < num_int_curves1);          
          num_new_int_pts1[k]++;
          assert(num_new_int_pts1_other < MAX_NUM_NEW_INT_PTS);
          new_int_pts1_other[num_new_int_pts1_other]    = int_pts1[j];
          new_int_pts1_other[num_new_int_pts1_other]->ref_count++;
          ind_int_curves1_other[num_new_int_pts1_other] = k;
          num_new_int_pts1_other++;
        }
        
        //  Add intersections on p's domain to int_curves2.
        
        for (k = 0;
             k < num_int_curves2 && !int_curves2[k]->add_pt(int_pts2[j], 1);
             k++)
          ;
        
        assert(k < num_int_curves2);
        num_new_int_pts2[k]++;
        num_new_trim_pts2[k]++;
        assert(num_new_int_pts2_this < MAX_NUM_NEW_INT_PTS);
        new_int_pts2_this[num_new_int_pts2_this]    = int_pts2[j];
        new_int_pts2_this[num_new_int_pts2_this]->ref_count++;
        ind_int_curves2_this[num_new_int_pts2_this] = k;
        num_new_int_pts2_this++;
      }
      
//      for (j = 0; j < num_int_curves1; j++)
//      {
//        cerr << "   *int_curves1[" << j << "]->poly = " << endl << *int_curves1[j]->poly << endl << flush;
//        cerr << "   *int_curves1[" << j << "] = " << endl << *int_curves1[j] << endl << flush;
//      }
//      cerr << " ------------------------------------- " << endl << flush;
//      for (j = 0; j < num_int_curves2; j++)
//      {
//        cerr << "   *int_curves2[" << j << "]->poly = " << endl << *int_curves2[j]->poly << endl << flush;
//        cerr << "   *int_curves2[" << j << "] = " << endl << *int_curves2[j] << endl << flush;
//      }
//      cerr << " ------------------------------------- " << endl << flush;
    }
    else  //  if (!num_int_pts2)
    {
      num_int_pts1 = 0;
      int_pts1     = 0;
    }
    
    if (num_int_pts2 > 0)
    {
      for (j = 0; j < num_int_pts2; j++)
        if (!--int_pts2[j]->ref_count)
          delete int_pts2[j];
      
      delete [] int_pts2;
    }
    
    if (num_int_pts1 > 0)
    {
      for (j = 0; j < num_int_pts1; j++)
        if (!--int_pts1[j]->ref_count)
          delete int_pts1[j];
      
      delete [] int_pts1;
    }
  }
  
//  for (i = 0; i < num_trim_curves; i++)
//  {
//    cerr << "   *trim_curves[" << i << "]->poly = " << endl << *trim_curves[i]->poly << flush;
//    cerr << "   *trim_curves[" << i << "] = " << *trim_curves[i] << endl << flush;
//  }
//  cerr << " ------------------------------------- " << endl << flush;
//  for (i = 0; i < p.num_trim_curves; i++)
//  {
//    cerr << "   *p.trim_curves[" << i << "]->poly = " << endl << *p.trim_curves[i]->poly << flush;
//    cerr << "   *p.trim_curves[" << i << "] = " << *p.trim_curves[i] << endl << flush;
//  }
//  cerr << " ------------------------------------- " << endl << flush;
//  for (i = 0; i < num_int_curves1; i++)
//  {
//    cerr << "   *int_curves1[" << i << "]->poly = " << endl << *int_curves1[i]->poly << flush;
//    cerr << "   *int_curves1[" << i << "] = " << *int_curves1[i] << endl << flush;
//  }
//  cerr << " ------------------------------------- " << endl << flush;
//  for (i = 0; i < num_int_curves2; i++)
//  {
//    cerr << "   *int_curves2[" << i << "]->poly = " << endl << *int_curves2[i]->poly << flush;
//    cerr << "   *int_curves2[" << i << "] = " << *int_curves2[i] << endl << flush;
//  }
//  cerr << " ------------------------------------- " << endl << flush;
  
  //  4. Determine direction-correspondences.
  
  int           f_0, f_1, f_2;
  K_POINT2D*    x1_0;
  K_POINT2D*    x2_0;
  unsigned long c_0;
  K_POINT2D*    x1_1;
  K_POINT2D*    x2_1;
  unsigned long c_1;
  K_BOXCO2      b;
  unsigned long split_dir1;
  bigrational   split_v1;
  K_RATPOLY     split_poly1;
  K_POINT2D**   split_pts1;
  unsigned long num_split_pts1;
  K_POINT2D*    x1_2;
  K_SURF        split_surf1;
  K_RATPOLY     split_poly2;
  K_POINT2D**   split_pts2;
  unsigned long num_split_pts2;
  K_POINT2D*    x2_2;
  K_POINT2D**   pts_sort;
  int           dir1, dir2;
  bigrational   split_v1_inc_step;
  K_POINT2D**   loop_pts1;
  K_POINT2D**   loop_pts2;
  int**         dir_corr;
  
  if (num_int_curves1 > 0)
  {
    dir_corr = new int* [num_int_curves1];
    
    for (i = 0; i < num_int_curves1; i++)
    {
      if (num_int_curves2 > 0)
        dir_corr[i] = new int [num_int_curves2];
      else
        dir_corr[i] = 0;
      
      for (j = 0; j < num_int_curves2; j++)
        dir_corr[i][j] = 0;
    }
  }
  else
    dir_corr = 0;
  
//  cerr << " Debug3: kpatch: intersect: ---------- " << endl << flush;
//  for (i = 0; i < num_int_curves1; i++)
//    cerr << "   num_new_int_pts1[" << i << "] = " << num_new_int_pts1[i] << endl << flush;
//  cerr << " ------------------------------------- " << endl << flush;
  
  for (i = 0; i < num_int_curves1; i++)
  //  Find 3 pts lie on curves on this domain.
  //  Also, find 3 pts lie on curves on p's domain.
  //  Sort 3 pts on each domain.
  //  Determine direction correspondence by seeing whether or not
  //    the middle one on this domain is also the middle one on p's domain.
    if (num_new_int_pts1[i] > 0)
    {
      //  intersections have been added to int_curves1[i]
      
      f_0 = f_1 = 0;
      
      for (j = 0; j < int_curves1[i]->num_segments; j++)
      {
//        cerr << "   3pts: j = " << j << ", f_0 = " << f_0 << ", f_1 = " << f_1 << endl << flush;
        
        if (!f_0)
        //  Find the 1st pts; *x1_0 on this domain and *x2_0 on p's domain.
        {
          for (k = 0;
               k < num_new_int_pts1_this
               &&
               !new_int_pts1_this[k]->equiv(*int_curves1[i]->segments[j]->end);
               k++)
            ;
          //  Every new_int_pts1_this[k] is added to some curve,
          //    i.e., it is an endpt of some segment of the curve.
          //  Thus, it is safe to call equiv instead of equal.
          
          if (k < num_new_int_pts1_this)
          //  if *new_int_pts1_this[k] is on *int_curves[i]->segments[j]
          {
            x2_0 = pt_surf_assoc.find_assoc_pt(new_int_pts1_this[k],
                                               surf, p.surf);
            
            if (p.in_dom(*x2_0))
            //  if *new_int_pts1_this[k] is associated some pt on p's domain
            {
              x1_0 = new_int_pts1_this[k];
              
//              cerr << "   3pts: *x1_0 = " << *x1_0 << ", *x2_0 = " << *x2_0 << endl << flush;
              
              for (l = 0;
                   l < num_new_int_pts2_other
                   &&
                   !new_int_pts2_other[l]->equiv(*x2_0);
                   l++)
                ;
              //  Every new_int_pts2_other[l] is added to some curve, and
              //  Every pt in pt_surf_assoc is also added to some curve.
              //  Thus, it is safe to call equiv instead of equal.
              
              assert(l < num_new_int_pts2_other);
              c_0 = ind_int_curves2_other[l];
              
              if (!dir_corr[i][c_0])
              //  dir_corr[i][c_0] has not yet been determined.
              //  We use *x1_0 and *x2_0 to determine dir_corr[i][c_0].
                f_0 = 1;
              else  //  if (dir_corr[i][c_0])
              //  dir_corr[i][c_0] has already been determined.
              //  Discard *x1_0 and *x2_0.
                f_0 = 0;
            }
            else  //  if (p.in_dom(*x2_0))
            //  *x2_0 does not lie on p's domain. Discard it.
              f_0 = 0;
          }
          else  //  if (k == num_new_int_pts1_this)
          //  if *new_int_pts1_this[k] is not on *int_curves[i]->segments[j]
          {
            for (k = 0;
                 k < num_new_int_pts1_other
                 &&
                 !new_int_pts1_other[k]->
                   equiv(*int_curves1[i]->segments[j]->end);
                 k++)
              ;
            //  Every new_int_pts1_other[k] is added to some curve,
            //    i.e., it is an endpt of some segment of the curve.
            //  Thus, it is safe to call equiv instead of equal.
            
            if (k < num_new_int_pts1_other)
            //  if *new_int_pts1_other[k] is on *int_curves[i]->segments[j]
            {
              x2_0 = pt_surf_assoc.find_assoc_pt(new_int_pts1_other[k],
                                                 surf, p.surf);
              
              if (p.in_dom(*x2_0))
              //  if *new_int_pts1_other[k] is associated some pt on p's domain
              {
                x1_0 = new_int_pts1_other[k];
                
//                cerr << "   3pts: *x1_0 = " << *x1_0 << ", *x2_0 = " << *x2_0 << endl << flush;
                
                for (l = 0;
                     l < num_new_int_pts2_this
                     &&
                     !new_int_pts2_this[l]->equiv(*x2_0);
                     l++)
                  ;
                //  Every new_int_pts2_this[l] is added to some curve, and,
                //  Every pt in pt_surf_assoc is also added to some curve.
                //  Thus, it is safe to call equiv instead of equal.
                
                assert(l < num_new_int_pts2_this);
                c_0 = ind_int_curves2_this[l];
                
                if (!dir_corr[i][c_0])
                //  dir_corr[i][c_0] has not yet been determined.
                //  We use *x1_0 and *x2_0 to determine dir_corr[i][c_0].
                  f_0 = 1;
                else  //  if (dir_corr[i][c_0])
                //  dir_corr[i][c_0] has already been determined.
                //  Discard *x1_0 and *x2_0.
                  f_0 = 0;
              }
              else  //  if (!p.in_dom(*x2_0))
              //  *x2_0 does not lie on p's domain. Discard it.
                f_0 = 0;
            }
            else  //  if (k == num_new_int_pts1_this)
            //  No pt on *int_curves1[i]->segments[j] can be used
            //    for *x1_0 or *x2_0.
              f_0 = 0;
          }
        }
        else  //  if (f_0)
        //  The 1st pts, x1_0 and x2_0, have already been found.
        //  Find the 2nd pts, *x1_1 on this domain and *x2_1 on p's domain.
        {
          for (k = 0;
               k < num_new_int_pts1_this
               &&
               !new_int_pts1_this[k]->equiv(*int_curves1[i]->segments[j]->end);
               k++)
            ;
          //  Every new_int_pts1_this[k] is added to some curve,
          //    i.e., it is an endpt of some segment of the curve.
          //  Thus, it is safe to call equiv instead of equal.
          
          if (k < num_new_int_pts1_this)
          //  if *new_int_pts1_this[k] is on *int_curves[i]->segments[j]
          {
            x2_1 = pt_surf_assoc.find_assoc_pt(new_int_pts1_this[k],
                                               surf, p.surf);
            
            if (p.in_dom(*x2_1))
            //  if *new_int_pts1_this[k] is associated some pt on p's domain
            {
              x1_1 = new_int_pts1_this[k];
              
//              cerr << "   3pts: *x1_1 = " << *x1_1 << ", *x2_1 = " << *x2_1 << endl << flush;
              
              for (l = 0;
                   l < num_new_int_pts2_other
                   &&
                   !new_int_pts2_other[l]->equiv(*x2_1);
                   l++)
                ;
              //  Every new_int_pts2_other[l] is added to some curve, and
              //  Every pt in pt_surf_assoc is also added to some curve.
              //  Thus, it is safe to call equiv instead of equal.
              
              assert(l < num_new_int_pts2_other);
              c_1 = ind_int_curves2_other[l];
              f_1 = 1;
            }
            else  //  if (!p.in_dom(*x2_1))
            //  *x2_1 does not lie on p's domain. Discard it.
              f_1 = 0;
          }
          else  //  if (k == num_new_int_pts1_this)
          //  if *new_int_pts1_this[k] is not on *int_curves[i]->segments[j]
          {
            for (k = 0;
                 k < num_new_int_pts1_other
                 &&
                 !new_int_pts1_other[k]->
                   equiv(*int_curves1[i]->segments[j]->end);
                 k++)
              ;
            //  Every new_int_pts1_other[k] is added to some curve,
            //    i.e., it is an endpt of some segment of the curve.
            //  Thus, it is safe to call equiv instead of equal.
            
            if (k < num_new_int_pts1_other)
            //  if *new_int_pts1_other[k] is on *int_curves[i]->segments[j]
            {
              x2_1 = pt_surf_assoc.find_assoc_pt(new_int_pts1_other[k],
                                                 surf, p.surf);
              
              if (p.in_dom(*x2_1))
              //  if *new_int_pts1_other[k] is associated some pt on p's domain
              {
                x1_1 = new_int_pts1_other[k];
                
//                cerr << "   3pts: *x1_1 = " << *x1_1 << ", *x2_1 = " << *x2_1 << endl << flush;
                
                for (l = 0;
                     l < num_new_int_pts2_this
                     &&
                     !new_int_pts2_this[l]->equiv(*x2_1);
                     l++)
                  ;
                //  Every new_int_pts2_this[l] is added to some curve, and
                //  Every pt in pt_surf_assoc is also added to some curve.
                //  Thus, it is safe to call equiv instead of equal.
                
                assert(l < num_new_int_pts2_this);
                c_1 = ind_int_curves2_this[l];
                f_1 = 1;
              }
              else  //  if (!p.in_dom(*x2_1))
              //  *x2_1 does not lie on p's domain. Discard it.
                f_1 = 0;
            }
            else  //  if (k == num_new_int_pts1_other)
            //  No pt on *int_curves1[i]->segments[j] can be used
            //    for *x1_1 or *x2_1.
              f_1 = 0;
          }
          
          if (f_1)
          //  The 1st and 2nd pts have been found.
          //  Generate the 3rd pt, *x1_2, on this domain s.t.
          //    it lies between *x1_0 and *x1_1.
          //  Find the 3rd pt *x2_2 on p's domain by inverting *x1_2.
          {
//            cerr << "   3pts: c_0 = " << c_0 << ", c_1 = " << c_1 << endl << flush;
            assert(c_0 == c_1);
            
            f_2 = 0;
            
            while (!f_2)
            {
//              cerr << "   3pts: *x1_0 = " << *x1_0 << endl << flush;
//              cerr << "   3pts: *x1_1 = " << *x1_1 << endl << flush;
//              cerr << "   3pts: -------------------- " << endl << flush;
              
              b = x1_0->bbox().merge(x1_1->bbox());
              
              if (b.high[0] - b.low[0] > b.high[1] - b.low[1])
              {
                split_dir1 = 0;
                split_v1   = (b.high[0] + b.low[0]) / 2;
              }
              else  //  if (b.high[0] - b.low[0] <= b.high[1] - b.low[1])
              {
                split_dir1 = 1;
                split_v1   = (b.high[1] + b.low[1]) / 2;
              }
              
              while (!(num_split_pts1 =
                         int_curves1[i]->find_intersections(
                         K_RATPOLY(2, split_dir1, split_v1), split_pts1, 1)))
              {
                //  b is too wide. Shrink pts until b is appropriately wide.
                
                x1_0->shrink(shrink_step, shrink_step);
                x1_1->shrink(shrink_step, shrink_step);
                
                b = x1_0->bbox().merge(x1_1->bbox());
                
                if (b.high[0] - b.low[0] > b.high[1] - b.low[1])
                {
                  split_dir1 = 0;
                  split_v1   = (b.high[0] + b.low[0]) / 2;
                }
                else  //  if (b.high[0] - b.low[0] <= b.high[1] - b.low[1])
                {
                  split_dir1 = 1;
                  split_v1   = (b.high[1] + b.low[1]) / 2;
                }
              }
              
              x1_2 = split_pts1[0];
              f_2  = 1;
              
//              cerr << "   3pts: *x1_0 = " << *x1_0 << endl << flush;
//              cerr << "   3pts: *x1_1 = " << *x1_1 << endl << flush;
//              cerr << "   3pts: *x1_2 = " << *x1_2 << endl << flush;
//              cerr << "   3pts: -------------------- " << endl << flush;
              
              if (x1_0->overlap(*x1_2))
              {
//                cerr << "   3pts: *x1_0 and *x1_2 overlap. " << endl << flush;
                
                //  *x1_0 and *x1_2 overlap.
                //  Cut *x1_0 by *x1_2 and recompute *x1_2.
                
                x1_0->cut_s(x1_2->get_low_s());
                x1_0->cut_s(x1_2->get_high_s());
                x1_0->cut_t(x1_2->get_low_t());
                x1_0->cut_t(x1_2->get_high_t());
                f_2 = 0;
              }
              
              if (x1_1->overlap(*x1_2))
              {
//                cerr << "   3pts: *x1_1 and *x1_2 overlap. " << endl << flush;
                
                //  *x1_1 and *x1_2 overlap.
                //  Cut *x1_1 by *x1_2 and recompute *x1_2.
                
                x1_1->cut_s(x1_2->get_low_s());
                x1_1->cut_s(x1_2->get_high_s());
                x1_1->cut_t(x1_2->get_low_t());
                x1_1->cut_t(x1_2->get_high_t());
                f_2 = 0;
              }
            }
            
            //  Compute a surface, split_surf1, which splits *this->surf.
            //  Compute an intersection, split_poly2,
            //     of split_surf1 and *p.surf.
            //  Compute intersections of *int_curves2[c_0] and split_poly2.
            //   *x2_2 is one of them which matches to *x1_2.
            
            split_surf1    = surf->split_surf(split_v1, split_dir1);
//            split_poly2    = split_surf1.Impl->subst_param_expr(*p.surf->X,
//                                                                *p.surf->Y,
//                                                                *p.surf->Z,
//                                                                *p.surf->W);
//            num_split_pts2 = get_all_pts(*int_curves2[c_0]->poly, split_poly2,
//                                         split_pts2,
//                                         0);
            num_split_pts2 = get_all_int_pts(*int_curves2[c_0]->poly,
                                             split_surf1, *p.surf,
                                             split_pts2);
            assert(num_split_pts2 > 0);
            match_pts(split_pts1, 1, surf, split_pts2, num_split_pts2, p.surf);
            x2_2           = split_pts2[0];
            
            //  Determine direction correspondences.
            
            //  *x1_2 does not necessarily lie between *x1_0 and *x1_1
            //    (although it does in most cases).
            //  We must sort them.
            
            pts_sort = new K_POINT2D* [3];
            pts_sort[0] = x1_0;
            pts_sort[0]->ref_count++;
            pts_sort[1] = x1_1;
            pts_sort[1]->ref_count++;
            pts_sort[2] = x1_2;
            pts_sort[2]->ref_count++;
            
//            cerr << "   3pts: before sort: ---------- " << endl << flush;
//            cerr << "                *x1_0 = " << *x1_0 << endl << flush;
//            cerr << "                *x1_1 = " << *x1_1 << endl << flush;
//            cerr << "                *x1_2 = " << *x1_2 << endl << flush;
//            cerr << "         *pts_sort[0] = " << *pts_sort[0] << endl << flush;
//            cerr << "         *pts_sort[1] = " << *pts_sort[1] << endl << flush;
//            cerr << "         *pts_sort[2] = " << *pts_sort[2] << endl << flush;
//            cerr << "   3pts: ----------------------- " << endl << flush;
            
            int_curves1[i]->sort_pts(pts_sort, 3);
            
//            cerr << "   3pts: after sort: ----------- " << endl << flush;
//            cerr << "                *x1_0 = " << *x1_0 << endl << flush;
//            cerr << "                *x1_1 = " << *x1_1 << endl << flush;
//            cerr << "                *x1_2 = " << *x1_2 << endl << flush;
//            cerr << "         *pts_sort[0] = " << *pts_sort[0] << endl << flush;
//            cerr << "         *pts_sort[1] = " << *pts_sort[1] << endl << flush;
//            cerr << "         *pts_sort[2] = " << *pts_sort[2] << endl << flush;
//            cerr << "   3pts: ----------------------- " << endl << flush;
            
            if (pts_sort[0] == x1_0 && pts_sort[1] == x1_1)
            //                      && pts_sort[2] == x1_2)
              dir1 = 1;
            else if (pts_sort[0] == x1_0 && pts_sort[1] == x1_2)
            //                           && pts_sort[2] == x1_1)
              dir1 = - 1;
            else if (pts_sort[0] == x1_1 && pts_sort[1] == x1_2)
            //                           && pts_sort[2] == x1_0)
              dir1 = 1;
            else if (pts_sort[0] == x1_1 && pts_sort[1] == x1_0)
            //                           && pts_sort[2] == x1_2)
              dir1 = - 1;
            else if (pts_sort[0] == x1_2 && pts_sort[1] == x1_0)
            //                           && pts_sort[2] == x1_1)
              dir1 = 1;
            else  //  if (pts_sort[0] == x1_2 && pts_sort[1] == x1_1
                  //                          && pts_sort[2] == x1_0)
              dir1 = - 1;
            
//            cerr << "   3pts: dir1 = " << dir1 << endl << flush;
            
            for (k = 0; k < 3; k++)
              if (!--pts_sort[k]->ref_count)
                delete pts_sort[k];
            
            delete [] pts_sort;
            
            //  Sort *x2_0, *x2_1 and *x2_2.
            
            if (int_curves2[c_0]->contains(*x2_2))
            //  if *x2_2 lies on p's domain
            {
              pts_sort    = new K_POINT2D* [3];
              pts_sort[0] = x2_0;
              pts_sort[0]->ref_count++;
              pts_sort[1] = x2_1;
              pts_sort[1]->ref_count++;
              pts_sort[2] = x2_2;
              pts_sort[2]->ref_count++;
              
              int_curves2[c_0]->sort_pts(pts_sort, 3);
              
              if (pts_sort[0] == x2_0 && pts_sort[1] == x2_1)
              //                      && pts_sort[2] == x2_2)
                dir2 = 1;
              else if (pts_sort[0] == x2_0 && pts_sort[1] == x2_2)
              //                           && pts_sort[2] == x2_1)
                dir2 = - 1;
              else if (pts_sort[0] == x2_1 && pts_sort[1] == x2_2)
              //                           && pts_sort[2] == x2_0)
                dir2 = 1;
              else if (pts_sort[0] == x2_1 && pts_sort[1] == x2_0)
              //                           && pts_sort[2] == x2_2)
                dir2 = - 1;
              else if (pts_sort[0] == x2_2 && pts_sort[1] == x2_0)
              //                           && pts_sort[2] == x2_1)
                dir2 = 1;
              else  //  if (pts_sort[0] == x2_2 && pts_sort[1] == x2_1
                    //                          && pts_sort[2] == x2_0)
                dir2 = - 1;
              
              for (k = 0; k < 3; k++)
                if (!--pts_sort[k]->ref_count)
                  delete pts_sort[k];
              
              delete [] pts_sort;
            }
            else  // if (!int_curves2[c_0]->contains(*x2_2))
            //  if *x2_2 does not lie on p's domain
            {
              k = int_curves2[c_0]->locate_pt_seg_end(*x2_0);
              assert(k < int_curves2[c_0]->num_segments);
              l = int_curves2[c_0]->locate_pt_seg_end(*x2_1);
              assert(l < int_curves2[c_0]->num_segments);
              
              if (k > l)
              //  *x2_1, *x2_0, *x2_2 or *x2_2, *x2_1, *x2_0
                dir2 = - 1;
              else  //  if (k <= l)
              //  *x2_0, *x2_1, *x2_2 or *x2_2, *x2_0, *x2_1
                dir2 = 1;
              
//              cerr << "   3pts: dir2 = " << dir2 << endl << flush;
            }
            
            dir_corr[i][c_0] = dir1 * dir2;
//            cerr << "   int_curves1[" << i << "]->segments[" << j << "] " << endl << flush;
//            cerr << "   dir1 = " << dir1 << ", dir2 = " << dir2 << endl << flush;
//            cerr << "   dir_corr[" << i << "][" << c_0 << "] = " << dir_corr[i][c_0] << endl << flush;
//            cerr << " ------------------------------------- " << endl << flush;
            
            for (k = 0; k < num_split_pts1; k++)
              if (!--split_pts1[k]->ref_count)
                delete split_pts1[k];
            
            delete [] split_pts1;
            
            for (k = 0; k < num_split_pts2; k++)
              if (!--split_pts2[k]->ref_count)
                delete split_pts2[k];
            
            delete [] split_pts2;
            
            f_0 = f_1 = 0;
          }
        }
      }
    }
    else  //  if (!num_new_int_pts1[i])
    //  if *int_curves1[i] does not have any intersection
      if (int_curves1[i]->is_closed())
      //  *int_curves1[i] is a loop.
      //  See whether or not there is any corresponding loop on p's domain.
      //  Generate a pt on the loop, *int_curves1[i].
      //  If we can invert the pt to p's domain then
      //    there exists a curve on which the inverted pt lies.and
      //    the curve is a corresponding loop.
      //  Otherwise, there is no corresponding loop on p's domain.
      {
        b = int_curves1[i]->bbox();
        
        if (b.high[0] - b.low[0] > b.high[1] - b.low[1])
        {
          split_dir1        = 0;
          split_v1_inc_step = (b.high[0] - b.low[0]) / 10;
//          split_v1          = (b.high[0] + b.low[0]) / 2 + split_v1_inc_step;
          split_v1          = b.low[0] + split_v1_inc_step;
        }
        else  //  if (b.high[0] - b.low[0] <= b.high[1] - b.low[1])
        {
          split_dir1        = 1;
          split_v1_inc_step = (b.high[1] - b.low[1]) / 10;
//          split_v1          = (b.high[1] + b.low[1]) / 2 + split_v1_inc_step;
          split_v1          = b.low[1] + split_v1_inc_step;
        }
        
        //  Generate pts on the loop, *int_curves[i].
        
        while (!(num_int_pts1 =
                   int_curves1[i]->find_intersections(
                     K_RATPOLY(2, split_dir1, split_v1), int_pts1, 1)))
          split_v1 += split_v1_inc_step;
        
        //  Compute a surface, split_surf1 which splits *this->surf.
        //  Compute an intersection, split_poly2,
        //     of split_surf1 and *p.surf.
        //  Compute intersections of int_curve2 and split_poly2.
        //  Match pts on the loop, *int_curves1[i], and
        //        the intersections of int_curve2 and split_poly2.
        
        split_surf1 = surf->split_surf(split_v1, split_dir1);
//        split_poly2 = split_surf1.Impl->subst_param_expr(*p.surf->X,
//                                                         *p.surf->Y,
//                                                         *p.surf->Z,
//                                                         *p.surf->W);
//        num_int_pts2 = get_all_pts(int_curve2, split_poly2, int_pts2, 0);
        num_int_pts2 = get_all_int_pts(int_curve2, split_surf1, *p.surf,
                                       int_pts2);
        assert(num_int_pts2 > 0);
        match_pts(int_pts1, num_int_pts1, surf,
                  int_pts2, num_int_pts2, p.surf);
        
        if (p.contains(*int_pts2[0]))
        //  *int_pts2[0] is on p's domain.
        //  Thus, the loop containing *int_pts2[0] entirely lie on p's domain.
        {
          //  Find 3 pts on the loop on this domain and
          //       3 pts on the loop on p's domain.
          
          loop_pts1 = new K_POINT2D* [3];
          loop_pts2 = new K_POINT2D* [3];
          
          //  Find the loop which contains *int_pts2[0].
          
          for (j = 0;
               j < num_int_curves2 && !int_curves2[j]->contains(*int_pts2[0]);
               j++)
            ;
          
          for (k = 0; k < 3 && k < num_int_pts1; k++)
          {
            loop_pts1[k] = int_pts1[k];
            loop_pts1[k]->ref_count++;
            loop_pts2[k] = int_pts2[k];
            loop_pts2[k]->ref_count++;
          }
          
          while (k < 3)
          //  We have not yet found enough pts.
          {
            for (l = 0; l < num_int_pts1; l++)
              if (!--int_pts1[l]->ref_count)
                delete int_pts1[l];
            
            delete [] int_pts1;
            
            for (l = 0; l < num_int_pts2; l++)
              if (!--int_pts2[l]->ref_count)
                delete int_pts2[l];
            
            delete [] int_pts2;
            
            while (!(num_int_pts1 =
                       int_curves1[i]->find_intersections(
                         K_RATPOLY(2, split_dir1, split_v1), int_pts1, 1)))
              split_v1 += split_v1_inc_step;
            
            split_surf1 = surf->split_surf(split_v1, split_dir1);
//            split_poly2 = split_surf1.Impl->subst_param_expr(*p.surf->X,
//                                                             *p.surf->Y,
//                                                             *p.surf->Z,
//                                                             *p.surf->W);
//            num_int_pts2 = get_all_int_pts(int_curve2, split_poly2,
//                                           int_pts2,
//                                           0);
            num_int_pts2 = get_all_int_pts(int_curve2, split_surf1, *p.surf,
                                           int_pts2);
            assert(num_int_pts2 > 0);
            match_pts(int_pts1, num_int_pts1, surf,
                      int_pts2, num_int_pts2, p.surf);
            
            for (l = 0; k < 3 && l < num_int_pts1; k++, l++)
            {
              loop_pts1[k] = int_pts1[l];
              loop_pts1[k]->ref_count++;
              loop_pts2[k] = int_pts2[l];
              loop_pts2[k]->ref_count++;
            }
          }
          
          //  Determine direction correspondences.
          
          //  Sort loop_pts1's
          //    to determine the direction of the loop on this domain.
          
          pts_sort    = new K_POINT2D* [3];
          pts_sort[0] = loop_pts1[0];
          pts_sort[0]->ref_count++;
          pts_sort[1] = loop_pts1[1];
          pts_sort[1]->ref_count++;
          pts_sort[2] = loop_pts1[2];
          pts_sort[2]->ref_count++;
          
          int_curves1[i]->sort_pts(pts_sort, 3);
          
          if (pts_sort[0] == loop_pts1[0] && pts_sort[1] == loop_pts1[1])
          //                              && pts_sort[2] == loop_pts1[2])
            dir1 = 1;
          else if (pts_sort[0] == loop_pts1[0] && pts_sort[1] == loop_pts1[2])
          //                                   && pts_sort[2] == loop_pts1[1])
            dir1 = - 1;
          else if (pts_sort[0] == loop_pts1[1] && pts_sort[1] == loop_pts1[2])
          //                                   && pts_sort[2] == loop_pts1[0])
            dir1 = 1;
          else if (pts_sort[0] == loop_pts1[1] && pts_sort[1] == loop_pts1[0])
          //                                   && pts_sort[2] == loop_pts1[2])
            dir1 = - 1;
          else if (pts_sort[0] == loop_pts1[2] && pts_sort[1] == loop_pts1[0])
          //                                   && pts_sort[2] == loop_pts1[1])
            dir1 = 1;
          else  //  if (pts_sort[0] == loop_pts1[2] &&
                //      pts_sort[1] == loop_pts1[1] &&
                //      pts_sort[2] == loop_pts1[0])
            dir1 = - 1;
          
          for (k = 0; k < 3; k++)
            if (!--pts_sort[k]->ref_count)
              delete pts_sort[k];
          
          delete [] pts_sort;
          
          //  Sort loop_pts2's
          //    to determine the direction of the loop on p's domain.
          
          pts_sort    = new K_POINT2D* [3];
          pts_sort[0] = loop_pts2[0];
          pts_sort[0]->ref_count++;
          pts_sort[1] = loop_pts2[1];
          pts_sort[1]->ref_count++;
          pts_sort[2] = loop_pts2[2];
          pts_sort[2]->ref_count++;
          
          int_curves2[j]->sort_pts(pts_sort, 3);
          
          if (pts_sort[0] == loop_pts2[0] && pts_sort[1] == loop_pts2[1])
          //                              && pts_sort[2] == loop_pts2[2])
            dir2 = 1;
          else if (pts_sort[0] == loop_pts2[0] && pts_sort[1] == loop_pts2[2])
          //                                   && pts_sort[2] == loop_pts2[1])
            dir2 = - 1;
          else if (pts_sort[0] == loop_pts2[1] && pts_sort[1] == loop_pts2[2])
          //                                   && pts_sort[2] == loop_pts2[0])
            dir2 = 1;
          else if (pts_sort[0] == loop_pts2[1] && pts_sort[1] == loop_pts2[0])
          //                                   && pts_sort[2] == loop_pts2[2])
            dir2 = - 1;
          else if (pts_sort[0] == loop_pts2[2] && pts_sort[1] == loop_pts2[0])
          //                                   && pts_sort[2] == loop_pts2[1])
            dir2 = 1;
          else  //  if (pts_sort[0] == loop_pts2[2] &&
                //      pts_sort[1] == loop_pts2[1] &&
                //      pts_sort[2] == looo_pts2[0])
            dir2 = - 1;
          
          for (k = 0; k < 3; k++)
            if (!--pts_sort[k]->ref_count)
              delete pts_sort[k];
          
          delete [] pts_sort;
          
          dir_corr[i][j] = dir1 * dir2;
          
//          cerr << "   int_curves1[" << i << "] is closed" << endl << flush;
//          cerr << "   dir1 = " << dir1 << ", dir2 = " << dir2 << endl << flush;
//          cerr << "   dir_corr[" << i << "][" << j << "] = " << dir_corr[i][j] << endl << flush;
//          cerr << " ------------------------------------- " << endl << flush;
          
          //  Add those pts used to determine direction correspondences
          //    to the loops.
          
          for (k = 0; k < 3; k++)
          {
            int_curves1[i]->add_pt(loop_pts1[k]);
            int_curves2[j]->add_pt(loop_pts2[k]);
          }
          
//          for (unsigned long k1 = 0; k1 < num_int_curves1; k1++)
//          {
//            cerr << "   *int_curves1[" << k1 << "]->poly = " << endl << int_curves1[k1]->poly << endl << flush;
//            cerr << "   *int_curves1[" << k1 << "] = " << endl << *int_curves1[k1] << endl << flush;
//          }
//          cerr << " ------------------------------------- " << endl << flush;
//          for (unsigned long k1 = 0; k1 < num_int_curves2; k1++)
//          {
//            cerr << "   *int_curves2[" << k1 << "]->poly = " << endl << int_curves2[k1]->poly << endl << flush;
//            cerr << "   *int_curves2[" << k1 << "] = " << endl << *int_curves2[k1] << endl << flush;
//          }
//          cerr << " ------------------------------------- " << endl << flush;
          
          //  Rotate the loops
          //    because the added pts ate artificial endpts and
          //            the loops should not start from them.
          
          k = int_curves1[i]->locate_pt_seg_end(*loop_pts1[0]);
          assert(k < int_curves1[i]->num_segments);
          int_curves1[i]->rotate_closed_curve(k + 1);
          
          k = int_curves2[j]->locate_pt_seg_end(*loop_pts2[0]);
          assert(k < int_curves2[j]->num_segments);
          int_curves2[j]->rotate_closed_curve(k + 1);
          
          for (k = 0; k < 3; k++)
            if (!--loop_pts1[k]->ref_count)
              delete loop_pts1[k];
          
          delete [] loop_pts1;
          
          for (k = 0; k < 3; k++)
            if (!--loop_pts2[k]->ref_count)
              delete loop_pts2[k];
          
          delete [] loop_pts2;
        }
        
        for (j = 0; j < num_int_pts1; j++)
          if (!--int_pts1[j]->ref_count)
            delete int_pts1[j];
        
        delete [] int_pts1;
        
        for (j = 0; j < num_int_pts2; j++)
          if (!--int_pts2[j]->ref_count)
            delete int_pts2[j];
        
        delete [] int_pts2;
      }
  
//  for (i = 0; i < num_int_curves1; i++)
//    for (j = 0; j < num_int_curves2; j++)
//      cerr << "   dir_corr[" << i << "][" << j << "] = " << dir_corr[i][j] << endl << flush;
//  cerr << " ------------------------------------- " << endl << flush;
  
  //  5. Mark segments good or bad.
  
  K_POINT2D*    other_pt;
  unsigned long c_2, s_1, s_2;
  unsigned long o_c, o_s;
  int           in_on_out_region, is_fw;
  int           g;
  int**         good1;
  int**         good2;
  
  if (num_int_curves1 > 0)
  {
    good1 = new int* [num_int_curves1];
    
    for (i = 0; i < num_int_curves1; i++)
    {
      good1[i] = new int [int_curves1[i]->num_segments];
      
      for (j = 0; j < int_curves1[i]->num_segments; j++)
        good1[i][j] = 1;
    }
  }
  else
    good1 = 0;
  
  if (num_int_curves2 > 0)
  {
    good2 = new int* [num_int_curves2];
    
    for (i = 0; i < num_int_curves2; i++)
    {
      good2[i] = new int [int_curves2[i]->num_segments];
      
      for (j = 0; j < int_curves2[i]->num_segments; j++)
        good2[i][j] = 1;
    }
  }
  else
    good2 = 0;
  
//  cerr << " Debug4: kpatch: intersect: ---------- " << endl << flush;
  
  for (i = 0; i < num_int_curves1; i++)
  {
    //  See if there exists a curve lying entirely in/outside the trim region,
    //    i.e., a curve having no intersection with the trim_curves.
    //  If such a curve exists and it also lies outside the trim region then
    //    mark all its segments bad.
    
//    cerr << "   clipping: 00: i = " << i << ", num_int_curves1 = " << num_int_curves1 << ", num_new_trim_pts1[" << i << "] = " << num_new_trim_pts1[i] << endl << flush;
    
    if (!num_new_trim_pts1[i])
    //  if *int_curves1[i] has no intersection with the trim_curves
      if (int_curves1[i]->is_closed())
      //  *int_curves1[i] is a loop.
      //  See if the loop lies entirely in the trim region;
      //  Pick some pts and see if each of them is IN/OUTSIDE the trim region.
      //  We knew that *int_curves1[i] does not hit the trim_curves.
      //  Thus, it is safe to call contains.
        if (!contains(*int_curves1[i]->start()))
        //  *int_curves1[i] lies entirely OUTSIDE the trim region.
        //  Mark all the segments bad and
        //    remove the corresponding loop on p's domain.
        {
          for (j = 0; j < int_curves1[i]->num_segments; j++)
            good1[i][j] = 0;
          
          if (!num_new_int_pts1[i])
          {
            for (j = 0; j < num_int_curves2; j++)
              if (dir_corr[i][j])
              //  *int_curves1[i] and *int_curves2[j] correspond and
              //    do not have intersections corresponding with each other.
              //  Mark all the segments of *int_curves2[j] bad.
                for (k = 0; k < int_curves2[j]->num_segments; k++)
                  good2[j][k] = 0;
          }
          else  //  if (num_new_int_pts1[i] > 0)
            for (j = 0; j < num_new_int_pts1_other; j++)
              if (ind_int_curves1_other[j] == i)
              //  *new_int_pts1_other[j] lies on *int_curves1[i] and
              //    corresponds to some pt, *other_pt, on p's domain.
              //  Find a curve on p's domain on which *other_pt lies.
              {
                other_pt = pt_surf_assoc.find_assoc_pt(new_int_pts1_other[j],
                                                       surf, p.surf);
                
                for (k = 0;
                     k < num_new_int_pts2_this
                     &&
                     !new_int_pts2_this[k]->equiv(*other_pt);
                     k++)
                  ;
                //  Every new_int_pts2_this[k] is added to some curve, and
                //  Every pt in pt_surf_assoc is also added to some curve.
                //  Thus, it is safe to call equiv instead of equal.
                
                assert(k < num_new_int_pts2_this);
                c_2 = ind_int_curves2_this[k];
                s_2 = int_curves2[c_2]->locate_pt_seg_end(*other_pt);
                
                //  On p's domain, march along the curve to mark bad segments.
                //  Go forward from segments[s_2 + 1] as long as
                //    their START's are NOT new_int_pts2_this/other's, and
                //  go backward from segments[s_2] as long as
                //    their END's are NOT new_int_pts2_this/other's.
                
                if (s_2 < int_curves2[c_2]->num_segments - 1)
                //  Go forward from segments[s_2 + 1] as long as
                //    their START's are NOT new_int_pts2_this/other's.
                  int_curves2[c_2]->mk_seg_fw(good2[c_2], 0,
                                              s_2 + 1,
                                              new_int_pts2_this,
                                              num_new_int_pts2_this,
                                              new_int_pts2_other,
                                              num_new_int_pts2_other);
                else  //  if (s_2 == int_curves2[c_2]->num_segments - 1)
                //  s_2 is out of range.
                //  Continue only if *int_curves2[c_2] is a loop.
                  if (int_curves2[c_2]->is_closed())
                  //  Go forward from segments[0] as long as
                  //    their START's are NOT new_int_pts2_this/other's.
                    int_curves2[c_2]->mk_seg_fw(good2[c_2], 0,
                                                0,
                                                new_int_pts2_this,
                                                num_new_int_pts2_this,
                                                new_int_pts2_other,
                                                num_new_int_pts2_other);
                
                //  Go backward from segments[s_2] as long as
                //    their END' are NOT new_int_pts2_this/other's.
                
                int_curves2[c_2]->mk_seg_bw(good2[c_2], 0,
                                            s_2,
                                            new_int_pts2_this,
                                            num_new_int_pts2_this,
                                            new_int_pts2_other,
                                            num_new_int_pts2_other);
              }
        }
        else  //  if (contains(*int_curves1[i]->start()))
        //  *int_curves1[i] lies entirely INSIDE the trim region.
        //  Unless there exists a corresponding loop on p's domain
        //    make all the segments bad.
        {
          for (j = 0; j < num_int_curves2 && !dir_corr[i][j]; j++)
            ;
          
          if (j == num_int_curves2)
            for (k = 0; k < int_curves1[i]->num_segments; k++)
              good1[i][k] = 0;
        }
      else  //  if (!int_curves1[i]->is_closed())
      //  *int_curves1[i] intersect with *this outside of the trim region.
      //  Mark all the segments bad.
      //  Also, mark all the corresponding segments on p's domain bad.
      {
        for (j = 0; j < int_curves1[i]->num_segments; j++)
          good1[i][j] = 0;
        
        for (j = 0; j < num_new_int_pts1_other; j++)
          if (ind_int_curves1_other[j] == i)
          {
            other_pt = pt_surf_assoc.find_assoc_pt(new_int_pts1_other[j],
                                                   surf, p.surf);
            
            for (k = 0;
                 k < num_new_int_pts2_this
                 &&
                 !new_int_pts2_this[k]->equiv(*other_pt);
                 k++)
              ;
            //  Every new_int_pts2_this[k] is added to some curve, and
            //  Every pt in pt_surf_assoc is also added to some curve.
            //  Thus, it is safe to call equiv instead of equal.
            
            assert(k < num_new_int_pts2_this);
            c_2 = ind_int_curves2_this[k];
            s_2 = int_curves2[c_2]->locate_pt_seg_end(*other_pt);
            
            //  On p's domain, march along the curve to mark bad segments.
            //  Go forward from segments[s_2 + 1] as long as
            //    their START's are NOT new_int_pts2_this/other's, and
            //  go backward from segments[s_2] as long as
            //    their END's are NOT new_int_pts2_this/other's.
            
            if (s_2 < int_curves2[c_2]->num_segments - 1)
            //  Go forward from segments[s_2 + 1] as long as
            //    their START's are NOT new_int_pts2_this/other's.
              int_curves2[c_2]->mk_seg_fw(good2[c_2], 0,
                                          s_2 + 1,
                                          new_int_pts2_this,
                                          num_new_int_pts2_this,
                                          new_int_pts2_other,
                                          num_new_int_pts2_other);
            else  //  if (s_2 == int_curves2[c_2]->num_segments - 1)
            //  s_2 is out of range.
            //  Continue only if *int_curves2[c_2] is a loop.
              if (int_curves2[c_2]->is_closed())
              //  Go forward from segments[0] as long as
              //    their START's are NOT new_int_pts2_this/other's.
                int_curves2[c_2]->mk_seg_fw(good2[c_2], 0,
                                            0,
                                            new_int_pts2_this,
                                            num_new_int_pts2_this,
                                            new_int_pts2_other,
                                            num_new_int_pts2_other);
            
            //  Go backward from segments[s_2] as long as
            //    their END' are NOT new_int_pts2_this/other's.
            
            int_curves2[c_2]->mk_seg_bw(good2[c_2], 0,
                                        s_2,
                                        new_int_pts2_this,
                                        num_new_int_pts2_this,
                                        new_int_pts2_other,
                                        num_new_int_pts2_other);
          }
      }
    else  //  if (num_new_trim_pts1[i] > 0)
    //  if *int_curves1[i] has intersections with the trim_curves
    {
      //  Determine whether we are currently IN/OUTSIDE the trim region.
      
//      cerr << "   clipping: 01: int_curves1[" << i << "]->is_closed() = " << int_curves1[i]->is_closed() << endl << flush;
      
      if (int_curves1[i]->is_closed())
      //  *int_curves1[i] is a loop.
      //  Find a pt which is CERTAINLY IN/OUTSIDE the loop.
      {
        in_on_out_region = 0;
        
        for (j = 0;
             j < int_curves1[i]->num_segments
             &&
             !(in_on_out_region = in_on_out(*int_curves1[i]->start()));
             j++)
          int_curves1[i]->rotate_closed_curve(1);
        
        while (!in_on_out_region)
        {
          int_curves1[i]->start()->shrink(shrink_step, shrink_step);
          
          if (!(in_on_out_region = in_on_out(*int_curves1[i]->start())))
            int_curves1[i]->rotate_closed_curve(1);
        }
      }
      else  //  if (!int_curves1[i]->is_closed())
      //  All the trim_curves are defined inside *this.
      //  Since *int_curves1[i] is not a loop, in particular,
      //        *int_curves1[i]->start is at the boundary of *this,
      //    mark it bad.
        in_on_out_region = - 1;
      
      assert(in_on_out_region);
      
      //  On this domain, march along the curve to mark segments bad.
      
      for (j = 0; j < int_curves1[i]->num_segments; j++)
      {
//        cerr << "   clipping: 02: j = " << j << ", in_on_out_region = " << in_on_out_region << endl << flush;
//        cerr << "   good1[" << i << "] = " << flush;
//        for (unsigned long j1 = 0; j1 < int_curves1[i]->num_segments - 1; j1++)
//          cerr << good1[i][j1] << ", " << flush;
//        cerr << good1[i][int_curves1[i]->num_segments - 1] << endl << flush;
//        cerr << " ------------------------------------- " << endl << flush;
        
        if (in_on_out_region < 0)
          good1[i][j] = 0;
        
//        cerr << "   good1[" << i << "] = " << flush;
//        for (unsigned long j1 = 0; j1 < int_curves1[i]->num_segments - 1; j1++)
//          cerr << good1[i][j1] << ", " << flush;
//        cerr << good1[i][int_curves1[i]->num_segments - 1] << endl << flush;
//        cerr << " ------------------------------------- " << endl << flush;
        
        for (k = 0;
             k < num_new_int_pts1_this
             &&
             !new_int_pts1_this[k]->equiv(*int_curves1[i]->segments[j]->end);
             k++)
          ;
        //  Every new_int_pts1_this[k] is added to some curve,
        //    i.e., it is an endpt of some segment of the curve.
        //  Thus, it is safe to call equiv instead of equal.
        
//        cerr << "   clipping: 03: k = " << k << ", num_new_int_pts1_this = " << num_new_int_pts1_this << endl << flush;
        
        if (k < num_new_int_pts1_this)
        //  if *new_int_pts1_this[k] lies on some trim_curve
        {
          in_on_out_region = - in_on_out_region;
          other_pt         = pt_surf_assoc.find_assoc_pt(new_int_pts1_this[k],
                                                         surf, p.surf);
          
//          cerr << "   clipping: 04: other_pt = " << *other_pt << endl << flush;
          
          //  *new_int_pts1_this[k] lies on some trim_curve on this patch.
          //  Unless there is a degeneracy,
          //    *other_pt  does not lie on any trim_curve on p's patch.
          //  Thus, it is safe to call contains.
          if (p.contains(*other_pt))
          {
            for (l = 0;
                 l < num_new_int_pts2_other
                 &&
                 !new_int_pts2_other[l]->equiv(*other_pt);
                 l++)
              ;
            //  Every new_int_pts2_other[k] is added to some curve, and
            //  Every pt in pt_surf_assoc is also added to some curve.
            //  Thus, it is safe to call equiv instead of equal.
            
            assert(l < num_new_int_pts2_other);
            c_2   = ind_int_curves2_other[l];
            s_2   = int_curves2[c_2]->locate_pt_seg_end(*other_pt);
            is_fw = dir_corr[i][c_2];
//            cerr << "   is_fw = " << is_fw << " = dir_corr[" << i << "][" << c_2 << "]" << endl << flush;
//            cerr << " ------------------------------------- " << endl << flush;
            assert(is_fw);
            
            //  On p's domain, march along the curve to mark bad segments.
            //  If we are currently inside the region and
            //        march to the same direction on both domains then
            //    go forward from segments[s_2 + 1] as long as
            //      their START's are NOT new_int_pts2_this/other's.
            //  If we are currently outside the region and
            //     march to the opposite directions then
            //    go forward from segments[s_2 + 1] as long as
            //      their START's are NOT new_int_pts2_this/other's.
            //  If we are currently inside the region and
            //     march to the same direction on both domains then
            //    go backward from segments[s_2] as long as
            //      their END's are NOT new_int_pts2_this/other's.
            //  If we are currently outside the region and
            //     march to the opposite directions then
            //    go backward from segments[s_2] as long as
            //      their END's are NOT new_int_pts2_this/other's.
            
            if (in_on_out_region * is_fw > 0)
            //  Go backward from segments[s_2] as long as
            //     their END's are NOT new_int_pts2_this/other's.
              int_curves2[c_2]->mk_seg_bw(good2[c_2], 0,
                                          s_2,
                                          new_int_pts2_this,
                                          num_new_int_pts2_this,
                                          new_int_pts2_other,
                                          num_new_int_pts2_other);
            else  //  if (in_on_out_region * is_fw < 0)
            //  Go forward from segments[s_2 + 1] as long as
            //     their START's are NOT new_int_pts2_this/other's.
              int_curves2[c_2]->mk_seg_fw(good2[c_2], 0,
                                          s_2 + 1,
                                          new_int_pts2_this,
                                          num_new_int_pts2_this,
                                          new_int_pts2_other,
                                          num_new_int_pts2_other);
          }
          else  //  if (!p.contains(*other_pt))
          {
            //  On this domain, march along the curve to mark bad segments.
            //  Go forward from segments[j + 1] as long as
            //    their START's are NOT new_int_pts1_this/other's, and
            //  Go backward from segments[j] as long as
            //    their END's are NOT new_int_pts1_this/other's.
            
//            cerr << "   clipping: 05: j = " << j << ", int_curves1[" << i << "]->num_segments - 1 = " << int_curves1[i]->num_segments - 1 << endl << flush;
            
            if (j < int_curves1[i]->num_segments - 1)
            //  Go forward from segments[j + 1] as long as
            //    their START's are NOT new_int_pts1_this/other's, and
            {
//              cerr << "   good1[" << i << "] = " << flush;
//              for (unsigned long j1 = 0; j1 < int_curves1[i]->num_segments - 1; j1++)
//                cerr << good1[i][j1] << ", " << flush;
//              cerr << good1[i][int_curves1[i]->num_segments - 1] << endl << flush;
//              cerr << " ------------------------------------- " << endl << flush;
              
              int_curves1[i]->mk_seg_fw(good1[i], 0,
                                        j + 1,
                                        new_int_pts1_this,
                                        num_new_int_pts1_this,
                                        new_int_pts1_other,
                                        num_new_int_pts1_other);
              
//              cerr << "   good1[" << i << "] = " << flush;
//              for (unsigned long j1 = 0; j1 < int_curves1[i]->num_segments - 1; j1++)
//                cerr << good1[i][j1] << ", " << flush;
//              cerr << good1[i][int_curves1[i]->num_segments - 1] << endl << flush;
//              cerr << " ------------------------------------- " << endl << flush;
            }
            
            //  Go backward from segments[j] as long as
            //    their END's are NOT new_int_pts1_this/other's.
            
//            cerr << "   good1[" << i << "] = " << flush;
//            for (unsigned long j1 = 0; j1 < int_curves1[i]->num_segments - 1; j1++)
//              cerr << good1[i][j1] << ", " << flush;
//            cerr << good1[i][int_curves1[i]->num_segments - 1] << endl << flush;
//            cerr << " ------------------------------------- " << endl << flush;
            
            int_curves1[i]->mk_seg_bw(good1[i], 0,
                                      j,
                                      new_int_pts1_this,
                                      num_new_int_pts1_this,
                                      new_int_pts1_other,
                                      num_new_int_pts1_other);
            
//            cerr << "   good1[" << i << "] = " << flush;
//            for (unsigned long j1 = 0; j1 < int_curves1[i]->num_segments - 1; j1++)
//              cerr << good1[i][j1] << ", " << flush;
//            cerr << good1[i][int_curves1[i]->num_segments - 1] << endl << flush;
//            cerr << " ------------------------------------- " << endl << flush;
          }
        }
      }
    }
  }
  
//  cerr << " Debug4.5: kpatch: intersect: ---------- " << endl << flush;
  
  for (i = 0; i < num_int_curves2; i++)
  {
    //  See if there exists a curve lying entirely in/outside the trim region,
    //    i.e., a curve having no intersection with the trim_curves.
    //  If such a curve exists and it also lies outside the trim region then
    //    mark all its segments bad.
    
//    cerr << "   clipping: 10: i = " << i << ", num_int_curves2 = " << num_int_curves2 << ", num_new_trim_pts2[" << i << "] = " << num_new_trim_pts2[i] << endl << flush;
    
    if (!num_new_trim_pts2[i])
    //  if *int_curves2[i] has no intersection with the trim_curves
      if (int_curves2[i]->is_closed())
      //  *int_curves2[i] is a loop.
      //  See if the loop lies entirely in the trim region;
      //  Pick some pts and see if each of them is in/outside the trim region.
      //  We knew that *int_curves2[i] does not hit the trim_curves.
      //  Thus, it is safe to call contains.
        if (!p.contains(*int_curves2[i]->start()))
        //  *int_curves2[i] lies entirely outside the trim region.
        //  Mark all the segments bad and
        //    remove the corresponding loop on this domain.
        {
          for (j = 0; j < int_curves2[i]->num_segments; j++)
            good2[i][j] = 0;
          
          if (!num_new_int_pts2[i])
          {
            for (j = 0; j < num_int_curves1; j++)
              if (dir_corr[j][i])
              //  *int_curves1[i] and *int_curves2[j] correspond and
              //    do not have intersections corresponding with each other.
              //  Mark all the segments of *int_curves2[j] bad.
                for (k = 0; k < int_curves1[j]->num_segments; k++)
                  good1[j][k] = 0;
          }
          else  //  if (num_new_int_pts2[i] > 0)
            for (j = 0; j < num_new_int_pts2_other; j++)
              //  *new_int_pts2_other[j] lies on *int_curves2[i] and
              //    corresponds to some pt, *other_pt, on this domain.
              //  Find a curve on this domain on which *other_pt lies.
              if (ind_int_curves2_other[j] == i)
              {
                other_pt = pt_surf_assoc.find_assoc_pt(new_int_pts2_other[j],
                                                       p.surf, surf);
                
                for (k = 0;
                     k < num_new_int_pts1_this
                     &&
                     !new_int_pts1_this[k]->equiv(*other_pt);
                     k++)
                  ;
                //  Every new_int_pts1_this[k] is added to some curve, and
                //  Every pt in pt_surf_assoc is also added to some curve.
                //  Thus, it is safe to call equiv instead of equal.
                
                assert(k < num_new_int_pts1_this);
                c_1 = ind_int_curves1_this[k];
                s_1 = int_curves1[c_1]->locate_pt_seg_end(*other_pt);
                
                //  On this domain, march along the curve to mark bad segments.
                //  Go forward from segments[s_1 + 1] as long as
                //    their START's are NOT new_int_pts1_this/other's, and
                //  go backward from segments[s_1] as long as
                //    their END's are NOT new_int_pts1_this/other's.
                
                if (s_1 < int_curves1[c_1]->num_segments - 1)
                //  Go forward from segments[s_1 + 1] as long as
                //    their START's are NOT new_int_pts1_this/other's.
                  int_curves1[c_1]->mk_seg_fw(good1[c_1], 0,
                                              s_1 + 1,
                                              new_int_pts1_this,
                                              num_new_int_pts1_this,
                                              new_int_pts1_other,
                                              num_new_int_pts1_other);
                else  //  if (s_1 == int_curves1[c_1]->num_segments - 1)
                //  s_1 is out of range.
                //  Continue only if *int_curves1[c_1] is a loop.
                  if (int_curves1[c_1]->is_closed())
                  //  Go forward from segments[0] as long as
                  //    their START's are NOT new_int_pts1_this/other's.
                    int_curves1[c_1]->mk_seg_fw(good1[c_1], 0,
                                                0,
                                                new_int_pts1_this,
                                                num_new_int_pts1_this,
                                                new_int_pts1_other,
                                                num_new_int_pts1_other);
                
                //  Go backward from segments[s_1] as long as
                //    their END' are NOT new_int_pts1_this/other's.
                
                int_curves1[c_1]->mk_seg_bw(good1[c_1], 0,
                                            s_1,
                                            new_int_pts1_this,
                                            num_new_int_pts1_this,
                                            new_int_pts1_other,
                                            num_new_int_pts1_other);
              }
        }
        else  //  if (p.contains(*int_curves2[i]->start()))
        //  *int_curves2[i] lies entirely INSIDE the trim region.
        //  Unless there exists a corresponding loop on this domain
        //    make all the segments bad.
        {
          for (j = 0; j < num_int_curves1 && !dir_corr[j][i]; j++)
            ;
          
          if (j == num_int_curves1)
            for (k = 0; k < int_curves2[i]->num_segments; k++)
              good2[i][k] = 0;
        }
      else  //  if (!int_curves2[i]->is_closed())
      //  *int_curves2[i] intersect with *this outside of the trim region.
      //  Mark all the segments bad.
      //  Also, mark all the corresponding segments on this domain bad.
      {
        for (j = 0; j < int_curves2[i]->num_segments; j++)
          good2[i][j] = 0;
        
        for (j = 0; j < num_new_int_pts2_other; j++)
          if (ind_int_curves2_other[j] == i)
          {
            other_pt = pt_surf_assoc.find_assoc_pt(new_int_pts2_other[j],
                                                   p.surf, surf);
            
            for (k = 0;
                 k < num_new_int_pts1_this
                 &&
                 !new_int_pts1_this[k]->equiv(*other_pt);
                 k++)
              ;
            //  Every new_int_pts1_this[k] is added to some curve, and
            //  Every pt in pt_surf_assoc is also added to some curve.
            //  Thus, it is safe to call equiv instead of equal.
            
            assert(k < num_new_int_pts1_this);
            c_1 = ind_int_curves1_this[k];
            s_1 = int_curves1[c_1]->locate_pt_seg_end(*other_pt);
            
            //  On this domain, march along the curve to mark bad segments.
            //  Go forward from segments[s_1 + 1] as long as
            //    their START's are NOT new_int_pts1_this/other's, and
            //  go backward from segments[s_1] as long as
            //    their END's are NOT new_int_pts1_this/other's.
            
            if (s_1 < int_curves1[c_1]->num_segments - 1)
            //  Go forward from segments[s_1 + 1] as long as
            //    their START's are NOT new_int_pts1_this/other's.
              int_curves1[c_1]->mk_seg_fw(good1[c_1], 0,
                                          s_1 + 1,
                                          new_int_pts1_this,
                                          num_new_int_pts1_this,
                                          new_int_pts1_other,
                                          num_new_int_pts1_other);
            else  //  if (s_1 == int_curves1[c_1]->num_segments - 1)
            //  s_1 is out of range.
            //  Continue only if *int_curves1[c_1] is a loop.
              if (int_curves1[c_1]->is_closed())
              //  Go forward from segments[0] as long as
              //    their START's are NOT new_int_pts1_this/other's.
                int_curves1[c_1]->mk_seg_fw(good1[c_1], 0,
                                            0,
                                            new_int_pts1_this,
                                            num_new_int_pts1_this,
                                            new_int_pts1_other,
                                            num_new_int_pts1_other);
            
            //  Go backward from segments[s_1] as long as
            //    their END' are NOT new_int_pts1_this/other's.
            
            int_curves1[c_1]->mk_seg_bw(good1[c_1], 0,
                                        s_1,
                                        new_int_pts1_this,
                                        num_new_int_pts1_this,
                                        new_int_pts1_other,
                                        num_new_int_pts1_other);
          }
      }
    else  //  if (num_new_trim_pts2[i] > 0)
    //  if *int_curves2[i] has intersections with the trim_curves
    {
//      cerr << "   clipping: 11: int_curves2[" << i << "]->is_closed() = " << int_curves2[i]->is_closed() << endl << flush;
      
      if (int_curves2[i]->is_closed())
      //  *int_curves2[i] is a loop.
      //  Find a pt which is CERTAINLY IN/OUTSIDE the loop.
      {
        in_on_out_region = 0;
        
//        cerr << "   good2[" << i << "] = " << flush;
//        for (unsigned long j1 = 0; j1 < int_curves2[i]->num_segments - 1; j1++)
//          cerr << good2[i][j1] << ", " << flush;
//        cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//        cerr << " ------------------------------------- " << endl << flush;
        
        for (j = 0;
             j < int_curves2[i]->num_segments
             &&
             !(in_on_out_region = p.in_on_out(*int_curves2[i]->start()));
             j++)
        {
          int_curves2[i]->rotate_closed_curve(1);
          
          //  Some entries of good2 may have been 0.
          //  Thus, we must rotate good2.
          
          g = good2[i][0];
          
          for (k = 1; k < int_curves2[i]->num_segments; k++)
            good2[i][k - 1] = good2[i][k];
          
          good2[i][int_curves2[i]->num_segments - 1] = g;
          
//          cerr << "   good2[" << i << "] = " << flush;
//          for (unsigned long j1 = 0; j1 < int_curves2[i]->num_segments - 1; j1++)
//            cerr << good2[i][j1] << ", " << flush;
//          cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//          cerr << " ------------------------------------- " << endl << flush;
        }
        
        while (!in_on_out_region)
        {
          int_curves2[i]->start()->shrink(shrink_step, shrink_step);
          
          if (!(in_on_out_region = in_on_out(*int_curves2[i]->start())))
          {
            int_curves2[i]->rotate_closed_curve(1);
            
            //  Some entries of good2 may have been 0.
            //  Thus, we must rotate good2.
            
            g = good2[i][0];
            
            for (k = 1; k < int_curves2[i]->num_segments; j++)
              good2[i][k - 1] = good2[i][k];
            
            good2[i][int_curves2[i]->num_segments - 1] = g;
            
//            cerr << "   good2[" << i << "] = " << flush;
//            for (unsigned long j1 = 0; j1 < int_curves2[i]->num_segments - 1; j1++)
//              cerr << good2[i][j1] << ", " << flush;
//            cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//            cerr << " ------------------------------------- " << endl << flush;
          }
        }
      }
      else  //  if (!int_curves2[i]->is_closed())
      //  All the trim_curves are defined inside p.
      //  Since *int_curves2[i] is not a loop, in particular,
      //        *int_curves2[i]->start is at the boundary of *this,
      //    mark it bad.
        in_on_out_region = - 1;
      
      assert(in_on_out_region);
      
      //  On p's domain, march along the curve to mark segments bad.
      
      for (j = 0; j < int_curves2[i]->num_segments; j++)
      {
//        cerr << "   clipping: 12: j = " << j << ", int_curves2[" << i << "]->num_segments = " << int_curves2[i]->num_segments << endl << flush;
//        cerr << "   good2[" << i << "] = " << flush;
//        for (unsigned long j1 = 0; j1 < int_curves2[i]->num_segments - 1; j1++)
//          cerr << good2[i][j1] << ", " << flush;
//        cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//        cerr << " ------------------------------------- " << endl << flush;
        
        if (in_on_out_region < 0)
          good2[i][j] = 0;
        
//        cerr << "   good2[" << i << "] = " << flush;
//        for (unsigned long j1 = 0; j1 < int_curves2[i]->num_segments - 1; j1++)
//          cerr << good2[i][j1] << ", " << flush;
//        cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//        cerr << " ------------------------------------- " << endl << flush;
        
        for (k = 0;
             k < num_new_int_pts2_this
             &&
             !new_int_pts2_this[k]->equiv(*int_curves2[i]->segments[j]->end);
             k++)
          ;
        //  Every new_int_pts2_this[k] is added to some curve,
        //    i.e., it is an endpt of some segment of the curve.
        //  Thus, it is safe to call equiv instead of equal.
        
//        cerr << "   clipping: 13: k = " << k << ", num_new_int_pts2_this = " << num_new_int_pts2_this << endl << flush;
        
        if (k < num_new_int_pts2_this)
        //  if *new_int_pts2_this[k] lies on some trim_curve
        {
          in_on_out_region = - in_on_out_region;
          other_pt         = pt_surf_assoc.find_assoc_pt(new_int_pts2_this[k],
                                                         p.surf, surf);
          
//          cerr << "   clipping: 14: other_pt = " << *other_pt << endl << flush;
          
          //  *new_int_pts2_this[k] lies on some trim_curve on p's patch.
          //  Unless there is a degeneracy,
          //    *other_pt  does not lie on any trim_curve on this patch.
          //  Thus, it is safe to call contains.
          if (contains(*other_pt))
          {
            for (l = 0;
                 l < num_new_int_pts1_other
                 &&
                 !new_int_pts1_other[l]->equiv(*other_pt);
                 l++)
              ;
            //  Every new_int_pts1_other[k] is added to some curve, and
            //  Every pt in pt_surf_assoc is also added to some curve.
            //  Thus, it is safe to call equiv instead of equal.
            
            assert(l < num_new_int_pts1_other);
            c_1   = ind_int_curves1_other[l];
            s_1   = int_curves1[c_1]->locate_pt_seg_end(*other_pt);
            is_fw = dir_corr[c_1][i];
            
//            cerr << "   is_fw = " << is_fw << " = dir_corr[" << c_1 << "][" << i << "]" << endl << flush;
//            cerr << " ------------------------------------- " << endl << flush;
            assert(is_fw);
            
            //  On this domain, march along the curve to mark bad segments.
            //  If we are currently inside the region and
            //        march to the same direction on both domains then
            //    go forward from segments[s_1 + 1] as long as
            //      their START's are NOT new_int_pts1_this/other's.
            //  If we are currently outside the region and
            //     march to the opposite directions then
            //    go forward from segments[s_1 + 1] as long as
            //      their START's are NOT new_int_pts1_this/other's.
            //  If we are currently inside the region and
            //     march to the same direction on both domains then
            //    go backward from segments[s_1] as long as
            //      their END's are NOT new_int_pts1_this/other's.
            //  If we are currently outside the region and
            //     march to the opposite directions then
            //    go backward from segments[s_1] as long as
            //      their END's are NOT new_int_pts1_this/other's.
            
            if (in_on_out_region * is_fw > 0)
            //  Go backward from segments[s_1] as long as
            //     their END's are NOT new_int_pts1_this/other's.
              int_curves1[c_1]->mk_seg_bw(good1[c_1], 0,
                                          s_1,
                                          new_int_pts1_this,
                                          num_new_int_pts1_this,
                                          new_int_pts1_other,
                                          num_new_int_pts1_other);
            else  //  if (in_on_out_region * is_fw < 0)
            //  Go forward from segments[s_1 + 1] as long as
            //     their START's are NOT new_int_pts2_this/other's.
              int_curves1[c_1]->mk_seg_fw(good1[c_1], 0,
                                          s_1 + 1,
                                          new_int_pts1_this,
                                          num_new_int_pts1_this,
                                          new_int_pts1_other,
                                          num_new_int_pts1_other);
          }
          else  //  if (!contains(*other_pt))
          {
            //  On this domain, march along the curve to mark bad segments.
            //  Go forward from segments[j + 1] as long as
            //    their START's are NOT new_int_pts2_this/other's, and
            //  Go backward from segments[j] as long as
            //    their END's are NOT new_int_pts2_this/other's.
            
//            cerr << "   clipping: 15: j = " << j << ", int_curves2[" << i << "]->num_segments - 1 = " << int_curves2[i]->num_segments - 1 << endl << flush;
            
            if (j < int_curves2[i]->num_segments - 1)
            //  Go forward from segments[j + 1] as long as
            //    their START's are NOT new_int_pts1_this/other's, and
            {
//              cerr << "   good2[" << i << "] = " << flush;
//              for (unsigned long j1 = 0; j1 < int_curves2[i]->num_segments - 1; j1++)
//                cerr << good2[i][j1] << ", " << flush;
//              cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//              cerr << " ------------------------------------- " << endl << flush;
              
              int_curves2[i]->mk_seg_fw(good2[i], 0,
                                        j + 1,
                                        new_int_pts2_this,
                                        num_new_int_pts2_this,
                                        new_int_pts2_other,
                                        num_new_int_pts2_other);
              
//              cerr << "   good2[" << i << "] = " << flush;
//              for (unsigned long j1 = 0; j1 < int_curves2[i]->num_segments - 1; j1++)
//                cerr << good2[i][j1] << ", " << flush;
//              cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//              cerr << " ------------------------------------- " << endl << flush;
            }
            
            //  Go backward from segments[j] as long as
            //    their END's are NOT new_int_pts2_this/other's.
            
            int_curves2[i]->mk_seg_bw(good2[i], 0,
                                      j,
                                      new_int_pts2_this,
                                      num_new_int_pts2_this,
                                      new_int_pts2_other,
                                      num_new_int_pts2_other);
          }
        }
      }
    }
  }
  
  //  6. List new intersection curves.
  
  int           in_curve, good_closed;
  long          head1, head2, tail1, tail2, prev_head2;
  int           f_3;
  K_CURVE**     new_int_curves1;
  K_CURVE**     new_int_curves2;
  unsigned long num_new_int_curves;
  
//  cerr << " Debug5: kpatch: intersect: ---------- " << endl << flush;
  
//  for (i = 0; i < num_int_curves1; i++)
//  {
//    cerr << "   *int_curves1[" << i << "]->poly = " << endl << *int_curves1[i]->poly << flush;
//    cerr << "   *int_curves1[" << i << "] = " << endl << *int_curves1[i] << endl << flush;
//    cerr << "   good1[" << i << "] = " << flush;
//    for (j = 0; j < int_curves1[i]->num_segments - 1; j++)
//      cerr << good1[i][j] << ", " << flush;
//    cerr << good1[i][int_curves1[i]->num_segments - 1] << endl << flush;
//    cerr << " ------------------------------------- " << endl << flush;
//  }
//  for (i = 0; i < num_int_curves2; i++)
//  {
//    cerr << "   *int_curves2[" << i << "]->poly = " << endl << *int_curves2[i]->poly << flush;
//    cerr << "   *int_curves2[" << i << "] = " << endl << *int_curves2[i] << endl << flush;
//    cerr << "   good2[" << i << "] = " << flush;
//    for (j = 0; j < int_curves2[i]->num_segments - 1; j++)
//      cerr << good2[i][j] << ", " << flush;
//    cerr << good2[i][int_curves2[i]->num_segments - 1] << endl << flush;
//    cerr << " ------------------------------------- " << endl << flush;
//  }
  
  // Determine how may new int_curves there are.
  
  num_new_int_curves = 0;
  
  for (i = 0; i < num_int_curves1; i++)
  {
    good_closed = 0;
    
    if (int_curves1[i]->is_closed())
    {
      in_curve = good1[i][int_curves1[i]->num_segments - 1];
      
      if (in_curve)
        good_closed = 1;
    }
    else  //  if (!int_curves1[i]->is_closed())
    //  Since *int_curves1[i] is not a loop, in particular,
    //        *int_curves1[i]->start is at the boundary of *this,
    //    we are out.
      in_curve = 0;
    
    for (j = 0; j < int_curves1[i]->num_segments; j++)
      if (!in_curve && good1[i][j])
      {
        in_curve = 1;
        num_new_int_curves++;
      }
      else if (in_curve && !good1[i][j])
      {
        in_curve    = 0;
        good_closed = 0;
      }
    
    if (good_closed)
    //  *int_curves1[i] is a loop and al its segments are marked good.
      num_new_int_curves++;
  }
  
//  cerr << "   num_new_int_curves = " << num_new_int_curves << endl << flush;
//  cerr << " ------------------------------------- " << endl << flush;
  
  //  Generate new_int_curves on both this and p's domains.
  
  new_int_curves1 = new K_CURVE* [num_new_int_curves];
  new_int_curves2 = new K_CURVE* [num_new_int_curves];
  i               = 0;
  
  for (j = 0; j < num_int_curves1; j++)
  {
    in_curve = 0;
    
    if (good1[j][0])
    {
      //  See whether or not *int_curves1[j] is a good loop.
      
      for (tail1 = int_curves1[j]->num_segments;
           tail1 > 0 && good1[j][tail1 - 1];
           tail1--)
        ;
      
      if (tail1 > 0)
      {
        head1       = tail1 % int_curves1[j]->num_segments;
        good_closed = 0;
      }
      else  //  if (!tail1)
        good_closed = 1;
      
      in_curve = 1;
    }
    else  //  if (!good1[j][0])
    {
      tail1       = int_curves1[j]->num_segments;
      in_curve    = 0;
      good_closed = 0;
    }
    
    if (good_closed)
    //  *int_curves1[j] is a loop and all its segments are marked good.
    {
      for (k = 0; k < num_int_curves2 && !dir_corr[j][k]; k++)
        ;
      
      assert(k < num_int_curves2);
      
      new_int_curves1[i] = new K_CURVE(*int_curves1[j]);
      new_int_curves1[i]->ref_count++;
      new_int_curves2[i] = new K_CURVE(*int_curves2[k]);
      new_int_curves2[i]->ref_count++;
      new_int_curves1[i]->assoc(new_int_curves2[i], dir_corr[j][k]);
      i++;
    }
    else  //  if (!good_closed)
    {
      for (k = 1; k < tail1; k++)
        if (!in_curve && good1[j][k])
        {
          in_curve = 1;
          head1    = k;
        }
        else if (in_curve && !good1[j][k])
        {
          in_curve = 0;
          
          //  head1 has already been set.
          //  k is incremented at least once until an execution reaches here.
          //  Thus, head1 < k.
          
          new_int_curves1[i] = new K_CURVE(*int_curves1[j], head1, k - 1);
          new_int_curves1[i]->ref_count++;
          
          //  Find *new_int_curevs2[i] on p's domain.
          
          f_3 = 0;
          
          for (l = 0; !f_3 && l < num_int_curves2; l++)
          {
//            cerr << "   good2[" << l << "] = " << flush;
//            for (unsigned long l1 = 0; l1 < int_curves2[l]->num_segments - 1; l1++)
//              cerr << good2[l][l1] << ", " << flush;
//            cerr << good2[l][int_curves2[l]->num_segments - 1] << endl << flush;
//            cerr << " ------------------------------------- " << endl << flush;
            
            for (head2 = 0, prev_head2 = int_curves2[l]->num_segments - 1;
                 !f_3 && head2 < int_curves2[l]->num_segments;
                 prev_head2 = ++head2 - 1)
            {
              if (good2[l][head2] && !good2[l][prev_head2])
              {
                other_pt =
                  pt_surf_assoc.
                  find_assoc_pt(int_curves2[l]->segments[head2]->start,
                                p.surf, surf);
                
                if (other_pt->equiv(*int_curves1[j]->segments[head1]->start)
                    ||
                    other_pt->equiv(*int_curves1[j]->segments[k - 1]->end))
                  f_3 = 1;
//                else
//                  cerr << "   new_int_curves: good2 but bad1. " << endl << flush;
              }
              
              if (f_3)
              {
                for (tail2 = (head2 + 1) % int_curves2[l]->num_segments;
                     tail2 != head2 && good2[l][tail2];
                     tail2 = (tail2 + 1) % int_curves2[l]->num_segments)
                  ;
                
                if (!tail2)
                  tail2 = int_curves2[l]->num_segments;
                
                //  head2 < tail2
                
                new_int_curves2[i] =
                  new K_CURVE(*int_curves2[l], head2, tail2 - 1);
                new_int_curves2[i]->ref_count++;
                new_int_curves1[i]->assoc(new_int_curves2[i], dir_corr[j][l]);
                i++;
              }
            }
          }
          
          assert(f_3);
        }
      
      if (in_curve)
      //  *int_curves1[j] is a loop and
      //    not all of its segments are marked good.
      {
        new_int_curves1[i] =
          new K_CURVE(*int_curves1[j], head1, tail1 - 1);
        new_int_curves1[i]->ref_count++;
        
        //  Find *new_int_curevs2[i] on p's domain.
        
        f_3 = 0;
        
        for (l = 0; !f_3 && l < num_int_curves2; l++)
        {
//          cerr << "   good2[" << l << "] = " << flush;
//          for (unsigned long l1 = 0; l1 < int_curves2[l]->num_segments - 1; l1++)
//            cerr << good2[l][l1] << ", " << flush;
//          cerr << good2[l][int_curves2[l]->num_segments - 1] << endl << flush;
//          cerr << " ------------------------------------- " << endl << flush;
          
          for (head2 = 0, prev_head2 = int_curves2[l]->num_segments - 1;
               !f_3 && head2 < int_curves2[l]->num_segments;
               prev_head2 = ++head2 - 1)
          {
            if (good2[l][head2] && !good2[l][prev_head2])
            {
              other_pt =
                pt_surf_assoc.
                find_assoc_pt(int_curves2[l]->segments[head2]->start,
                              p.surf, surf);
              
              if (other_pt->equiv(*int_curves1[j]->segments[head1]->start)
                  ||
                  other_pt->equiv(*int_curves1[j]->segments[tail1 - 1]->end))
                f_3 = 1;
//              else
//                cerr << "   new_int_curves: good2 but bad1. " << endl << flush;
            }
            
            if (f_3)
            {
              for (tail2 = (head2 + 1) % int_curves2[l]->num_segments;
                   tail2 != head2 && good2[l][tail2];
                   tail2 = (tail2 + 1) % int_curves2[l]->num_segments)
                ;
              
              if (!tail2)
                tail2 = int_curves2[l]->num_segments;
              
              new_int_curves2[i] =
                new K_CURVE(*int_curves2[l], head2, tail2 - 1);
              new_int_curves2[i]->ref_count++;
              new_int_curves1[i]->assoc(new_int_curves2[i], dir_corr[j][l]);
              i++;
            }
          }
        }
        
        assert(f_3);
      }
    }
  }
  
//  cerr << "   ID = " << ID << ", p.ID = " << p.ID << endl << flush;
//  for (i = 0; i < num_new_int_curves; i++)
//  {
//    cerr << "   *new_int_curves1[" << i << "]->poly = " << endl << *new_int_curves1[i]->poly << flush;
//    cerr << "   *new_int_curves1[" << i << "] = " << endl << *new_int_curves1[i] << endl << flush;
//  }
//  for (i = 0; i < num_new_int_curves; i++)
//  {
//    cerr << "   *new_int_curves2[" << i << "]->poly = " << endl << *new_int_curves2[i]->poly << flush;
//    cerr << "   *new_int_curves2[" << i << "] = " << endl << *new_int_curves2[i] << endl << flush;
//  }
  
  //  Add new intersection curves to patches.
  
  unsigned long num_int_curves_proto;
  K_CURVE**     int_curves_proto;
  K_SURF**      adj_int_surfs_proto;
  K_PATCH**     adj_int_patches_proto;
  
  if (num_new_int_curves > 0)
  {
    num_int_curves_proto  = num_int_curves + num_new_int_curves;
    int_curves_proto      = new K_CURVE* [num_int_curves_proto];
    adj_int_surfs_proto   = new K_SURF* [num_int_curves_proto];
    adj_int_patches_proto = new K_PATCH* [num_int_curves_proto];
    
    for (i = 0; i < num_int_curves; i++)
    {
      int_curves_proto[i]      = int_curves[i];
      int_curves_proto[i]->ref_count++;
      adj_int_surfs_proto[i]   = adj_int_surfs[i];
      adj_int_surfs_proto[i]->ref_count++;
      adj_int_patches_proto[i] = adj_int_patches[i];
      adj_int_patches_proto[i]->ref_count++;
    }
    
    for (i = 0; i < num_new_int_curves; i++)
    {
      int_curves_proto[i + num_int_curves]      = new_int_curves1[i];
      int_curves_proto[i + num_int_curves]->ref_count++;
      adj_int_surfs_proto[i + num_int_curves]   = p.surf;
      adj_int_surfs_proto[i + num_int_curves]->ref_count++;
      adj_int_patches_proto[i + num_int_curves] = &p;
      adj_int_patches_proto[i + num_int_curves]->ref_count++;
    }
    
    for (i = 0; i < num_int_curves; i++)
    {
      if (!--int_curves[i]->ref_count)
        delete int_curves[i];
      
      if (!--adj_int_surfs[i]->ref_count)
        delete adj_int_surfs[i];
      
      if (!--adj_int_patches[i]->ref_count)
        delete adj_int_patches[i];
    }
    
    delete [] int_curves;
    delete [] adj_int_surfs;
    delete [] adj_int_patches;
    
    num_int_curves  = num_int_curves_proto;
    int_curves      = int_curves_proto;
    adj_int_surfs   = adj_int_surfs_proto;
    adj_int_patches = adj_int_patches_proto;
  }
  
  if (num_new_int_curves > 0)
  {
    num_int_curves_proto  = p.num_int_curves + num_new_int_curves;
    int_curves_proto      = new K_CURVE* [num_int_curves_proto];
    adj_int_surfs_proto   = new K_SURF* [num_int_curves_proto];
    adj_int_patches_proto = new K_PATCH* [num_int_curves_proto];
    
    for (i = 0; i < p.num_int_curves; i++)
    {
      int_curves_proto[i]      = p.int_curves[i];
      int_curves_proto[i]->ref_count++;
      adj_int_surfs_proto[i]   = p.adj_int_surfs[i];
      adj_int_surfs_proto[i]->ref_count++;
      adj_int_patches_proto[i] = p.adj_int_patches[i];
      adj_int_patches_proto[i]->ref_count++;
    }
    
    for (i = 0; i < num_new_int_curves; i++)
    {
      int_curves_proto[i + p.num_int_curves]      = new_int_curves2[i];
      int_curves_proto[i + p.num_int_curves]->ref_count++;
      adj_int_surfs_proto[i + p.num_int_curves]   = surf;
      adj_int_surfs_proto[i + p.num_int_curves]->ref_count++;
      adj_int_patches_proto[i + p.num_int_curves] = this;
      adj_int_patches_proto[i + p.num_int_curves]->ref_count++;
    }
    
    for (i = 0; i < p.num_int_curves; i++)
    {
      if (!--p.int_curves[i]->ref_count)
        delete p.int_curves[i];
      
      if (!--p.adj_int_surfs[i]->ref_count)
        delete p.adj_int_surfs[i];
      
      if (!--p.adj_int_patches[i]->ref_count)
        delete p.adj_int_patches[i];
    }
    
    delete [] p.int_curves;
    delete [] p.adj_int_surfs;
    delete [] p.adj_int_patches;
    
    p.num_int_curves  = num_int_curves_proto;
    p.int_curves      = int_curves_proto;
    p.adj_int_surfs   = adj_int_surfs_proto;
    p.adj_int_patches = adj_int_patches_proto;
  }
  
  if (num_new_int_curves > 0)
  {
    for (i = 0; i < num_new_int_curves; i++)
      if (!--new_int_curves1[i]->ref_count)
        delete new_int_curves1[i];
    
    delete [] new_int_curves1;
    
    for (i = 0; i < num_new_int_curves; i++)
      if (!--new_int_curves2[i]->ref_count)
        delete new_int_curves2[i];
    
    delete [] new_int_curves2;
  }
  
  if (num_int_curves1 > 0)
  {
    for (i = 0; i < num_int_curves1; i++)
      delete [] good1[i];
    
    delete [] good1;
  }
  
  if (num_int_curves2 > 0)
  {
    for (i = 0; i < num_int_curves2; i++)
      delete [] good2[i];
    
    delete [] good2;
  }
  
  if (num_int_curves1 > 0)
  {
    for (i = 0; i < num_int_curves1; i++)
      if (num_int_curves2 > 0)
        delete [] dir_corr[i];
    
    delete [] dir_corr;
  }
  
  for (i = 0; i < num_new_int_pts1_this; i++)
    if (!--new_int_pts1_this[i]->ref_count)
      delete new_int_pts1_this[i];
  
  delete [] new_int_pts1_this;
  delete [] new_int_pts1_other;
  delete [] ind_int_curves1_this;
  delete [] ind_int_curves1_other;
  
  for (i = 0; i < num_new_int_pts2_this; i++)
    if (!--new_int_pts2_this[i]->ref_count)
      delete new_int_pts2_this[i];
  
  delete [] new_int_pts2_this;
  delete [] new_int_pts2_other;
  delete [] ind_int_curves2_this;
  delete [] ind_int_curves2_other;
  
  if (num_int_curves1 > 0)
  {
    delete [] num_new_int_pts1;
    delete [] num_new_trim_pts1;
  }
  
  if (num_int_curves2 > 0)
  {
    delete [] num_new_int_pts2;
    delete [] num_new_trim_pts2;
  }
  
  if (num_int_curves1 > 0)
  {
    for (i = 0; i < num_int_curves1; i++)
      if (!--int_curves1[i]->ref_count)
        delete int_curves1[i];
    
    delete [] int_curves1;
  }
  
  if (num_int_curves2 > 0)
  {
    for (i = 0; i < num_int_curves2; i++)
      if (!--int_curves2[i]->ref_count)
        delete int_curves2[i];
    
    delete [] int_curves2;
  }
  
  return 0;
}

int K_PATCH :: split_trim_curves(const unsigned long t, const unsigned long s)
{
  assert(t < num_trim_curves);               //  num_trim_curves > t >= 0
  assert(s < trim_curves[t]->num_segments);
  
  unsigned long i, j, k;
  
  if (s < trim_curves[t]->num_segments - 1)
  {
    int           other_dir;
    unsigned long t2, s2;
    K_PATCH*      other_patch;
    K_CURVE**     new_this_trim_curves;
    K_SURF**      new_this_adj_surfs;
    K_PATCH**     new_this_adj_patches;
    K_CURVE**     new_other_trim_curves;
    K_SURF**      new_other_adj_surfs;
    K_PATCH**     new_other_adj_patches;
    
    other_dir       = trim_curves[t]->dir_in_other_dom;
    other_patch     = adj_patches[t];
    
    new_this_trim_curves = new K_CURVE* [num_trim_curves + 1];
    new_this_adj_surfs   = new K_SURF* [num_trim_curves + 1];
    new_this_adj_patches = new K_PATCH* [num_trim_curves + 1];
    
    if (other_dir)
    {
      new_other_trim_curves = new K_CURVE* [other_patch->num_trim_curves + 1];
      new_other_adj_surfs   = new K_SURF* [other_patch->num_trim_curves + 1];
      new_other_adj_patches = new K_PATCH* [other_patch->num_trim_curves + 1];
    }
    
    //  Identify curves on *this and on the other patch.
    
    if (other_dir)
    {
      K_POINT2D**   this_seg_endpts;
      K_POINT2D**   other_seg_endpts;
      unsigned long num_other_segments;
      
      for (t2 = 0;
           t2 < other_patch->num_trim_curves
           &&
           other_patch->trim_curves[t2] != trim_curves[t]->curve_in_other_dom;
           t2++)
        ;
      
      assert(t2 < other_patch->num_trim_curves);
      
      this_seg_endpts    = new K_POINT2D* [1];
      this_seg_endpts[0] = trim_curves[t]->segments[s]->end;
      this_seg_endpts[0]->ref_count++;
      
      assert(num_other_segments = other_patch->trim_curves[t2]->num_segments);
      other_seg_endpts = new K_POINT2D* [num_other_segments];
      
      for (i = 0; i < num_other_segments; i++)
      {
        other_seg_endpts[i] = other_patch->trim_curves[t2]->segments[i]->end;
        other_seg_endpts[i]->ref_count++;
      }
      
      match_pts(this_seg_endpts, 1, surf,
                other_seg_endpts, num_other_segments, other_patch->surf);
      
      for (s2 = 0;
           s2 < other_patch->trim_curves[t2]->num_segments
           &&
           !other_patch->trim_curves[t2]->segments[s2]->end->
             equiv(*other_seg_endpts[0]);
           s2++)
        ;
      
      assert(s2 < other_patch->trim_curves[t2]->num_segments);
      
      if (!--this_seg_endpts[0]->ref_count)
        delete this_seg_endpts[0];
      
      delete [] this_seg_endpts;
      
      for (i = 0; i < s2; i++)
        if (!--other_seg_endpts[i]->ref_count)
          delete other_seg_endpts[i];
      
      delete [] other_seg_endpts;
    }
    
    //  Set the new curves on *this.
    
    for (i = 0; i < t; i++)
    {
//      new_this_trim_curves[i] = new K_CURVE(*trim_curves[i]);
      new_this_trim_curves[i] = trim_curves[i];
      new_this_trim_curves[i]->ref_count++;
//      new_this_trim_curves[i]->assoc(trim_curves[i]->curve_in_other_dom,
//                                     trim_curves[i]->dir_in_other_dom);
      new_this_adj_surfs[i]   = adj_surfs[i];
      new_this_adj_surfs[i]->ref_count++;
      new_this_adj_patches[i] = adj_patches[i];
      new_this_adj_patches[i]->ref_count++;
    }
    
    trim_curves[t]->split(s,
                          *(new_this_trim_curves[t] = new K_CURVE),
                          *(new_this_trim_curves[t + 1] = new K_CURVE));
    new_this_trim_curves[t]->ref_count++;
    new_this_trim_curves[t + 1]->ref_count++;
    new_this_adj_surfs[t]   = new_this_adj_surfs[t + 1]   = adj_surfs[t];
    new_this_adj_surfs[t]->ref_count++;
    new_this_adj_surfs[t + 1]->ref_count++;
    new_this_adj_patches[t] = new_this_adj_patches[t + 1] = adj_patches[t];
    new_this_adj_patches[t]->ref_count++;
    new_this_adj_patches[t + 1]->ref_count++;
    
    for (i = t + 1; i < num_trim_curves; i++)
    {
//      new_this_trim_curves[i + 1] = new K_CURVE(*trim_curves[i]);
      new_this_trim_curves[i + 1] = trim_curves[i];
      new_this_trim_curves[i + 1]->ref_count++;
//      new_this_trim_curves[i + 1]->assoc(trim_curves[i]->curve_in_other_dom,
//                                         trim_curves[i]->dir_in_other_dom);
      new_this_adj_surfs[i + 1]   = adj_surfs[i];
      new_this_adj_surfs[i + 1]->ref_count++;
      new_this_adj_patches[i + 1] = adj_patches[i];
      new_this_adj_patches[i + 1]->ref_count++;
    }
    
    //  Set the new curves on the other patch.
    
    if (other_dir)
    {      
      for (i = 0; i < t2; i++)
      {
//        new_other_trim_curves[i] = new K_CURVE(*other_patch->trim_curves[i]);
        new_other_trim_curves[i] = other_patch->trim_curves[i];
        new_other_trim_curves[i]->ref_count++;
//        new_other_trim_curves[i]->
//          assoc(other_patch->trim_curves[i]->curve_in_other_dom,
//                other_patch->trim_curves[i]->dir_in_other_dom);
        new_other_adj_surfs[i]   = other_patch->adj_surfs[i];
        new_other_adj_surfs[i]->ref_count++;
        new_other_adj_patches[i] = other_patch->adj_patches[i];
        new_other_adj_patches[i]->ref_count++;
      }
      
      other_patch->trim_curves[t2]->
        split(s2,
              *(new_other_trim_curves[t2] = new K_CURVE),
              *(new_other_trim_curves[t2 + 1] = new K_CURVE));
      new_other_trim_curves[t2]->ref_count++;
      new_other_trim_curves[t2 + 1]->ref_count++;
      new_other_adj_surfs[t2] = new_other_adj_surfs[t2 + 1]
                              = other_patch->adj_surfs[t2];
      new_other_adj_surfs[t2]->ref_count++;
      new_other_adj_surfs[t2 + 1]->ref_count++;
      new_other_adj_patches[t2] = new_other_adj_patches[t2 + 1]
                                = other_patch->adj_patches[t2];
      new_other_adj_patches[t2]->ref_count++;
      new_other_adj_patches[t2 + 1]->ref_count++;
      
      for (i = t2 + 1; i < other_patch->num_trim_curves; i++)
      {
//        new_other_trim_curves[i + 1] = new K_CURVE(*other_patch->trim_curves[i]);
        new_other_trim_curves[i + 1] = other_patch->trim_curves[i];
        new_other_trim_curves[i + 1]->ref_count++;
//        new_other_trim_curves[i + 1]->
//          assoc(other_patch->trim_curves[i]->curve_in_other_dom,
//                other_patch->trim_curves[i]->dir_in_other_dom);
        new_other_adj_surfs[i + 1]   = other_patch->adj_surfs[i];
        new_other_adj_surfs[i + 1]->ref_count++;
        new_other_adj_patches[i + 1] = other_patch->adj_patches[i];
        new_other_adj_patches[i + 1]->ref_count++;
      }
    }
    
    //  Set direction correspondences.
    
    if (other_dir > 0)
    {
      new_this_trim_curves[t]->assoc(new_other_trim_curves[t2], 1);
      new_this_trim_curves[t + 1]->assoc(new_other_trim_curves[t2 + 1], 1);
    }
    else if (other_dir < 0)
    {
      new_this_trim_curves[t]->assoc(new_other_trim_curves[t2 + 1], - 1);
      new_this_trim_curves[t + 1]->assoc(new_other_trim_curves[t2], - 1);
    }
    
    //  Update the other patch.
    //  CAUTION:  We update the other patch first,
    //            since other_patch = adj_patch[t] might be deleted
    //            when we update *this.
    
    if (other_dir)
    {
      for (i = 0; i < other_patch->num_trim_curves; i++)
      {
        if (!--other_patch->trim_curves[i]->ref_count)
          delete other_patch->trim_curves[i];
        
        if (other_patch->adj_surfs[i]
            &&
            !--other_patch->adj_surfs[i]->ref_count)
          delete other_patch->adj_surfs[i];
        
        if (other_patch->adj_patches[i]
            &&
            !--other_patch->adj_patches[i]->ref_count)
          delete other_patch->adj_patches[i];
      }
      
      //  num_trim_curves > 0 && other_patch = adj_patches[t]
      //    => other_patch->num_trim_curves > 0
      
      delete [] other_patch->trim_curves;
      delete [] other_patch->adj_surfs;
      delete [] other_patch->adj_patches;
      
      other_patch->num_trim_curves++;
      other_patch->trim_curves = new_other_trim_curves;
      other_patch->adj_surfs   = new_other_adj_surfs;
      other_patch->adj_patches = new_other_adj_patches;
    }
    
    //  Update *this.
    
    for (i = 0; i < num_trim_curves; i++)
    {
      if (!--trim_curves[i]->ref_count)
        delete trim_curves[i];
      
      if (adj_surfs[i] && !--adj_surfs[i]->ref_count)
        delete adj_surfs[i];
      
      if (adj_patches[i] && !--adj_patches[i]->ref_count)
        delete adj_patches[i];
    }
    
    delete [] trim_curves;  //  num_trim_curves > 0
    delete [] adj_surfs;    //  num_trim_curves > 0
    delete [] adj_patches;  //  num_trim_curves > 0
    
    num_trim_curves++;
    trim_curves = new_this_trim_curves;
    adj_surfs   = new_this_adj_surfs;
    adj_patches = new_this_adj_patches;
  }
  
  return 0;
}

unsigned long K_PATCH :: merge_curves()
{
  long          i, j, k;
  unsigned long num_merged_proto, num_ic_closed;
  
  //  Delete old informations.
  
  if (num_merged > 0)
  {
    for (i = 0; i < num_merged; i++)
    {
      for (j = 0; j < len_ic[i]; j++)
        if (!--ic[i][j]->ref_count)
          delete ic[i][j];
      
      if (len_ic[i] > 0)
        delete [] ic[i];
    }
    
    delete [] ic;
    delete [] len_ic;
    delete [] ic_closed;
  }
  
  num_merged = 0;
  
  //  Group trim curves and intersection curves to form partitions.
  //  Group intersection curves into contiguous sections.
  //  The endpts of intersection curves which match up with each other
  //  will NOT have the same location in memory!!!
  //  Must see if points are the same.
  
  if ((num_merged_proto = num_int_curves) > 0)
  {
    K_CURVE***     ic_proto;
    unsigned long* len_ic_proto;
    unsigned long* start_on_trim;
    unsigned long* end_on_trim;
    K_POINT2D**    int_endpts;
    int            cc, cs;
    long           num_potential, last;
    K_BOXCO2       b1, b2;
    
    ic_proto     = new K_CURVE** [num_merged_proto];
    len_ic_proto = new unsigned long [num_merged_proto];
    
    start_on_trim = new unsigned long [num_merged_proto];
    end_on_trim   = new unsigned long [num_merged_proto];
    
    for (i = 0; i < num_merged_proto; i++)
    {
      ic_proto[i]     = new K_CURVE* [num_merged_proto];
      ic_proto[i][0]  = int_curves[i];
      ic_proto[i][0]->ref_count++;
      len_ic_proto[i] = 1;
      
      start_on_trim[i] = 0;
      end_on_trim[i]   = 0;
    }
    
    //  Build a list of endpts of intersection curves, and
    //  break up trimming curves.
    
    int_endpts = new K_POINT2D* [2 * num_int_curves];
    
    for (i = 0; i < num_int_curves; i++)
    {
      int_endpts[2 * i]     = int_curves[i]->start();
      int_endpts[2 * i]->ref_count++;
      int_endpts[2 * i + 1] = int_curves[i]->end();
      int_endpts[2 * i + 1]->ref_count++;
    }
    
    //  Set start_on_trim and end_on_trim.
    
    for (i = 0; i < num_trim_curves; i++)
    {
//      cerr << endl << " kpatch: merge_curves: ---------- " << endl << flush;
//      cerr << " kpatch: merge_curves: ID = " << ID << endl << flush;
//      cerr << " kpatch: merge_curves: i = " << i << ", num_trim_curves = " << num_trim_curves << endl << flush;
//      for (j = 0; j < num_trim_curves; j++)
//        cerr << " kpatch: merge_curves: *trim_curves[" << j << "] = " << endl << *trim_curves[j] << endl << flush;
//      cerr << " kpatch: merge_curves: ---------- " << endl << flush;
//      cerr << " kpatch: merge_curves: int_curve = " << endl << *int_curves[0]->poly << endl << flush;
//      for (k = 0; k < num_int_curves; k++)
//        cerr << " kpatch: merge_curves: *int_curves[" << k << "] = " << endl << *int_curves[k] << endl << flush;
//      for (k = 0; k < 2 * num_int_curves; k++)
//        cerr << " kpatch: merge_curves: *int_endpts[" << k << "] = " << *int_endpts[k] << endl << flush;
//      cerr << " kpatch: merge_curves: 1: ========== " << endl << endl << flush;
      
      for (j = 0; j < trim_curves[i]->num_segments - 1; j++)
      {
        for (k = 0;
             k < 2 * num_int_curves
             &&
             !trim_curves[i]->segments[j]->end->equal(*int_endpts[k]);
             k++)
          ;
        
        if (k < 2 * num_int_curves)
        {
          split_trim_curves(i, j);
          
          if (!(k % 2))
            start_on_trim[k / 2] = 1;
          else  //  if (k % 2)
            end_on_trim[k / 2] = 1;
        }
      }
      
      for (k = 0;
           k < 2 * num_int_curves
           &&
           !trim_curves[i]->end()->equal(*int_endpts[k]);
           k++)
        ;
      
      if (k < 2 * num_int_curves)
      {
        if (!(k % 2))
          start_on_trim[k / 2] = 1;
        else  //  if (k % 2)
          end_on_trim[k / 2] = 1;
      }
    }
    
    for (i = 0; i < 2 * num_int_curves; i++)
      if (!--int_endpts[i]->ref_count)
        delete int_endpts[i];
    
    delete [] int_endpts;
    
    //  Match up int_endpts.
    
//    cerr << endl << " kpatch: merge_curves: ---------- " << endl << flush;
//    for (i = 0; i < num_merged_proto; i++)
//      cerr << "   end_on_trim[" << i << "] = " << end_on_trim[i] << ", start_on_trim[" << i << "] = " << start_on_trim[i] << endl << flush;
//    cerr << " kpatch: merge_curves: ========== " << endl << endl << flush;
    
    num_merged = num_merged_proto;
    i          = 0;
    
    while (i < num_merged)
    {
      //  Find all curves which potentially match up with this curve.
      
      if (!end_on_trim[i])
      //  Find all curves which potentially match up with this curve
      //  at this end pts.
      {
//        b1            = ic_proto[i][len_ic_proto[i] - 1]->end()->bbox();
        num_potential = last = 0;
        
        for (j = i + 1; j < num_merged; j++)
        {
          if (!start_on_trim[j])
          {
//            b2 = ic_proto[j][0]->start()->bbox();
//            
//            if (b1.overlap(b2))
            if (ic_proto[i][len_ic_proto[i] - 1]->end()->equal(*ic_proto[j][0]->start()))
            {
              num_potential++;
              last = j;
            }
          }
          
          if (!end_on_trim[j])
          {
//            b2 = ic_proto[j][len_ic_proto[j] - 1]->end()->bbox();
//            
//            if (b1.overlap(b2))
            if (ic_proto[i][len_ic_proto[i] - 1]->end()->equal(*ic_proto[j][len_ic_proto[j] - 1]->end()))
            {
              num_potential++;
              last = - j;
            }
          }
        }
        
        if (num_potential == 1)
          if (last > 0)
          //  Match this end pt to the start pt of the other curve.
          {
            for (j = 0; j < len_ic_proto[last]; j++)
            {
              ic_proto[i][len_ic_proto[i]] = ic_proto[last][j];
              ic_proto[i][len_ic_proto[i]]->ref_count++;
              len_ic_proto[i]++;
            }
            
            end_on_trim[i] = end_on_trim[last];
            
            for (j = 0; j < len_ic_proto[num_merged - 1]; j++)
            {
              if (j < len_ic_proto[last])
                ic_proto[last][j]->ref_count--;
              
              ic_proto[last][j] = ic_proto[num_merged - 1][j];
              ic_proto[last][j]->ref_count++;
            }
            
            len_ic_proto[last]  = len_ic_proto[num_merged - 1];
            start_on_trim[last] = start_on_trim[num_merged - 1];
            end_on_trim[last]   = end_on_trim[num_merged - 1];
            num_merged--;
          }
          else  //  if (last < 0)
          //  Match this end pt to the end pt of the other curve.
          {
            last *= - 1;
            
            for (j = len_ic_proto[last] - 1; j >= 0; j--)
            {
              ic_proto[last][j]->reverse();
              ic_proto[i][len_ic_proto[i]] = ic_proto[last][j];
              ic_proto[i][len_ic_proto[i]]->ref_count++;
              len_ic_proto[i]++;
            }
            
            end_on_trim[i] = start_on_trim[last];
            
            for (j = 0; j < len_ic_proto[num_merged - 1]; j++)
            {
              if (j < len_ic_proto[last])
                ic_proto[last][j]->ref_count--;
              
              ic_proto[last][j] = ic_proto[num_merged - 1][j];
              ic_proto[last][j]->ref_count++;
            }
            
            len_ic_proto[last]  = len_ic_proto[num_merged - 1];
            start_on_trim[last] = start_on_trim[num_merged - 1];
            end_on_trim[last]   = end_on_trim[num_merged - 1];
            num_merged--;
          }
        else  //  if (!num_potential || num_potential > 1)
          if (!num_potential)
          //  Double-check this curve is closed.
          {
//            b2 = ic_proto[i][0]->start()->bbox();
//            assert(b1.overlap(b2));
            assert(ic_proto[i][len_ic_proto[i] - 1]->end()->
                   equal(*ic_proto[i][0]->start()));
            i++;
          }
          else  //  if (num_potential > 1)
            ic_proto[i][len_ic_proto[i] - 1]->end()->
              shrink(shrink_step, shrink_step);
      }
      else if (!start_on_trim[i])
      //  Find all curves which potentially match up with this curve
      //  at this start pts.
      {
//        b1            = ic_proto[i][0]->start()->bbox();
        num_potential = last = 0;
        
        for (j = i + 1; j < num_merged; j++)
        {
          if (!start_on_trim[j])
          {
//            b2 = ic_proto[j][0]->start()->bbox();
//            
//            if (b1.overlap(b2))
            if (ic_proto[i][0]->start()->equal(*ic_proto[j][0]->start()))
            {
              num_potential++;
              last = j;
            }
          }
          
          if (!end_on_trim[j])
          {
//            b2 = ic_proto[j][len_ic_proto[j]- 1]->end()->bbox();
//            
//            if (b1.overlap(b2))
            if (ic_proto[i][0]->start()->equal(*ic_proto[j][len_ic_proto[j] - 1]->end()))
            {
              num_potential++;
              last = - j;
            }
          }
        }
        
        if (num_potential == 1)
          if (last > 0)
          //  Match this start pt to the start pt of the other curve.
          {
            for (j = len_ic_proto[i] - 1; j >= 0; j--)
            {
              if (j + len_ic_proto[last] < len_ic_proto[i])
                ic_proto[i][j + len_ic_proto[last]]->ref_count--;
              
              ic_proto[i][j + len_ic_proto[last]] = ic_proto[i][j];
              ic_proto[i][j + len_ic_proto[last]]->ref_count++;
            }
            
            for (j = 0; j < len_ic_proto[last]; j++)
            {
              ic_proto[last][j]->reverse();
              ic_proto[i][j]->ref_count--;
              ic_proto[i][j] = ic_proto[last][j];
              ic_proto[i][j]->ref_count++;
            }
            
            len_ic_proto[i]  += len_ic_proto[last];
            start_on_trim[i]  = end_on_trim[last];
            
            for (j = 0; j < len_ic_proto[num_merged - 1]; j++)
            {
              if (j < len_ic_proto[last])
                ic_proto[last][j]->ref_count--;
              
              ic_proto[last][j] = ic_proto[num_merged - 1][j];
              ic_proto[last][j]->ref_count++;
            }
            
            len_ic_proto[last]  = len_ic_proto[num_merged - 1];
            start_on_trim[last] = start_on_trim[num_merged - 1];
            end_on_trim[last]   = end_on_trim[num_merged - 1];
            num_merged--;
          }
          else  //  if (last < 0)
          //  Match this start pt to the end pt of the other curve.
          {
            last *= - 1;
            
            for (j = len_ic_proto[i] - 1; j >= 0; j--)
            {
              if (j + len_ic_proto[last] < len_ic_proto[i])
                ic_proto[i][j + len_ic_proto[last]]->ref_count--;
              
              ic_proto[i][j + len_ic_proto[last]] = ic_proto[i][j];
              ic_proto[i][j + len_ic_proto[last]]->ref_count++;
            }
            
            for (j = 0; j < len_ic_proto[last]; j++)
            {
              ic_proto[i][j]->ref_count--;
              ic_proto[i][j] = ic_proto[last][j];
              ic_proto[i][j]->ref_count++;
              len_ic_proto[i]++;
            }
            
            start_on_trim[i] = start_on_trim[last];
            
            for (j = 0; j < len_ic_proto[num_merged - 1]; j++)
            {
              if (j < len_ic_proto[last])
                ic_proto[last][j]->ref_count--;
              
              ic_proto[last][j] = ic_proto[num_merged - 1][j];
              ic_proto[last][j]->ref_count++;
            }
            
            len_ic_proto[last]  = len_ic_proto[num_merged - 1];
            start_on_trim[last] = start_on_trim[num_merged - 1];
            end_on_trim[last]   = end_on_trim[num_merged - 1];
            num_merged--;
          }
        else  //  if (!num_potential || num_potential > 1)
        {
          assert(num_potential > 1);
          ic_proto[i][0]->start()->shrink(shrink_step, shrink_step);
        }
      }
      else
        i++;
    }
    
    assert(num_merged > 0);
    
    //  Update ic.
    
    ic        = new K_CURVE** [num_merged];
    len_ic    = new unsigned long [num_merged];
    ic_closed = new int [num_merged];
    
    for (i = 0; i < num_merged; i++)
    {
      ic[i] = new K_CURVE* [len_ic[i] = len_ic_proto[i]];
      
      for (j = 0; j < len_ic[i]; j++)
      {
        ic[i][j] = ic_proto[i][j];
        ic[i][j]->ref_count++;
      }
    }
    
    //  Count closed curves.
    
    num_ic_closed = 0;
    
    for (i = 0; i < num_merged; i++)
      if (!start_on_trim[i] && !end_on_trim[i])
      {
        ic_closed[i] = 1;
        num_ic_closed++;
      }
      else  //  if (start_on_trim[i] || end_on_trm[i])
        ic_closed[i] = 0;
    
    for (i = 0; i < num_merged_proto; i++)
    {
      for (j = 0; j < len_ic_proto[i]; j++)
        if (!--ic_proto[i][j]->ref_count)
          delete ic_proto[i][j];
      
      delete [] ic_proto[i];  //  num_merged_proto > 0
    }
    
    delete [] ic_proto;      //  num_merged_proto > 0
    delete [] len_ic_proto;  //  num_merged_proto > 0
    
    delete [] start_on_trim;
    delete [] end_on_trim;
  }
  else  //  if (!num_merged)
  {
    ic        = 0;
    len_ic    = 0;
    ic_closed = 0;
    
    num_ic_closed = 0;
  }
  
  return num_ic_closed;
}

#define MAX_NUM_PTS_SPLIT 128

unsigned long K_PATCH :: split_loops(K_PATCH**& new_patches)
{
  assert(num_trim_curves > 0);
  assert(num_merged > 0);
  
  long           i, j, k;
  unsigned long  loop;
  K_BOXCO2       b;
  unsigned int   dir_split;
  bigrational    v_split;
  K_POINT2D**    pts_split;
  unsigned long  num_pts_split;
  K_SEGMENT*     segments_split[1];
  K_RATPOLY      poly_split;
  K_CURVE        curve_split;
  K_SURF*        new_surf;
  K_POINT2D**    int_pts1;
  unsigned long  num_int_pts1;
  K_POINT2D**    int_pts2;
  unsigned long  num_int_pts2;
  int            go_on;
  unsigned long  i_c1, i_s1, i_p1, i_c2, i_s2, i_p2;
  K_CURVE*       this_curve;
  K_CURVE*       new_this_curve;
  K_CURVE**      new_this_trim_curves;
  K_SURF**       new_this_adj_surfs;
  K_PATCH**      new_this_adj_patches;
  unsigned long  num_new_this_trim_curves;
  K_CURVE**      new_other_trim_curves;
  K_SURF**       new_other_adj_surfs;
  K_PATCH**      new_other_adj_patches;
  K_CURVE**      curves_split;
  unsigned long  num_curves_split;
  K_CURVE        cs1, cs2;
  K_PATCH**      sub_patches;
  unsigned long  num_sub_patches;
  int*           used;
  int            done;
  K_CURVE**      this_trim_curves;
  K_SURF**       this_adj_surfs;
  K_PATCH**      this_adj_patches;
  unsigned long  num_this_trim_curves;
  K_POINT2D*     this_start;
  K_POINT2D*     this_end;
  int            last_on_split;
  K_POINT2D*     int_start;
  K_POINT2D*     int_end;
  K_CURVE**      new_this_int_curves;
  K_SURF**       new_this_adj_int_surfs;
  K_PATCH**      new_this_adj_int_patches;
  unsigned long  num_new_this_int_curves;
  K_CURVE**      new_other_int_curves;
  K_SURF**       new_other_adj_int_surfs;
  K_PATCH**      new_other_adj_int_patches;
  unsigned long  num_new_other_int_curves;
  K_POINT2D      x;
  K_CURVE**      this_int_curves;
  K_SURF**       this_adj_int_surfs;
  K_PATCH**      this_adj_int_patches;
  unsigned long  num_this_int_curves;
  K_PATCH***     new_patches_proto;
  unsigned long* num_new_patches_proto;
  unsigned long  num_new_patches;
  
  loop = num_merged;
  
  for (i = 0; i < num_merged && loop == num_merged; i++)
    if (ic_closed[i])
      loop = i;
  
  if (loop < num_merged)
  {
    pts_split     = new K_POINT2D* [MAX_NUM_PTS_SPLIT];
    num_pts_split = 0;
    
    b = ic[loop][0]->bbox();
    
    for (i = 1; i < len_ic[loop]; i++)
      b = b.merge(ic[loop][i]->bbox());
    
    if (b.high[0] - b.low[0] > b.high[1] - b.low[1])
    {
      dir_split         = 0;
      v_split           = (b.low[0] + b.high[0]) / 2;
      segments_split[0] = new K_SEGMENT(new K_POINT2D(v_split, low_t),
                                        new K_POINT2D(v_split, high_t));
      segments_split[0]->ref_count++;
      poly_split        = K_RATPOLY(2, 0, v_split);
      curve_split       = K_CURVE(poly_split, segments_split, 1);
    }
    else  //  if (b.high[0] - b.low[0] <= b.high[1] - b.low[1])
    {
      dir_split         = 1;
      v_split           = (b.low[1] + b.high[1]) / 2;
      segments_split[0] = new K_SEGMENT(new K_POINT2D(low_s, v_split),
                                        new K_POINT2D(high_s, v_split));
      segments_split[0]->ref_count++;
      poly_split        = K_RATPOLY(2, 1, v_split);
      curve_split       = K_CURVE(poly_split, segments_split, 1);
    }
    
    new_surf = new K_SURF(surf->split_surf(v_split, dir_split));
    
    //  Split trim curves.
    
    new_this_trim_curves     = new K_CURVE* [MAX_NUM_PTS_SPLIT];
    new_this_adj_surfs       = new K_SURF* [MAX_NUM_PTS_SPLIT];
    new_this_adj_patches     = new K_PATCH* [MAX_NUM_PTS_SPLIT];
    num_new_this_trim_curves = 0;
    
    for (i = 0; i < num_trim_curves; i++)
    {
      if ((num_int_pts1 =
           trim_curves[i]->find_intersections(poly_split, int_pts1, 1)) > 0)
      {
//        num_int_pts2 = get_all_pts(*trim_curves[i]->curve_in_other_dom->poly,
//                                   new_surf->Impl->subst_param_expr(
//                                     *adj_patches[i]->surf->X,
//                                     *adj_patches[i]->surf->Y,
//                                     *adj_patches[i]->surf->Z,
//                                     *adj_patches[i]->surf->W),
//                                   int_pts2,
//                                   0);
        num_int_pts2 =
          get_all_int_pts(*trim_curves[i]->curve_in_other_dom->poly,
                          *new_surf, *adj_patches[i]->surf,
                          int_pts2);
        assert(num_int_pts1 <= num_int_pts2);
        match_pts(int_pts1, num_int_pts1, surf,
                  int_pts2, num_int_pts2, adj_patches[i]->surf);
        
        for (j = 0; j < num_int_pts1; j++)
        {
          //  Add int_pts1[j] & int_pts2[j] to trim_curves.
          
          if (!int_pts1[j]->equal(*trim_curves[i]->end()))
          {
            pt_surf_assoc.record_as_assoc(int_pts1[j], surf,
                                          int_pts2[j], adj_patches[i]->surf);
            
            assert(trim_curves[i]->add_pt(int_pts1[j], 0));
            assert(trim_curves[i]->curve_in_other_dom->add_pt(int_pts2[j], 1));
          }
          
          //  Add int_pts1[j] to pts_split.
          
          assert(num_pts_split < MAX_NUM_PTS_SPLIT);
          pts_split[num_pts_split] = int_pts1[j];
          pts_split[num_pts_split]->ref_count++;
          num_pts_split++;
        }
        
        this_curve = trim_curves[i];
        this_curve->ref_count++;
        
        for (j = 0; j < num_int_pts1; j++)
        {
          //  Locate int_pts[i_p1] on this_curve->segments[i_s1].
          
          i_p1 = num_int_pts1;
          i_s1 = 0;
          
          while (i_p1 == num_int_pts1 && i_s1 < this_curve->num_segments - 1)
          {
            for (i_p1 = 0;
                 i_p1 < num_int_pts1
                 &&
                 !int_pts1[i_p1]->equal(*this_curve->segments[i_s1]->end);
                 i_p1++)
              ;
            
            if (i_p1 == num_int_pts1)
              i_s1++;
          }
          
          if (i_s1 < this_curve->num_segments - 1)
          //  && if (i_p1 < num_int_pts1)
          {
            //  Split this domain.
            
            this_curve->
              split(i_s1,
                    *(new_this_trim_curves[num_new_this_trim_curves] =
                      new K_CURVE),
                    *(new_this_curve = new K_CURVE));
            new_this_trim_curves[num_new_this_trim_curves]->ref_count++;
            new_this_curve->ref_count++;
            new_this_adj_surfs[num_new_this_trim_curves]   = adj_surfs[i];
            new_this_adj_surfs[num_new_this_trim_curves]->ref_count++;
            new_this_adj_patches[num_new_this_trim_curves] = adj_patches[i];
            new_this_adj_patches[num_new_this_trim_curves]->ref_count++;
//            num_new_this_trim_curves++;  //  Do this later.
            
            //  Split the other domain.
            
            for (i_s2 = 0;
                 i_s2 < this_curve->curve_in_other_dom->num_segments
                 &&
                 !this_curve->curve_in_other_dom->segments[i_s2]->end->
                   equal(*int_pts2[i_p1]);
                 i_s2++)
              ;
            
            for (i_c2 = 0;
                 i_c2 < adj_patches[i]->num_trim_curves
                 &&
                 adj_patches[i]->trim_curves[i_c2] !=
                   this_curve->curve_in_other_dom;
                 i_c2++)
              ;
            
            assert(i_s2 < this_curve->curve_in_other_dom->num_segments);
            assert(i_c2 < adj_patches[i]->num_trim_curves);
            
            new_other_trim_curves =
              new K_CURVE* [adj_patches[i]->num_trim_curves + 1];
            new_other_adj_surfs =
              new K_SURF* [adj_patches[i]->num_trim_curves + 1];
            new_other_adj_patches =
              new K_PATCH* [adj_patches[i]->num_trim_curves + 1];
            
            for (k = 0; k < i_c2; k++)
            {
              new_other_trim_curves[k] = adj_patches[i]->trim_curves[k];
              new_other_trim_curves[k]->ref_count++;
              new_other_adj_surfs[k]   = adj_patches[i]->adj_surfs[k];
              new_other_adj_surfs[k]->ref_count++;
              new_other_adj_patches[k] = adj_patches[i]->adj_patches[k];
              new_other_adj_patches[k]->ref_count++;
            }
            
            adj_patches[i]->trim_curves[i_c2]->
              split(i_s2,
                    *(new_other_trim_curves[i_c2] = new K_CURVE),
                    *(new_other_trim_curves[i_c2 + 1] = new K_CURVE));
            new_other_trim_curves[i_c2]->ref_count++;
            new_other_trim_curves[i_c2 + 1]->ref_count++;
            new_other_adj_surfs[i_c2] = new_other_adj_surfs[i_c2 + 1] =
              adj_patches[i]->adj_surfs[i_c2];
            new_other_adj_surfs[i_c2]->ref_count++;
            new_other_adj_surfs[i_c2 + 1]->ref_count++;
            new_other_adj_patches[i_c2] = new_other_adj_patches[i_c2 + 1] =
              adj_patches[i]->adj_patches[i_c2];
            new_other_adj_patches[i_c2]->ref_count++;
            new_other_adj_patches[i_c2 + 1]->ref_count++;
            
            for (k = i_c2 + 1; k < adj_patches[i]->num_trim_curves; k++)
            {
              new_other_trim_curves[k + 1] = adj_patches[i]->trim_curves[k];
              new_other_trim_curves[k + 1]->ref_count++;
              new_other_adj_surfs[k + 1]   = adj_patches[i]->adj_surfs[k];
              new_other_adj_surfs[k + 1]->ref_count++;
              new_other_adj_patches[k + 1] = adj_patches[i]->adj_patches[k];
              new_other_adj_patches[k + 1]->ref_count++;
            }
            
            for (k = 0; k < adj_patches[i]->num_trim_curves; k++)
            {
              if (!--adj_patches[i]->trim_curves[k]->ref_count)
                delete adj_patches[i]->trim_curves[k];
              
              if (!--adj_patches[i]->adj_surfs[k]->ref_count)
                delete adj_patches[i]->adj_surfs[k];
              
              if (!--adj_patches[i]->adj_patches[k]->ref_count)
                delete adj_patches[i]->adj_patches[k];
            }
            
            delete [] adj_patches[i]->trim_curves;
            delete [] adj_patches[i]->adj_surfs;
            delete [] adj_patches[i]->adj_patches;
            
            adj_patches[i]->trim_curves = new_other_trim_curves;
            adj_patches[i]->adj_surfs   = new_other_adj_surfs;
            adj_patches[i]->adj_patches = new_other_adj_patches;
            adj_patches[i]->num_trim_curves++;
            
            //  Associate curves in both domains.
            
            assert(this_curve->dir_in_other_dom);
            
            if (this_curve->dir_in_other_dom > 0)
            {
              new_this_trim_curves[num_new_this_trim_curves]->
                assoc(adj_patches[i]->trim_curves[i_c2], 1);
              new_this_curve->assoc(adj_patches[i]->trim_curves[i_c2 + 1], 1);
            }
            else  //  if (this_curve->dir_in_other_dom < 0)
            {
              new_this_trim_curves[num_new_this_trim_curves]->
                assoc(adj_patches[i]->trim_curves[i_c2 + 1], - 1);
              new_this_curve->assoc(adj_patches[i]->trim_curves[i_c2], - 1);
            }
            
            num_new_this_trim_curves++;
            
            if (!--this_curve->ref_count)
              delete this_curve;
            
            this_curve = new_this_curve;
            this_curve->ref_count++;
          }
        }
        
        new_this_trim_curves[num_new_this_trim_curves] = this_curve;
        new_this_trim_curves[num_new_this_trim_curves]->ref_count++;
        new_this_adj_surfs[num_new_this_trim_curves]   = adj_surfs[i];
        new_this_adj_surfs[num_new_this_trim_curves]->ref_count++;
        new_this_adj_patches[num_new_this_trim_curves] = adj_patches[i];
        new_this_adj_patches[num_new_this_trim_curves]->ref_count++;
        num_new_this_trim_curves++;
        
        for (j = 0; j < num_int_pts1; j++)
          if (!--int_pts1[j]->ref_count)
            delete int_pts1[j];
        
        delete [] int_pts1;  //  num_int_pts1 > 0 => int_pts1 != 0
        
        if (num_int_pts2 > 0)
        {
          for (j = 0; j < num_int_pts2; j++)
            if (!--int_pts2[j]->ref_count)
              delete int_pts2[j];
          
          delete [] int_pts2;
        }
        
        if (!--this_curve->ref_count)
          delete this_curve;
      }
      else  //  if (!num_int_pts1)
      {
        new_this_trim_curves[num_new_this_trim_curves] = trim_curves[i];
        new_this_trim_curves[num_new_this_trim_curves]->ref_count++;
        new_this_adj_surfs[num_new_this_trim_curves]   = adj_surfs[i];
        new_this_adj_surfs[num_new_this_trim_curves]->ref_count++;
        new_this_adj_patches[num_new_this_trim_curves] = adj_patches[i];
        new_this_adj_patches[num_new_this_trim_curves]->ref_count++;
        num_new_this_trim_curves++;
      }
    }
    
    //  Generate curves on curve_split.
    
    assert(!(num_pts_split % 2));
    
    curve_split.sort_pts(pts_split, num_pts_split);
    curves_split     = new K_CURVE* [num_pts_split];
    num_curves_split = num_pts_split / 2;
    
    for (i = 0; i < num_curves_split; i++)
    {
      curve_split.add_pt(pts_split[2 * i]);
      curve_split.split(0, cs1, cs2);
      cs2.add_pt(pts_split[2 * i + 1]);
      cs2.split(0, *(curves_split[i] = new K_CURVE), curve_split);
      curves_split[i]->ref_count++;
      curves_split[i + num_curves_split] = new K_CURVE(*curves_split[i]);
      curves_split[i + num_curves_split]->ref_count++;
      curves_split[i + num_curves_split]->reverse();
      curves_split[i]->assoc(curves_split[i + num_curves_split], - 1);
    }
    
    if (!--segments_split[0]->ref_count)
      delete segments_split[0];
    
    if (num_pts_split > 0)
    {
      for (i = 0; i < num_pts_split; i++)
        if (!--pts_split[i]->ref_count)
          delete pts_split[i];
      
      delete [] pts_split;
    }
    
    //  Generate new patches by tracing trim curves.
    
    sub_patches     = new K_PATCH* [MAX_NUM_PTS_SPLIT];
    num_sub_patches = 0;
    
    used  = new int [num_new_this_trim_curves];
    
    for (i = 0; i < num_new_this_trim_curves; i++)
      used[i] = 0;
    
    done = 0;
    
    while (!done)
    {
      sub_patches[num_sub_patches] = new K_PATCH;
      sub_patches[num_sub_patches]->ref_count++;
      this_trim_curves             = new K_CURVE* [MAX_NUM_PTS_SPLIT];
      this_adj_surfs               = new K_SURF* [MAX_NUM_PTS_SPLIT];
      this_adj_patches             = new K_PATCH* [MAX_NUM_PTS_SPLIT];
      
      //  Find a trim_curve to start at.
      
      for (i = 0; i < num_new_this_trim_curves && used[i]; i++)
        ;
      
      this_trim_curves[0]  = new_this_trim_curves[i];
      this_trim_curves[0]->ref_count++;
      this_adj_surfs[0]    = new_this_adj_surfs[i];
      this_adj_surfs[0]->ref_count++;
      this_adj_patches[0]  = new_this_adj_patches[i];
      this_adj_patches[0]->ref_count++;
      this_start           = new_this_trim_curves[i]->start();
      this_end             = new_this_trim_curves[i]->end();
      num_this_trim_curves = 1;
      
      for (i_c2 = 0;
           i_c2 < new_this_adj_patches[i]->num_trim_curves
           &&
           new_this_adj_patches[i]->trim_curves[i_c2] !=
             new_this_trim_curves[i]->curve_in_other_dom;
           i_c2++)
        ;
      
      assert(i_c2 < new_this_adj_patches[i]->num_trim_curves);
      
      if (!--new_this_adj_patches[i]->adj_patches[i_c2]->ref_count)
        delete new_this_adj_patches[i]->adj_patches[i_c2];
      
      new_this_adj_patches[i]->adj_patches[i_c2] =
        sub_patches[num_sub_patches];
      new_this_adj_patches[i]->adj_patches[i_c2]->ref_count++;
      
      used[i]       = 1;
      go_on         = 1;
      last_on_split = 0;
      
      while (go_on)
      {
        //  Locate curves_split[j] starting at *this_end.
        
        for (j = 0;
             j < num_pts_split
             &&
             !curves_split[j]->start()->overlap(*this_end);
             j++)
          ;
        
        if (j < num_pts_split && !last_on_split)
        {
          this_trim_curves[num_this_trim_curves] = curves_split[j];
          this_trim_curves[num_this_trim_curves]->ref_count++;
          this_adj_surfs[num_this_trim_curves]   = new_surf;
          this_adj_surfs[num_this_trim_curves]->ref_count++;
          this_adj_patches[num_this_trim_curves] = 0;
          this_end                               = curves_split[j]->end();
          num_this_trim_curves++;
          
          last_on_split = 1;
        }
        else  //  if (j == num_pts_split || last_on_split)
        {
          for (i = 0;
               i < num_new_this_trim_curves
               &&
               (used[i]
                ||
                !new_this_trim_curves[i]->start()->overlap(*this_end));
               i++)
            ;
          
          this_trim_curves[num_this_trim_curves] = new_this_trim_curves[i];
          this_trim_curves[num_this_trim_curves]->ref_count++;
          this_adj_surfs[num_this_trim_curves]   = new_this_adj_surfs[i];
          this_adj_surfs[num_this_trim_curves]->ref_count++;
          this_adj_patches[num_this_trim_curves] = new_this_adj_patches[i];
          this_adj_patches[num_this_trim_curves]->ref_count++;
          this_end                               =
            new_this_trim_curves[i]->end();
          num_this_trim_curves++;
          
          for (i_c2 = 0;
               i_c2 < new_this_adj_patches[i]->num_trim_curves
               &&
               new_this_adj_patches[i]->trim_curves[i_c2] !=
                 new_this_trim_curves[i]->curve_in_other_dom;
               i_c2++)
            ;
          
          assert(i_c2 < new_this_adj_patches[i]->num_trim_curves);
          
          if (!--new_this_adj_patches[i]->adj_patches[i_c2]->ref_count)
            delete new_this_adj_patches[i]->adj_patches[i_c2];
          
          new_this_adj_patches[i]->adj_patches[i_c2] =
            sub_patches[num_sub_patches];
          new_this_adj_patches[i]->adj_patches[i_c2]->ref_count++;
          
          used[i]       = 1;
          last_on_split = 0;
        }
        
        if (this_start->equal(*this_end))
//        if (this_start->overlap(*this_end))
          go_on = 0;
      }
      
      //  Set members of sub_patches[num_sub_patches]
      //  which are initialized to 0.
      
      sub_patches[num_sub_patches]->surf = surf;
      sub_patches[num_sub_patches]->surf->ref_count++;
      sub_patches[num_sub_patches]->trim_curves =
        new K_CURVE* [num_this_trim_curves];
      sub_patches[num_sub_patches]->adj_surfs   =
        new K_SURF* [num_this_trim_curves];
      sub_patches[num_sub_patches]->adj_patches =
        new K_PATCH* [num_this_trim_curves];
      
      for (i = 0; i < num_this_trim_curves; i++)
      {
        sub_patches[num_sub_patches]->trim_curves[i] = this_trim_curves[i];
        sub_patches[num_sub_patches]->trim_curves[i]->ref_count++;
        sub_patches[num_sub_patches]->adj_surfs[i]   = this_adj_surfs[i];
        sub_patches[num_sub_patches]->adj_surfs[i]->ref_count++;
        
        if (sub_patches[num_sub_patches]->adj_patches[i] = this_adj_patches[i])
          sub_patches[num_sub_patches]->adj_patches[i]->ref_count++;
      }
      
      sub_patches[num_sub_patches]->num_trim_curves = num_this_trim_curves;
      
      for (i = 0; i < num_this_trim_curves; i++)
      {
        if (!--this_trim_curves[i]->ref_count)
          delete this_trim_curves[i];
        
        if (!--this_adj_surfs[i]->ref_count)
          delete this_adj_surfs[i];
        
        if (this_adj_patches[i] && !--this_adj_patches[i]->ref_count)
          delete this_adj_patches[i];
      }
      
      delete [] this_trim_curves;  //  MAX_NUM_PTS_SPLIT > 0
      delete [] this_adj_surfs;    //  MAX_NUM_PTS_SPLIT > 0
      delete [] this_adj_patches;  //  MAX_NUM_PTS_SPLIT > 0
      
      sub_patches[num_sub_patches]->set_range();
      num_sub_patches++;
      
      //  See whether or not done.
      
      for (i = 0; i < num_new_this_trim_curves && used[i]; i++)
        ;
      
      if (i == num_new_this_trim_curves)
        done = 1;
    }
    
    delete [] used;
    
    //  Update adj_patches of curves_split.
    
    for (i = 0; i < num_curves_split; i++)
    {
      //  Locate curves_split[i] on sub_patches[i_p1]->trim_curves[i_c1].
      
      i_p1 = 0;
      i_c1 = sub_patches[i_p1]->num_trim_curves;
      
      while (i_p1 < num_sub_patches
             &&
             i_c1 == sub_patches[i_p1]->num_trim_curves)
      {
        for (i_c1 = 0;
             i_c1 < sub_patches[i_p1]->num_trim_curves
             &&
             sub_patches[i_p1]->trim_curves[i_c1] != curves_split[i];
             i_c1++)
          ;
        
        if (i_c1 == sub_patches[i_p1]->num_trim_curves)
          i_c1 = sub_patches[++i_p1]->num_trim_curves;
      }
      
      //  Locate curves_split[i + num_curves_split]
      //  on sub_patches[i_p2]->trim_curves[i_c2].
      
      i_p2 = 0;
      i_c2 = sub_patches[i_p2]->num_trim_curves;
      
      while (i_p2 < num_sub_patches
             &&
             i_c2 == sub_patches[i_p2]->num_trim_curves)
      {
        for (i_c2 = 0;
             i_c2 < sub_patches[i_p2]->num_trim_curves
             &&
             sub_patches[i_p2]->trim_curves[i_c2] !=
               curves_split[i + num_curves_split];
             i_c2++)
          ;
        
        if (i_c2 == sub_patches[i_p2]->num_trim_curves)
          i_c2 = sub_patches[++i_p2]->num_trim_curves;
      }
      
      sub_patches[i_p1]->adj_patches[i_c1] = sub_patches[i_p2];
      sub_patches[i_p1]->adj_patches[i_c1]->ref_count++;
      sub_patches[i_p2]->adj_patches[i_c2] = sub_patches[i_p1];
      sub_patches[i_p2]->adj_patches[i_c2]->ref_count++;
    }
    
    for (i = 0; i < num_new_this_trim_curves; i++)
    {
      if (!--new_this_trim_curves[i]->ref_count)
        delete new_this_trim_curves[i];
      
      if (!--new_this_adj_surfs[i]->ref_count)
        delete new_this_adj_surfs[i];
      
      if (!--new_this_adj_patches[i]->ref_count)
        delete new_this_adj_patches[i];
    }
    
    delete [] new_this_trim_curves;  //  MAX_NUM_PTS_SPLIT > 0
    delete [] new_this_adj_surfs;    //  MAX_NUM_PTS_SPLIT > 0
    delete [] new_this_adj_patches;  //  MAX_NUM_PTS_SPLIT > 0
    
    //  Split intersection curves.
    
    new_this_int_curves      = new K_CURVE* [MAX_NUM_PTS_SPLIT];
    new_this_adj_int_surfs   = new K_SURF* [MAX_NUM_PTS_SPLIT];
    new_this_adj_int_patches = new K_PATCH* [MAX_NUM_PTS_SPLIT];
    num_new_this_int_curves  = 0;
    
    for (i = 0; i < num_int_curves; i++)
    {
      if (!dir_split)
        int_curves[i]->start()->cut_s(v_split);
      else  //  if (dir_split == 1)
        int_curves[i]->start()->cut_t(v_split);
      
      if (!dir_split
          &&
          (int_curves[i]->start()->type == 3
           ||
           int_curves[i]->start()->type == 4)
          ||
          dir_split == 1
          &&
          (int_curves[i]->start()->type == 2
           ||
           int_curves[i]->start()->type == 4))
      {
        for (k = 0; k < num_curves_split; k++)
          if (curves_split[k]->contains(*int_curves[i]->start()))
          {
            curves_split[k]->add_pt(int_curves[i]->start());
            curves_split[k + num_curves_split]->add_pt(int_curves[i]->start());
          }
      }
      
      num_int_pts1 =
        int_curves[i]->find_intersections(poly_split, int_pts1, 1);
      
      j = 0;
      
      while (j < num_int_pts1)
        if (int_pts1[j]->equal(*int_curves[i]->end()))
        //  int_pts1[j] has already been located at the end of int_curves[i].
        //  int_curves[i] does not have to be split.
        //  Add int_pts1[j] to some curves_split.
        {
          for (k = 0; k < num_curves_split; k++)
            if (curves_split[k]->contains(*int_curves[i]->end()))
            {
              curves_split[k]->add_pt(int_curves[i]->end());
              curves_split[k + num_curves_split]->add_pt(int_curves[i]->end());
            }
          
          if (!--int_pts1[j]->ref_count)
            delete int_pts1[j];
          
          for (k = j + 1; k < num_int_pts1; k++)
            int_pts1[k - 1] = int_pts1[k];
          
          num_int_pts1--;
        }
        else  //  if (!int_pts1[j]->equal(*int_curves[i]->end()))
          j++;
      
      if (num_int_pts1 > 0)
      {
//        num_int_pts2 = get_all_pts(*int_curves[i]->curve_in_other_dom->poly,
//                                   new_surf->Impl->subst_param_expr(
//                                     *adj_int_patches[i]->surf->X,
//                                     *adj_int_patches[i]->surf->Y,
//                                     *adj_int_patches[i]->surf->Z,
//                                     *adj_int_patches[i]->surf->W),
//                                   int_pts2,
//                                   0);
        num_int_pts2 =
          get_all_int_pts(*int_curves[i]->curve_in_other_dom->poly,
                          *new_surf, *adj_int_patches[i]->surf,
                          int_pts2);
        assert(num_int_pts1 <= num_int_pts2);
        match_pts(int_pts1, num_int_pts1, surf,
                  int_pts2, num_int_pts2, adj_int_patches[i]->surf);
        
        for (j = 0; j < num_int_pts1; j++)
        {
          //  Add int_pts1[j] & int_pts2[j] to int_curves.
          
          pt_surf_assoc.record_as_assoc(int_pts1[j], surf,
                                        int_pts2[j], adj_int_patches[i]->surf);
          
          assert(int_curves[i]->add_pt(int_pts1[j], 0));
          assert(int_curves[i]->curve_in_other_dom->add_pt(int_pts2[j], 1));
          
          for (k = 0; k < num_curves_split; k++)
            if (curves_split[k]->contains(*int_pts1[j]))
            //  Add int_pts1[j] to some curves_split.
            {
              curves_split[k]->add_pt(int_pts1[j]);
              curves_split[k + num_curves_split]->add_pt(int_pts1[j]);
            }
        }
        
        this_curve = int_curves[i];
        this_curve->ref_count++;
        
        for (j = 0; j < num_int_pts1; j++)
        {
          //  Locate int_pts[i_p1] on this_curve->segments[i_s1].
          
          i_p1 = num_int_pts1;
          i_s1 = 0;
          
          while (i_p1 == num_int_pts1 && i_s1 < this_curve->num_segments - 1)
          {
            for (i_p1 = 0;
                 i_p1 < num_int_pts1
                 &&
                 !int_pts1[i_p1]->equal(*this_curve->segments[i_s1]->end);
                 i_p1++)
              ;
            
            if (i_p1 == num_int_pts1)
              i_s1++;
          }
          
          assert(i_p1 < num_int_pts1);
          assert(i_s1 < this_curve->num_segments);
          
          this_curve->
            split(i_s1,
                  *(new_this_int_curves[num_new_this_int_curves] =
                    new K_CURVE),
                  *(new_this_curve = new K_CURVE));
          new_this_int_curves[num_new_this_int_curves]->ref_count++;
          new_this_adj_int_surfs[num_new_this_int_curves]   = adj_int_surfs[i];
          new_this_adj_int_surfs[num_new_this_int_curves]->ref_count++;
          new_this_adj_int_patches[num_new_this_int_curves] =
            adj_int_patches[i];
          new_this_adj_int_patches[num_new_this_int_curves]->ref_count++;
//        num_new_this_int_curves++;  //  Do this later.
          
          //  Split the other domain.
          
          for (i_s2 = 0;
               i_s2 < this_curve->curve_in_other_dom->num_segments
                 &&
                 !this_curve->curve_in_other_dom->segments[i_s2]->end->
                 equiv(*int_pts2[i_p1]);
               i_s2++)
            ;
          
          for (i_c2 = 0;
               i_c2 < adj_int_patches[i]->num_int_curves
                 &&
                 this_curve->curve_in_other_dom !=
                 adj_int_patches[i]->int_curves[i_c2];
               i_c2++)
            ;
          
          assert(i_s2 < this_curve->curve_in_other_dom->num_segments);
          assert(i_c2 < adj_int_patches[i]->num_int_curves);
          
          new_other_int_curves =
            new K_CURVE* [adj_int_patches[i]->num_int_curves + 1];
          new_other_adj_int_surfs =
            new K_SURF* [adj_int_patches[i]->num_int_curves + 1];
          new_other_adj_int_patches =
            new K_PATCH* [adj_int_patches[i]->num_int_curves + 1];
          
          for (k = 0; k < i_c2; k++)
          {
            new_other_int_curves[k]      = adj_int_patches[i]->int_curves[k];
            new_other_int_curves[k]->ref_count++;
            new_other_adj_int_surfs[k]   =
              adj_int_patches[i]->adj_int_surfs[k];
            new_other_adj_int_surfs[k]->ref_count++;
            new_other_adj_int_patches[k] =
              adj_int_patches[i]->adj_int_patches[k];
            new_other_adj_int_patches[k]->ref_count++;
          }
          
          adj_int_patches[i]->int_curves[i_c2]->
            split(i_s2,
                  *(new_other_int_curves[i_c2] = new K_CURVE),
                  *(new_other_int_curves[i_c2 + 1] = new K_CURVE));
          new_other_int_curves[i_c2]->ref_count++;
          new_other_int_curves[i_c2 + 1]->ref_count++;
          new_other_adj_int_surfs[i_c2]  = new_other_adj_int_surfs[i_c2 + 1] =
            adj_int_patches[i]->adj_int_surfs[i_c2];
          new_other_adj_int_surfs[i_c2]->ref_count++;
          new_other_adj_int_surfs[i_c2 + 1]->ref_count++;
          new_other_adj_int_patches[i_c2] =
          new_other_adj_int_patches[i_c2 + 1] =
            adj_int_patches[i]->adj_int_patches[i_c2];
          new_other_adj_int_patches[i_c2]->ref_count++;
          new_other_adj_int_patches[i_c2 + 1]->ref_count++;
          
          for (k = i_c2 + 1; k < adj_int_patches[i]->num_int_curves; k++)
          {
            new_other_int_curves[k + 1]      =
              adj_int_patches[i]->int_curves[k];
            new_other_int_curves[k + 1]->ref_count++;
            new_other_adj_int_surfs[k + 1]   =
              adj_int_patches[i]->adj_int_surfs[k];
            new_other_adj_int_surfs[k + 1]->ref_count++;
            new_other_adj_int_patches[k + 1] =
              adj_int_patches[i]->adj_int_patches[k];
            new_other_adj_int_patches[k + 1]->ref_count++;
          }
          
          for (k = 0; k < adj_int_patches[i]->num_int_curves; k++)
          {
            if (!--adj_int_patches[i]->int_curves[k]->ref_count)
              delete adj_int_patches[i]->int_curves[k];
            
            if (!--adj_int_patches[i]->adj_int_surfs[k]->ref_count)
              delete adj_int_patches[i]->adj_int_surfs[k];
            
            if (!--adj_int_patches[i]->adj_int_patches[k]->ref_count)
              delete adj_int_patches[i]->adj_int_patches[k];
          }
          
          delete [] adj_int_patches[i]->int_curves;
          delete [] adj_int_patches[i]->adj_int_surfs;
          delete [] adj_int_patches[i]->adj_int_patches;
          
          adj_int_patches[i]->int_curves      = new_other_int_curves;
          adj_int_patches[i]->adj_int_surfs   = new_other_adj_int_surfs;
          adj_int_patches[i]->adj_int_patches = new_other_adj_int_patches;
          adj_int_patches[i]->num_int_curves++;
          adj_int_patches[i]->merge_curves();
          
          //  Associate curves in both domains.
          
          assert(this_curve->dir_in_other_dom);
          
          if (this_curve->dir_in_other_dom > 0)
          {
            new_this_int_curves[num_new_this_int_curves]->
              assoc(adj_int_patches[i]->int_curves[i_c2], 1);
            new_this_curve->assoc(adj_int_patches[i]->int_curves[i_c2 + 1], 1);
          }
          else  //  if (this_curve->dir_in_other_dom < 0)
          {
            new_this_int_curves[num_new_this_int_curves]->
              assoc(adj_int_patches[i]->int_curves[i_c2 + 1], - 1);
            new_this_curve->assoc(adj_int_patches[i]->int_curves[i_c2], - 1);
          }
          
          num_new_this_int_curves++;
          
          if (!--this_curve->ref_count)
            delete this_curve;
          
          this_curve = new_this_curve;
          this_curve->ref_count++;
        }
        
        new_this_int_curves[num_new_this_int_curves]      = this_curve;
        new_this_int_curves[num_new_this_int_curves]->ref_count++;
        new_this_adj_int_surfs[num_new_this_int_curves]   = adj_int_surfs[i];
        new_this_adj_int_surfs[num_new_this_int_curves]->ref_count++;
        new_this_adj_int_patches[num_new_this_int_curves] = adj_int_patches[i];
        new_this_adj_int_patches[num_new_this_int_curves]->ref_count++;
        num_new_this_int_curves++;
        
        for (j = 0; j < num_int_pts1; j++)
          if (!--int_pts1[j]->ref_count)
            delete int_pts1[j];
        
        delete [] int_pts1;  //  num_int_pts1 > 0 => int_pts1 != 0
        
        if (num_int_pts2 > 0)
        {
          for (j = 0; j < num_int_pts2; j++)
            if (!--int_pts2[j]->ref_count)
              delete int_pts2[j];
          
          delete [] int_pts2;
        }
        
        if (!--this_curve->ref_count)
          delete this_curve;
      }
      else  //  if (!num_int_pts1)
      {
        new_this_int_curves[num_new_this_int_curves]      = int_curves[i];
        new_this_adj_int_surfs[num_new_this_int_curves]   = adj_int_surfs[i];
        new_this_adj_int_patches[num_new_this_int_curves] = adj_int_patches[i];
        num_new_this_int_curves++;
      }
    }
    
    //  Put intersection curves on patches.
    
    for (i = 0; i < num_new_this_int_curves; i++)
    {
      //  Find 1 & only 1 sub_patches[i_p1]
      //  on which new_this_int_curves[i] lies.
      
      x = new_this_int_curves[i]->pt_on();
      
      for (i_p1 = 0;
           i_p1 < num_sub_patches && !sub_patches[i_p1]->contains(x);
           i_p1++)
        ;
      
      assert(i_p1 < num_sub_patches);
      
      num_this_int_curves  = sub_patches[i_p1]->num_int_curves + 1;
      this_int_curves      = new K_CURVE* [num_this_int_curves];
      this_adj_int_surfs   = new K_SURF* [num_this_int_curves];
      this_adj_int_patches = new K_PATCH* [num_this_int_curves];
      
      for (j = 0; j < num_this_int_curves - 1; j++)
      {
        this_int_curves[j]       = sub_patches[i_p1]->int_curves[j];
        this_int_curves[j]->ref_count++;
        this_adj_int_surfs[j]    = sub_patches[i_p1]->adj_int_surfs[j];
        this_adj_int_surfs[j]->ref_count++;
        this_adj_int_patches[j]  = sub_patches[i_p1]->adj_int_patches[j];
        this_adj_int_patches[j]->ref_count++;
      }
      
      this_int_curves[num_this_int_curves - 1]      =
        new_this_int_curves[i];
      this_int_curves[num_this_int_curves - 1]->ref_count++;
      this_adj_int_surfs[num_this_int_curves - 1]   =
        new_this_adj_int_surfs[i];
      this_adj_int_surfs[num_this_int_curves - 1]->ref_count++;
      this_adj_int_patches[num_this_int_curves - 1] =
        new_this_adj_int_patches[i];
      this_adj_int_patches[num_this_int_curves - 1]->ref_count++;
      
      if (sub_patches[i_p1]->num_int_curves > 0)
      {
        for (j = 0; j < sub_patches[i_p1]->num_int_curves; j++)
        {
          if (!--sub_patches[i_p1]->int_curves[j]->ref_count)
            delete sub_patches[i_p1]->int_curves[j];
          
          if (!--sub_patches[i_p1]->adj_int_surfs[j]->ref_count)
            delete sub_patches[i_p1]->adj_int_surfs[j];
          
          if (!--sub_patches[i_p1]->adj_int_patches[j]->ref_count)
            delete sub_patches[i_p1]->adj_int_patches[j];
        }
        
        delete [] sub_patches[i_p1]->int_curves;
        delete [] sub_patches[i_p1]->adj_int_surfs;
        delete [] sub_patches[i_p1]->adj_int_patches;
      }
      
      sub_patches[i_p1]->int_curves      = this_int_curves;
      sub_patches[i_p1]->adj_int_surfs   = this_adj_int_surfs;
      sub_patches[i_p1]->adj_int_patches = this_adj_int_patches;
      sub_patches[i_p1]->num_int_curves  = num_this_int_curves;
      
      //  Update adj_int_patches of adj_int_patches.
      
      if (new_this_adj_int_patches[i]->num_int_curves > 0)
      {
        for (i_c2 = 0;
             i_c2 < new_this_adj_int_patches[i]->num_int_curves
             &&
             new_this_adj_int_patches[i]->int_curves[i_c2] !=
               new_this_int_curves[i]->curve_in_other_dom;
             i_c2++)
          ;
        
        assert(i_c2 < new_this_adj_int_patches[i]->num_int_curves);
        
        if (!--new_this_adj_int_patches[i]->adj_int_patches[i_c2]->ref_count)
          delete new_this_adj_int_patches[i]->adj_int_patches[i_c2];
        
        new_this_adj_int_patches[i]->adj_int_patches[i_c2] = sub_patches[i_p1];
        new_this_adj_int_patches[i]->adj_int_patches[i_c2]->ref_count++;
      }
    }
    
    for (i = 0; i < num_pts_split; i++)
      if (!--curves_split[i]->ref_count)
        delete curves_split[i];
    
    delete [] curves_split;
    
    for (i = 0; i < num_new_this_int_curves; i++)
    {
      if (!--new_this_int_curves[i]->ref_count)
        delete new_this_int_curves[i];
      
      if (!--new_this_adj_int_surfs[i]->ref_count)
        delete new_this_adj_int_surfs[i];
      
      if (!--new_this_adj_int_patches[i]->ref_count)
        delete new_this_adj_int_patches[i];
    }
    
    delete [] new_this_int_curves;       //  MAX_NUM_PTS_SPLIT > 0
    delete [] new_this_adj_int_surfs;    //  MAX_NUM_PTS_SPLIT > 0
    delete [] new_this_adj_int_patches;  //  MAX_NUM_PTS_SPLIT > 0
    
    //  Merge curves on sub_patches and recurse.
    
    num_new_patches_proto = new unsigned long [num_sub_patches];
    new_patches_proto     = new K_PATCH** [num_sub_patches];
    
    for (i = 0; i < num_sub_patches; i++)
    {
      sub_patches[i]->merge_curves();
      num_new_patches_proto[i] =
        sub_patches[i]->split_loops(new_patches_proto[i]);
    }
    
    num_new_patches = 0;
    
    for (i = 0; i < num_sub_patches; i++)
      num_new_patches += num_new_patches_proto[i];
    
    assert(num_new_patches > 0);
    
    new_patches = new K_PATCH* [num_new_patches];
    k           = 0;
    
    for (i = 0; i < num_sub_patches; i++)
      for (j = 0; j < num_new_patches_proto[i]; j++)
      {
        new_patches[k] = new_patches_proto[i][j];
        new_patches[k]->ref_count++;
        k++;
      }
    
    for (i = 0; i < num_sub_patches; i++)
    {
      for (j = 0; j < num_new_patches_proto[i]; j++)
        if (!--new_patches_proto[i][j]->ref_count)
          delete new_patches_proto[i][j];
      
      delete [] new_patches_proto[i];
    }
    
    delete [] new_patches_proto;      //  num_sub_patches >= 1
    delete [] num_new_patches_proto;  //  num_sub_patches >= 1
    
    for (i = 0; i < num_sub_patches; i++)
      if (!--sub_patches[i]->ref_count)
        delete sub_patches[i];
    
    delete [] sub_patches;  //  MAX_NUM_PTS_SPLIT > 0
  }
  else  //  if (loop == num_merged)
  {
    new_patches    = new K_PATCH* [num_new_patches = 1];
    new_patches[0] = this;
    new_patches[0]->ref_count++;
  }
  
  return num_new_patches;
}

#define MAX_NUM_SEGMENTS 1024

#ifdef CCW_OUTPUT
int K_PATCH :: Bezier_output(ostream& out_fs, const int is_fw0) const
{
  unsigned long i, j;
  unsigned long num_segments;
  double*       a[MAX_NUM_SEGMENTS];
  double        max_a_s, min_a_s, max_a_t;
  unsigned long max_s_segment, min_s_segment, max_t_segment;
  int           is_fw;
  double        d_l[2], d_h[2];
  
  surf->Bezier_output(out_fs, low_s, high_s, low_t, high_t);
  
//  out_fs << "172  172  172" << endl << flush;  //  color: grey
  out_fs << "0  127  127" << endl << flush;  //  color: green
  out_fs << "1" << endl << flush;              //  num_trim_loops: 1
  
  for (i = 0; i < MAX_NUM_SEGMENTS; i++)
    a[i] = new double [2];
  
  max_a_s       = - 1.0;
  min_a_s       = 10.0;
  max_a_t       = - 1.0;
  max_s_segment = min_s_segment = max_t_segment = 0;
  num_segments  = 0;
  
  for (i = 0; i < num_trim_curves; i++)
  {
    for (j = 0; j < 2; j++)
      trim_curves[i]->sub_divide(10, j);
    
    for (j = 0; j < trim_curves[i]->num_segments; j++)
    {
      assert(num_segments < MAX_NUM_SEGMENTS);
      trim_curves[i]->segments[j]->start->get_fp_approx(a[num_segments]);
      
      if (a[num_segments][0] > max_a_s)
      {
        max_a_s       = a[num_segments][0];
        max_t_segment = num_segments;
      }
      
      if (a[num_segments][0] < min_a_s)
      {
        min_a_s       = a[num_segments][0];
        min_s_segment = num_segments;
      }
      
      if (a[num_segments][1] > max_a_t)
      {
        max_a_t       = a[num_segments][1];
        max_t_segment = num_segments;
      }
      
      num_segments++;
    }   
  }
  
  assert(num_segments < MAX_NUM_SEGMENTS);
  trim_curves[0]->start()->get_fp_approx(a[num_segments]);
  
  assert(max_s_segment != min_s_segment
         &&
         max_s_segment != max_t_segment
         &&
         min_s_segment != max_t_segment);
  
  //  Orient CCW s.t. we visit segments in the following order:
  //    max_s_segment => max_t_segment => min_s_segment
  
  if (max_t_segment > max_s_segment)
    if (max_s_segment > min_s_segment || min_s_segment > max_t_segment)
      is_fw = 1;
    else
      is_fw = 0;
  else  //  if (max_s_segment > max_t_segment)
    if (max_t_segment > min_s_segment || min_s_segment > max_s_segment)
      is_fw = 1;
    else
      is_fw = 0;
  
  d_l[0] = low_s.as_double();
  d_h[0] = high_s.as_double();
  d_l[1] = low_t.as_double();
  d_h[1] = high_t.as_double();
  
  for (i = 0; i < num_segments + 1; i++)
    for (j = 0; j < 2; j++)
      a[i][j] = (a[i][j] - d_l[j]) / (d_h[j] - d_l[j]);
  
  out_fs << num_segments << endl << flush;
  
  if (is_fw)
    for (i = 0; i < num_segments; i++)
    {
      out_fs << "-1  -1  1  1" << endl << flush;
      out_fs << a[i][0] << "  " << a[i][1] << endl << flush;
      out_fs << a[i + 1][0] << "  " << a[i + 1][1] << endl << flush;
    }
  else  //  if (!is_fw)
    for (i = num_segments; i > 0; i--)
    {
      out_fs << "-1  -1  1  1" << endl << flush;
      out_fs << a[i][0] << "  " << a[i][1] << endl << flush;
      out_fs << a[i - 1][0] << "  " << a[i - 1][1] << endl << flush;
    }
  
  for (i = 0; i < MAX_NUM_SEGMENTS; i++)
    delete [] a[i];
  
  return 0;
}
#else
int K_PATCH :: Bezier_output(ostream& out_fs, const int is_fw0) const
{
  unsigned long i, j;
  unsigned long num_segments;
  double*       a[MAX_NUM_SEGMENTS];
  int           is_fw;
  double        d_l[2], d_h[2];
  
  surf->Bezier_output(out_fs, low_s, high_s, low_t, high_t);
  
//  out_fs << "172  172  172" << endl << flush;  //  color: grey
  out_fs << "0  127  127" << endl << flush;  //  color: green
  out_fs << "1" << endl << flush;              //  num_trim_loops: 1
  
  is_fw = is_fw0;
  
  for (i = 0; i < MAX_NUM_SEGMENTS; i++)
    a[i] = new double [2];
  
  num_segments = 0;
  
  for (i = 0; i < num_trim_curves; i++)
  {
    if (!is_fw)  //  s.t. sub_divide won't be applied twice.
      for (j = 0; j < 2; j++)
        trim_curves[i]->sub_divide(10, j);
    
    for (j = 0; j < trim_curves[i]->num_segments; j++)
    {
      assert(num_segments < MAX_NUM_SEGMENTS);
      trim_curves[i]->segments[j]->start->get_fp_approx(a[num_segments]);
      num_segments++;
    }   
  }
  
  assert(num_segments < MAX_NUM_SEGMENTS);
  trim_curves[0]->start()->get_fp_approx(a[num_segments]);
  
  d_l[0] = low_s.as_double();
  d_h[0] = high_s.as_double();
  d_l[1] = low_t.as_double();
  d_h[1] = high_t.as_double();
  
  for (i = 0; i < num_segments + 1; i++)
    for (j = 0; j < 2; j++)
      a[i][j] = (a[i][j] - d_l[j]) / (d_h[j] - d_l[j]);
  
  out_fs << num_segments << endl << flush;
  
  if (is_fw)
    for (i = 0; i < num_segments; i++)
    {
      out_fs << "-1 -1 1 1" << endl << flush;
      out_fs << a[i][0] << " " << a[i][1] << endl << flush;
      out_fs << a[i + 1][0] << " " << a[i + 1][1] << endl << flush;
    }
  else  //  if (!is_fw)
    for (i = num_segments; i > 0; i--)
    {
      out_fs << "-1 -1 1 1" << endl << flush;
      out_fs << a[i][0] << " " << a[i][1] << endl << flush;
      out_fs << a[i - 1][0] << " " << a[i - 1][1] << endl << flush;
    }
  
  for (i = 0; i < MAX_NUM_SEGMENTS; i++)
    delete [] a[i];
  
  return 0;
}
#endif

