#include <pt_surf_assoc.h>

PT_SURF_ASSOC :: PT_SURF_ASSOC()
{
  len   = 0;
  pts   = 0;
  surfs = 0;
  rep   = 0;
}

PT_SURF_ASSOC :: ~PT_SURF_ASSOC()
{
  unsigned long i;

  if (len > 0)
  {
    for (i = 0; i < len; i++)
    {
      if (!--pts[i]->ref_count)
        delete pts[i];

      if (!--surfs[i]->ref_count)
        delete surfs[i];
    }

    delete [] pts;
    delete [] surfs;
    delete [] rep;
  }
}

unsigned long PT_SURF_ASSOC :: find(const K_POINT2D* x, const K_SURF* y)
{
  unsigned long i;

  for (i = 0; i < len && (x != pts[i] || y != surfs[i]); i++)
    ;

  return i;
}

unsigned long PT_SURF_ASSOC :: locate(K_POINT2D* const x, K_SURF* const y)
{
  unsigned long i, j;
  K_POINT2D**   p;
  K_SURF**      s;
  long*         r;

  if ((i = find(x, y)) == len)
  {
    p = new K_POINT2D* [len + 1];
    s = new K_SURF* [len + 1];
    r = new long [len + 1];

    for (j = 0; j < len; j++)
    {
      p[j] = pts[j];
//      p[j]->ref_count++;
//
//      if (!--pts[j]->ref_count)
//        delete pts[j];
//
      s[j] = surfs[j];
//      s[j]->ref_count++;
//
//      if (!--surfs[j]->ref_count)
//        delete surfs[j];
//
      r[j] = rep[j];
    }

    p[len] = x;
    p[len]->ref_count++;
    s[len] = y;
    s[len]->ref_count++;
    r[len] = len;
    delete [] pts;
    delete [] surfs;
    delete [] rep;
    len++;
    pts    = p;
    surfs  = s;
    rep    = r;
  }

  return i;
}

unsigned long PT_SURF_ASSOC :: record_as_assoc(K_POINT2D* const x,
                                               K_SURF* const y,
                                               K_POINT2D* const z,
                                               K_SURF* const w)
{
  unsigned long i, j, k, q;
  unsigned long r;

  i = locate(x, y);
  j = locate(z, w);

  if (rep[i] < rep[j])
  {
    r = rep[i];
    q = rep[j];
  }
  else  //  if (rep[i] >= rep[j])
  {
    r = rep[j];
    q = rep[i];
  }

  if (r != q)
    for (k = 0; k < len; k++)
      if (rep[k] == q)
        rep[k] = r;

  return r;
}

K_POINT2D* PT_SURF_ASSOC :: find_assoc_pt(const K_POINT2D* x0, const K_SURF* y,
                                          const K_SURF* z)
{
  unsigned long i, j;
  K_POINT2D*    x;
  K_POINT2D*    w;

  w = 0;

  if ((i = find(x0, y)) < len)
  {
    for (j = 0; j < len && (rep[j] != rep[i] || surfs[j] != z); j++)
      ;

    if (j < len)
      w = pts[j];
  }

  for (x = x0->next; !w && x != x0; x = x->next)
  {
    if ((i = find(x, y)) < len)
      for (j = 0; j < len && (rep[j] != rep[i] || surfs[j] != z); j++)
        ;

    if (j < len)
      w = pts[j];
  }

  return w;
}

