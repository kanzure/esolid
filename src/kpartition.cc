//  file:    kpartition.cc
//  update:  01/15/03

#include <kpartition.h>
#include <kgraph.h>

static unsigned long PARTITION_ID_COUNT = 0;

K_PARTITION :: K_PARTITION()
  : ID(PARTITION_ID_COUNT++),
    from(0),
    num_trim_curves(0),
    trim_curves(0), adj_curves(0),
    adj_surfs(0),
    adj_partitions(0), rev_curves(0),
    ref_count(0)
{ }

K_PARTITION :: K_PARTITION(K_PATCH* const p)
  : ID(PARTITION_ID_COUNT++)
{
  unsigned long i;
  
  if (from = p)
    from->ref_count++;
  
  if ((num_trim_curves = p->num_trim_curves) > 0)
  {
    trim_curves    = new K_CURVE* [num_trim_curves];
    adj_curves     = new long [num_trim_curves];
    adj_surfs      = new K_SURF* [num_trim_curves];
    adj_partitions = new K_PARTITION* [num_trim_curves];
    rev_curves     = new int [num_trim_curves];
    
    for (i = 0; i < num_trim_curves; i++)
    {
      trim_curves[i]    = p->trim_curves[i];
      trim_curves[i]->ref_count++;
      adj_curves[i]     = i;
      adj_surfs[i]      = 0;
      adj_partitions[i] = 0;
      rev_curves[i]     = 0;
    }
  }
  else  //  if (!num_trim_curves)
  {
    trim_curves    = 0;
    adj_curves     = 0;
    adj_surfs      = 0;
    adj_partitions = 0;
    rev_curves     = 0;
  }
  
  ref_count = 0;
}

K_PARTITION :: K_PARTITION(const K_PARTITION& p)
  : ID(PARTITION_ID_COUNT++)
{
  unsigned long i;
  
  if (from = p.from)
    from->ref_count++;
  
  if ((num_trim_curves = p.num_trim_curves) > 0)
  {
    trim_curves    = new K_CURVE* [num_trim_curves];
    adj_curves     = new long [num_trim_curves];
    adj_surfs      = new K_SURF* [num_trim_curves];
    adj_partitions = new K_PARTITION* [num_trim_curves];
    rev_curves     = new int [num_trim_curves];
    
    for (i = 0; i < num_trim_curves; i++)
    {
      if (trim_curves[i] = p.trim_curves[i])
        trim_curves[i]->ref_count++;
      
      adj_curves[i] = p.adj_curves[i];
      
      if (adj_surfs[i] = p.adj_surfs[i])
        adj_surfs[i]->ref_count++;
      
      if (adj_partitions[i] = p.adj_partitions[i])
        adj_partitions[i]->ref_count++;
      
      rev_curves[i] = p.rev_curves[i];
    }
  }
  else  //  if (!num_trim_curves)
  {
    trim_curves    = 0;
    adj_curves     = 0;
    adj_surfs      = 0;
    adj_partitions = 0;
    rev_curves     = 0;
  }
  
  ref_count = 0;
}

K_PARTITION& K_PARTITION :: operator =(const K_PARTITION& p)
{
  if (this != &p)
  {
    unsigned long i;
    
    if (from && !--from->ref_count)
      delete from;
    
    if (num_trim_curves > 0)
    {
      for (i = 0; i < num_trim_curves; i++)
        if (!--trim_curves[i]->ref_count)
          delete trim_curves[i];
      
      delete [] trim_curves;
      delete [] adj_curves;
      
      for (i = 0; i < num_trim_curves; i++)
        if (adj_surfs[i] && !--adj_surfs[i]->ref_count)
          delete adj_surfs[i];
      
      delete [] adj_surfs;
      
      for (i = 0; i < num_trim_curves; i++)
        if (adj_partitions[i] && !--adj_partitions[i]->ref_count)
          delete adj_partitions[i];
      
      delete [] adj_partitions;
      delete [] rev_curves;
    }
    
    if (from = p.from)
      from->ref_count++;
    
    if ((num_trim_curves = p.num_trim_curves) > 0)
    {
      trim_curves    = new K_CURVE* [num_trim_curves];
      adj_curves     = new long [num_trim_curves];
      adj_surfs      = new K_SURF* [num_trim_curves];
      adj_partitions = new K_PARTITION* [num_trim_curves];
      rev_curves     = new int [num_trim_curves];
      
      for (i = 0; i < num_trim_curves; i++)
      {
        if (trim_curves[i] = p.trim_curves[i])
          trim_curves[i]->ref_count++;
        
        adj_curves[i] = p.adj_curves[i];
        
        if (adj_surfs[i] = p.adj_surfs[i])
          adj_surfs[i]->ref_count++;
        
        if (adj_partitions[i] = p.adj_partitions[i])
          adj_partitions[i]->ref_count++;
        
        rev_curves[i] = p.rev_curves[i];
      }
    }
    else  //  if (!num_trim_curves)
    {
      trim_curves    = 0;
      adj_curves     = 0;
      adj_surfs      = 0;
      adj_partitions = 0;
      rev_curves     = 0;
    }
  }
  
  return *this;
}

K_PARTITION :: ~K_PARTITION()
{
  unsigned long i;
  
  if (from && !--from->ref_count)
    delete from;
  
  if (num_trim_curves > 0)
  {
    for (i = 0; i < num_trim_curves; i++)
      if (!--trim_curves[i]->ref_count)
        delete trim_curves[i];
    
    delete [] trim_curves;
    delete [] adj_curves;
    
    for (i = 0; i < num_trim_curves; i++)
      if (adj_surfs[i] && !--adj_surfs[i]->ref_count)
        delete adj_surfs[i];
    
    delete [] adj_surfs;
    
    for (i = 0; i < num_trim_curves; i++)
      if (adj_partitions[i] && !--adj_partitions[i]->ref_count)
        delete adj_partitions[i];
    
    delete [] adj_partitions;
    delete [] rev_curves;
  }
}

ostream& K_PARTITION :: output(ostream& o) const
{
  o << " partition ID: " << ID << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const K_PARTITION& p)
{
  return p.output(o);
}

unsigned long gen_partitions(K_PATCH* const patch, K_PARTITION**& partitions)
{
  long           i, j, k;
  unsigned long  num_curves;
  K_POINT2D**    int_end_pts;
  long*          partitions_fw;
  long*          partitions_rev;  
  long**         list_trim;
  long**         list_int;
  unsigned long* num_curves_on_partition;
  int*           trim_used;
  int            done, go_on;
  unsigned long  i_s, i_c;
  K_POINT2D*     end_pt;
  unsigned long  num_partitions;
  
  if ((num_curves = patch->num_merged) > 0)
  {
    int_end_pts = new K_POINT2D* [2 * patch->num_int_curves];
    
    for (i = 0; i < num_curves; i++)
    {
      int_end_pts[2 * i]     = patch->ic[i][0]->start();
      int_end_pts[2 * i]->ref_count++;
      int_end_pts[2 * i + 1] = patch->ic[i][patch->len_ic[i] - 1]->end();
      int_end_pts[2 * i + 1]->ref_count++;
    }
    
    partitions_fw  = new long [num_curves];
    partitions_rev = new long [num_curves];
    
    for (i = 0; i < num_curves; i++)
    {
      partitions_fw[i]  = - 1;
      partitions_rev[i] = - 1;
    }
    
    list_trim               = new long* [num_curves + 1];
    list_int                = new long* [num_curves + 1];
    num_curves_on_partition = new unsigned long [num_curves + 1];
    
    for (i = 0; i < num_curves + 1; i++)
    {
      list_trim[i]               =
        new long [patch->num_trim_curves + num_curves];
      list_int[i]                =
        new long [patch->num_trim_curves + num_curves];
    }
    
    for (i = 0; i < num_curves + 1; i++)
    {
      for (j = 0; j < patch->num_trim_curves + num_curves; j++)
        list_trim[i][j] = list_int[i][j] = - 1;
      
      num_curves_on_partition[i] = 0;
    }
    
    trim_used = new int [patch->num_trim_curves];
    
    for (i = 0; i < patch->num_trim_curves; i++)
      trim_used[i] = 0;
    
    num_partitions = 0;
    done           = 0;
    
    while (!done)
    {
      for (i = 0; i < patch->num_trim_curves && trim_used[i]; i++)
        ;
      
      if (i < patch->num_trim_curves)
      {
        i_s          = i;
        go_on        = 1;
        
        while (go_on)
        {
          list_trim[num_partitions]
                   [num_curves_on_partition[num_partitions]] = i;
          trim_used[i] = 1;
          list_int[num_partitions]
                  [num_curves_on_partition[num_partitions]]  = - 1;
          num_curves_on_partition[num_partitions]++;
          
          for (j = 0;
               j < 2 * num_curves
               &&
               !int_end_pts[j]->equal(*patch->trim_curves[i]->end());
               j++)
            ;
          
          if (j < 2 * num_curves)
          //  The next curve is an intersection curve.
          {
            list_trim[num_partitions]
                     [num_curves_on_partition[num_partitions]] = - 1;
            list_int[num_partitions]
                    [num_curves_on_partition[num_partitions]]  = j / 2;
            num_curves_on_partition[num_partitions]++;
            
            if (j % 2)
            {
              end_pt                = patch->ic[j / 2][0]->start();
              partitions_rev[j / 2] = num_partitions;
            }
            else  //  if (!(j % 2))
            {
              end_pt               =
                patch->ic[j / 2][patch->len_ic[j / 2] - 1]->end();
              partitions_fw[j / 2] = num_partitions;
            }
            
            for (i = 0;
                 i < patch->num_trim_curves
                 &&
                 !end_pt->equal(*patch->trim_curves[i]->start());
                 i++)
              ;
            
            assert(i < patch->num_trim_curves);
          }
          else  //  if (j == 2 * num_curves)
          //  The next curve is a trimming curve.
            i = (i + 1) % patch->num_trim_curves;
          
          if (i == i_s)
            go_on = 0;
        }
        
        num_partitions++;
      }
      else  //  if (i == patch->num_trim_curves)
        done = 1;
    }
    
    assert(num_partitions > 0);
    
    partitions = new K_PARTITION* [num_partitions];
    
    for (i = 0; i < num_partitions; i++)
    {
      partitions[i] = new K_PARTITION;
      partitions[i]->ref_count++;
      
      partitions[i]->from = patch;
      partitions[i]->from->ref_count++;
      
      partitions[i]->num_trim_curves = 0;
      
      for (j = 0; j < num_curves_on_partition[i]; j++)
        if (list_trim[i][j] >= 0)
          partitions[i]->num_trim_curves++;
        else  //  if (list_int[i][j] >= 0)
          partitions[i]->num_trim_curves += patch->len_ic[list_int[i][j]];
      
      partitions[i]->trim_curves    =
        new K_CURVE* [partitions[i]->num_trim_curves];
      partitions[i]->adj_curves     =
        new long [partitions[i]->num_trim_curves];
      partitions[i]->adj_surfs      =
        new K_SURF* [partitions[i]->num_trim_curves];
      partitions[i]->adj_partitions =
        new K_PARTITION* [partitions[i]->num_trim_curves];
      partitions[i]->rev_curves     = new int [partitions[i]->num_trim_curves];
    }
    
    for (i = 0; i < num_partitions; i++)
    {
      i_c = 0;
      
      for (j = 0; j < num_curves_on_partition[i]; j++)
        if (list_trim[i][j] >= 0)
        {
          partitions[i]->trim_curves[i_c]    =
            patch->trim_curves[list_trim[i][j]];
          partitions[i]->trim_curves[i_c]->ref_count++;
          partitions[i]->adj_curves[i_c]     = list_trim[i][j];
          partitions[i]->adj_surfs[i_c]      = 0;
          partitions[i]->adj_partitions[i_c] = 0;
          partitions[i]->rev_curves[i_c]     = 0;
          i_c++;
        }
        else  //  if (list_trim[i][j] < 0 && list_int[i][j] >= 0)
          if (i == partitions_fw[list_int[i][j]])
            for (k = 0; k < patch->len_ic[list_int[i][j]]; k++)
            {
              partitions[i]->trim_curves[i_c]    =
                patch->ic[list_int[i][j]][k];
              partitions[i]->trim_curves[i_c]->ref_count++;
              partitions[i]->adj_curves[i_c]     = - 1;
              partitions[i]->adj_surfs[i_c]      = 0;
              partitions[i]->adj_partitions[i_c] =
                partitions[partitions_rev[list_int[i][j]]];
              partitions[i]->adj_partitions[i_c]->ref_count++;
              partitions[i]->rev_curves[i_c]     = 0;
              i_c++;
            }
          else  //  if (i != partitions_fw[list_int[i][j]])
            for (k = patch->len_ic[list_int[i][j]] - 1; k >= 0; k--)
            {
              partitions[i]->trim_curves[i_c]    =
                patch->ic[list_int[i][j]][k];
              partitions[i]->trim_curves[i_c]->ref_count++;
              partitions[i]->adj_curves[i_c]     = - 1;
              partitions[i]->adj_surfs[i_c]      = 0;
              partitions[i]->adj_partitions[i_c] =
                partitions[partitions_fw[list_int[i][j]]];
              partitions[i]->adj_partitions[i_c]->ref_count++;
              partitions[i]->rev_curves[i_c]     = 1;
              i_c++;
            }
    }
    
    //  Delete patch's old informations.
    
    for (i = 0; i < patch->num_merged; i++)
    {
      for (j = 0; j < patch->len_ic[i]; j++)
        if (!--patch->ic[i][j]->ref_count)
          delete patch->ic[i][j];
      
      if (patch->len_ic[i] > 0)
        delete [] patch->ic[i];
    }
    
    delete [] patch->ic;
    delete [] patch->len_ic;
    delete [] patch->ic_closed;
    
    patch->num_merged = 0;
    
    for (i = 0; i < 2 * num_curves; i++)
      if (!--int_end_pts[i]->ref_count)
        delete int_end_pts[i];
    
    delete [] int_end_pts;
    
    delete [] partitions_fw;
    delete [] partitions_rev;
    
    delete [] num_curves_on_partition;
    
    for (i = 0; i < num_curves + 1; i++)
    {
      delete [] list_trim[i];
      delete [] list_int[i];
    }
    
    delete [] list_trim;
    delete [] list_int;
    
    delete [] trim_used;
  }
  else  //  if (!num_curves)
  {
    partitions    = new K_PARTITION* [num_partitions = 1];
    partitions[0] = new K_PARTITION(patch);
    partitions[0]->ref_count++;
  }
  
  return num_partitions;
}

unsigned long gen_adjacency(K_PARTITION** const partitions,
                            const unsigned long num_partitions,
                            K_GRAPH*& G,
                            long*& C, K_GRAPH*& CG)
{
  assert(num_partitions > 0);
  
  unsigned long  i, j, k, l;
  unsigned long* V;
  unsigned long  f, t;
  K_PATCH*       p;
  long           c;
  K_PARTITION*   a;
  
  V = new unsigned long [num_partitions];
  
  for (i = 0; i < num_partitions; i++)
    V[i] = partitions[i]->ID;
  
  G = new K_GRAPH(num_partitions, V);
  
  for (i = 0; i < num_partitions; i++)
  {
    f = partitions[i]->ID;
    
    for (j = 0; j < partitions[i]->num_trim_curves; j++)
      if (a = partitions[i]->adj_partitions[j])
      {
        t = a->ID;
        G->add_edge(f, t, 1);
      }
      else  //  if (!partitions[i]->adj_partitions[j])
      {
        p = partitions[i]->from;
        c = partitions[i]->adj_curves[j];
        a = 0;
        
        for (k = 0; !a && k < num_partitions; k++)
          if (partitions[k]->from == p->adj_patches[c])
          {
            for (l = 0;
                 l < partitions[k]->num_trim_curves
                 &&
                 partitions[k]->trim_curves[l] !=
                   p->trim_curves[c]->curve_in_other_dom;
                 l++)
              ;
            
            if (l < partitions[k]->num_trim_curves)
              a = partitions[k];
          }
        
        if (a)
        {
          t = a->ID;
          G->add_edge(f, t, 0);
        }
      }
  }
  
  delete [] V;
  
  return G->get_connected_components(C, CG);
}

#define MAX_NUM_INT_PTS 128

K_POINT2D K_PARTITION :: get_pt_in() const
{
  assert(num_trim_curves > 0);
  
  unsigned long i, j, k;
  K_BOXCO2      b;
  bigrational   low_t, high_t, cut_t;
  K_RATPOLY     poly_cut;
  K_POINT2D**   int_pts_proto;
  unsigned long num_int_pts_proto;
  K_POINT2D**   int_pts;
  unsigned long num_int_pts;
  bigrational   v_s;
  
  b      = trim_curves[0]->bbox();
  low_t  = b.low[1];
  high_t = b.high[1];
  
  for (i = 1; i < num_trim_curves; i++)
  {
    b = trim_curves[i]->bbox();
    
    if (b.low[1] < low_t)
      low_t = b.low[1];
    
    if (b.high[1] > high_t)
      high_t = b.high[1];
  }
  
  cut_t       = low_t + 3 * (high_t - low_t) / 7;
  poly_cut    = K_RATPOLY(2, 1, cut_t);
  int_pts     = new K_POINT2D* [MAX_NUM_INT_PTS];
  num_int_pts = 0;
  
  for (i = 0; i < num_trim_curves; i++)
  {
    if (num_int_pts_proto =
        trim_curves[i]->find_intersections(poly_cut, int_pts_proto, 1) > 0)
    {
      for (j = 0; j < num_int_pts_proto; j++)
      {
        assert(num_int_pts < MAX_NUM_INT_PTS);
        int_pts[num_int_pts] = int_pts_proto[j];
        int_pts[num_int_pts]->ref_count++;
        num_int_pts++;
      }
      
      for (j = 0; j < num_int_pts_proto; j++)
        if (!--int_pts_proto[j]->ref_count)
          delete int_pts_proto[j];
      
      delete [] int_pts_proto;
    }
  }
  
  assert(num_int_pts >= 2);
  
  sort_s(int_pts, num_int_pts);
  v_s = (int_pts[0]->get_high_s() + int_pts[1]->get_low_s()) / 2;
  
  for (i = 0; i < num_int_pts; i++)
    if (!--int_pts[i]->ref_count)
      delete int_pts[i];
  
  delete [] int_pts;
  
  return K_POINT2D(v_s, cut_t);
}

unsigned long select_relevant_partitions(K_PARTITION** const partitions,
                                         const unsigned long num_partitions,
                                         K_GRAPH* const G,
                                         const long* const components,
                                         const int* const colors,
                                         const int color,
                                         K_PARTITION**& relevant_partitions,
                                         K_GRAPH*& SG)
{
  unsigned long  i;
  unsigned long* V_proto;
  K_PARTITION**  relevant_partitions_proto;
  unsigned long* V;
  unsigned long  num_relevant_partitions;
  
  V_proto                   = new unsigned long [num_partitions];
  relevant_partitions_proto = new K_PARTITION* [num_partitions];
  
  num_relevant_partitions = 0;
  
  for (i = 0; i < num_partitions; i++)
    if (colors[components[i]] == color)
    {
      V_proto[num_relevant_partitions] = partitions[i]->ID;
      relevant_partitions_proto[num_relevant_partitions] = partitions[i];
//      cerr << " kpartition: s_r_p: proto[" << num_relevant_partitions << "] => [" << i << "] " << endl << flush;
      relevant_partitions_proto[num_relevant_partitions]->ref_count++;
      num_relevant_partitions++;
    }
  
  V                   = new unsigned long [num_relevant_partitions];
  relevant_partitions = new K_PARTITION* [num_relevant_partitions];
  
  for (i = 0; i < num_relevant_partitions; i++)
  {
    V[i]                   = V_proto[i];
    relevant_partitions[i] = relevant_partitions_proto[i];
    relevant_partitions[i]->ref_count++;
  }
  
  for (i = 0; i < num_relevant_partitions; i++)
    if (!--relevant_partitions_proto[i]->ref_count)
      delete relevant_partitions_proto[i];
  
  delete [] V_proto;
  delete [] relevant_partitions_proto;
  
  SG = new K_GRAPH(G->gen_subgraph(num_relevant_partitions, V));
  
  delete [] V;
  
  return num_relevant_partitions;
}

