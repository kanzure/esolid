#include <config.h>
#ifdef _EXPERIMENT
#include <counter.h>
#endif

#include <ksolid.h>
#include <genbox.h>
#include <gencyl.h>
#include <genell.h>
#include <gentor.h>
#include <fpconversion.h>
#include <cstring>
#include <fstream>

K_SOLID :: K_SOLID()
{
  num_patches = 0;
  patches     = 0;
}

K_SOLID :: K_SOLID(K_PATCH* const p[], const unsigned long n)
{
  unsigned long i;
  
  if ((num_patches = n) > 0)
  {
    patches = new K_PATCH* [num_patches];
    
    for (i = 0; i < num_patches; i++)
    {
      patches[i] = p[i];
      patches[i]->ref_count++;
    }
  }
  else  //  if (!num_patches)
    patches = 0;
}

K_SOLID :: K_SOLID(const K_SOLID& s)
{
  unsigned long i;
  
  if ((num_patches = s.num_patches) > 0)
  {
    patches = new K_PATCH* [num_patches];
    
    for (i = 0; i < num_patches; i++)
    {
      patches[i] = s.patches[i];
      patches[i]->ref_count++;
    }
  }
  else  //  if (!num_patches)
    patches = 0;
}

K_SOLID& K_SOLID :: operator =(const K_SOLID& s)
{
  unsigned long i;
  
  if (this != &s)
  {
    if (num_patches > 0)
    {
      for (i = 0; i < num_patches; i++)
        if (!--patches[i]->ref_count)
          delete patches[i];
      
      delete [] patches;
    }
    
    if ((num_patches = s.num_patches) > 0)
    {
      patches = new K_PATCH* [num_patches];
      
      for (i = 0; i < num_patches; i++)
      {
        patches[i] = s.patches[i];
        patches[i]->ref_count++;
      }
    }
    else  //  if (!num_patches)
      patches = 0;
  }
  
  return *this;
}

K_SOLID :: ~K_SOLID()
{
  unsigned long i;
  
  if (num_patches > 0)
  {
    for (i = 0; i < num_patches; i++)
      if (!--patches[i]->ref_count)
        delete patches[i];
    
    delete [] patches;
  }
}

ostream& K_SOLID :: output(ostream& o) const
{
  unsigned long i;
  
  if (num_patches > 0)
  {
    for (i = 0; i < num_patches - 1; i++)
    {
      o << " patch[" << i << "] = " << *patches[i] << flush;
      o << " -------------------- " << endl << flush;
    }
    
    o << " patch[" << num_patches - 1 << "] = " << *patches[num_patches - 1] << flush;
  }
  else  //  if (!num_patches)
    o << " NULL";
  
  o << endl << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const K_SOLID& s)
{
  return s.output(o);
}

//  int K_SOLID :: classify_pt(const bigrational& x,
//                             const bigrational& y,
//                             const bigrational& z) const
//    returns IN  if (x, y, z) lies inside *this and
//            OUT if (x, y, z) lies outside *this.

int K_SOLID :: classify_pt(const bigrational& x,
                           const bigrational& y,
                           const bigrational& z) const
//  Ray 0: + x, 1: - x, 2: + y, 3: - y, 4: + z and 5: - z.
{
  assert(num_patches > 0);
  
  unsigned long  i, j;
  unsigned long  num_hits[6];
  unsigned long  total_num_hits;
  unsigned long* patches_hit[6];
  bigrational    low_s, high_s, low_t, high_t;
  K_BOX3D        b;
  unsigned long  ray_dir;
  K_RATPOLY      poly1, poly2;
  K_POINT2D**    int_pts;
  unsigned long  num_int_pts, num_int_pts_proto;
  int            separated;
  int            location;
  
  for (i = 0; i < 6; i++)
  {
    patches_hit[i] = new unsigned long [num_patches];  //  num_patches > 0
    num_hits[i]    = 0;
  }
  
  //  For each direction,
  //    count the number of patches POSSIBLY hit by a ray of the direction.
  
  for (i = 0; i < num_patches; i++)
  {
    patches[i]->set_range();
    low_s  = patches[i]->get_low_s();
    high_s = patches[i]->get_high_s();
    low_t  = patches[i]->get_low_t();
    high_t = patches[i]->get_high_t();
    b      = patches[i]->surf->get_range(low_s, high_s, low_t, high_t);
    
    if ((b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] < y)
        &&
        (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] > y)
        &&
        (b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] < z)
        &&
        (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] > z))
    {
      if (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] > x)
      {
        patches_hit[0][num_hits[0]] = i;
        num_hits[0]++;
      }
      
      if (b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] < x)
      {
        patches_hit[1][num_hits[1]] = i;
        num_hits[1]++;
      }
    }
    
    if ((b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] < x)
        &&
        (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] > x)
        &&
        (b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] < z)
        &&
        (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] > z))
    {
      if (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] > y)
      {
        patches_hit[2][num_hits[2]] = i;
        num_hits[2]++;
      }
      
      if (b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] < y)
      {
        patches_hit[3][num_hits[3]] = i;
        num_hits[3]++;
      }
    }
    
    if ((b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] < x)
        &&
        (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] > x)
        &&
        (b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] < y)
        &&
        (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] > y))
    {
      if (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] > z)
      {
        patches_hit[4][num_hits[4]] = i;
        num_hits[4]++;
      }
      
      if (b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] < z)
      {
        patches_hit[5][num_hits[5]] = i;
        num_hits[5]++;
      }
    }
  }
  
  //  Choose the ray with the smallest number of potential hits.
  
  ray_dir = 0;
  
  for (i = 1; i < 6; i++)
    if (num_hits[i] < num_hits[ray_dir])
      ray_dir = i;
  
  if (!num_hits[ray_dir])
  //  (x, y, z) must be outside of *this
    location = OUT;
  else  //  if (num_hits[i] >= num_hits[ray_dir] > 0)
  {
    total_num_hits = 0;
    
    //  For each patch potentially hit, see whether or not it is really hit.
    
    for (i = 0; i < num_hits[ray_dir]; i++)
    {
      //  Set polynomials to intersect.
      
      if (ray_dir == 0 || ray_dir == 1)
      {
        poly1 =
          *patches[patches_hit[ray_dir][i]]->surf->Y
          -
          y * *patches[patches_hit[ray_dir][i]]->surf->W;
        poly2 =
          *patches[patches_hit[ray_dir][i]]->surf->Z
          -
          z * *patches[patches_hit[ray_dir][i]]->surf->W;
      }
      else if (ray_dir == 2 || ray_dir == 3)
      {
        poly1 =
          *patches[patches_hit[ray_dir][i]]->surf->X
          -
          x * *patches[patches_hit[ray_dir][i]]->surf->W;
        poly2 =
          *patches[patches_hit[ray_dir][i]]->surf->Z
          -
          z * *patches[patches_hit[ray_dir][i]]->surf->W;
      }
      else  //  if (ray_dir == 4 || ray_dir == 5)
      {
        poly1 =
          *patches[patches_hit[ray_dir][i]]->surf->X
          -
          x * *patches[patches_hit[ray_dir][i]]->surf->W;
        poly2 =
          *patches[patches_hit[ray_dir][i]]->surf->Y
          -
          y * *patches[patches_hit[ray_dir][i]]->surf->W;
      }
      
      //  Find intersections of poly1 and poly2.
      
      patches[patches_hit[ray_dir][i]]->set_range();
      low_s  = patches[patches_hit[ray_dir][i]]->get_low_s();
      high_s = patches[patches_hit[ray_dir][i]]->get_high_s();
      low_t  = patches[patches_hit[ray_dir][i]]->get_low_t();
      high_t = patches[patches_hit[ray_dir][i]]->get_high_t();
      
      num_int_pts = num_int_pts_proto =
        get_pts(low_s, high_s, low_t, high_t, poly1, poly2, int_pts, 0, 0);
      
      //  Remove intersections outside of trim region.
      
      j = 0;
      
      while (j < num_int_pts)
      {
        if (!patches[patches_hit[ray_dir][i]]->contains(*int_pts[j]))
        {
          if (!--int_pts[j]->ref_count)
            delete int_pts[j];
          
          if (j < num_int_pts - 1)
//          {
            int_pts[j] = int_pts[num_int_pts - 1];
//            int_pts[j]->ref_count++;
//            
//            if (!--int_pts[num_int_pts - 1]->ref_count)
//              delete int_pts[num_int_pts - 1];
//          }
          else  //  if (j == num_int_pts - 1)
            j++;
          
          num_int_pts--;
        }
        else
          j++;
      }
      
      //  For each hit,
      //    see whether or not it is on the correct side of the ray.
      
      for (j = 0; j < num_int_pts; j++)
      {
        b         =
          patches[patches_hit[ray_dir][i]]->surf->
          get_range(int_pts[j]->get_low_s(), int_pts[j]->get_high_s(),
                    int_pts[j]->get_low_t(), int_pts[j]->get_high_t());
        separated = 1;
        
        if (((ray_dir == 0 || ray_dir == 1)
             &&
             ((b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] <= x)
              &&
              (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] >= x)))
            ||
            ((ray_dir == 2 || ray_dir == 3)
             &&
             ((b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] <= y)
              &&
              (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] >= y)))
            ||
            ((ray_dir == 4 || ray_dir == 5)
             &&
             ((b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] <= z)
              &&
              (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] >= z))))
          separated = 0;
        
        while (!separated)
        {
          int_pts[j]->shrink(shrink_step, shrink_step);
          b         =
          patches[patches_hit[ray_dir][i]]->surf->
            get_range(int_pts[j]->get_low_s(), int_pts[j]->get_high_s(),
                      int_pts[j]->get_low_t(), int_pts[j]->get_high_t());
          separated = 1;
          
          if (((ray_dir == 0 || ray_dir == 1)
               &&
               ((b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] <= x)
                &&
                (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] >= x)))
              ||
              ((ray_dir == 2 || ray_dir == 3)
               &&
               ((b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] <= y)
                &&
                (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] >= y)))
              ||
              ((ray_dir == 4 || ray_dir == 5)
               &&
               ((b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] <= z)
                &&
                (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] >= z))))
            separated = 0;
        }
        
        if ((ray_dir == 0
             &&
             (b.low_infty[0] > 0 || !b.low_infty[0] && b.low[0] > x))
            ||
            (ray_dir == 1
             &&
             (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] < x))
            ||
            (ray_dir == 2
             &&
             (b.low_infty[1] > 0 || !b.low_infty[1] && b.low[1] > y))
            ||
            (ray_dir == 3
             &&
             (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] < y))
            ||
            (ray_dir == 4
             &&
             (b.low_infty[2] > 0 || !b.low_infty[2] && b.low[2] > z))
            ||
            (ray_dir == 5
             &&
             (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] < z)))
          total_num_hits++;
      }
      
      for (j = 0; j < num_int_pts; j++)
        if (!--int_pts[j]->ref_count)
          delete int_pts[j];
      
      if (num_int_pts_proto > 0)
        delete [] int_pts;
    }
    
    if (total_num_hits % 2)
      location = IN;
    else
      location = OUT;
  }
  
  for (i = 0; i < 6; i++)
    delete [] patches_hit[i];
  
  return location;
}

//  int K_SOLID :: classify_pt(const bigrational_vector& pt) const
//    returns IN  if pt lies inside *this and
//            OUT if pt lies outside *this.

int K_SOLID :: classify_pt(const bigrational_vector& pt) const
//  Ray 0: + pt[0],    1: - pt[0],
//      2: + pt[1],    3: - pt[1],
//      4: + pt[2] and 5: - pt[2].
{
  assert(num_patches > 0);
  
  unsigned long  i, j;
  unsigned long  num_hits[6];
  unsigned long  total_num_hits;
  unsigned long* patches_hit[6];
  bigrational    low_s, high_s, low_t, high_t;
  K_BOX3D        b;
  unsigned long  ray_dir;
  K_RATPOLY      poly1, poly2;
  K_POINT2D**    int_pts;
  unsigned long  num_int_pts, num_int_pts_proto;
  int            separated;
  int            location;
  
  for (i = 0; i < 6; i++)
  {
    patches_hit[i] = new unsigned long [num_patches];  //  num_patches > 0
    num_hits[i]    = 0;
  }
  
  //  For each direction,
  //    count the number of patches POSSIBLY hit by a ray of the direction.
  
  for (i = 0; i < num_patches; i++)
  {
    patches[i]->set_range();
    low_s  = patches[i]->get_low_s();
    high_s = patches[i]->get_high_s();
    low_t  = patches[i]->get_low_t();
    high_t = patches[i]->get_high_t();
    b      = patches[i]->surf->get_range(low_s, high_s, low_t, high_t);
    
    if ((b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] < pt[1])
        &&
        (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] > pt[1])
        &&
        (b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] < pt[2])
        &&
        (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] > pt[2]))
    {
      if (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] > pt[0])
      {
        patches_hit[0][num_hits[0]] = i;
        num_hits[0]++;
      }
      
      if (b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] < pt[0])
      {
        patches_hit[1][num_hits[1]] = i;
        num_hits[1]++;
      }
    }
    
    if ((b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] < pt[0])
        &&
        (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] > pt[0])
        &&
        (b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] < pt[2])
        &&
        (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] > pt[2]))
    {
      if (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] > pt[1])
      {
        patches_hit[2][num_hits[2]] = i;
        num_hits[2]++;
      }
      
      if (b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] < pt[1])
      {
        patches_hit[3][num_hits[3]] = i;
        num_hits[3]++;
      }
    }
    
    if ((b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] < pt[0])
        &&
        (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] > pt[0])
        &&
        (b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] < pt[1])
        &&
        (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] > pt[1]))
    {
      if (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] > pt[2])
      {
        patches_hit[4][num_hits[4]] = i;
        num_hits[4]++;
      }
      
      if (b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] < pt[2])
      {
        patches_hit[5][num_hits[5]] = i;
        num_hits[5]++;
      }
    }
  }
  
  //  Choose the ray with the smallest number of potential hits.
  
  ray_dir = 0;
  
  for (i = 1; i < 6; i++)
    if (num_hits[i] < num_hits[ray_dir])
      ray_dir = i;
  
  if (!num_hits[ray_dir])
  //  pt must be outside of *this
    location = OUT;
  else  //  if (num_hits[i] >= num_hits[ray_dir] > 0)
  {
    total_num_hits = 0;
    
    //  For each patch potentially hit, see whether or not it is really hit.
    
    for (i = 0; i < num_hits[ray_dir]; i++)
    {
      //  Set polynomials to intersect.
      
      if (ray_dir == 0 || ray_dir == 1)
      {
        poly1 =
          *patches[patches_hit[ray_dir][i]]->surf->Y
          -
          pt[1] * *patches[patches_hit[ray_dir][i]]->surf->W;
        poly2 =
          *patches[patches_hit[ray_dir][i]]->surf->Z
          -
          pt[2] * *patches[patches_hit[ray_dir][i]]->surf->W;
      }
      else if (ray_dir == 2 || ray_dir == 3)
      {
        poly1 =
          *patches[patches_hit[ray_dir][i]]->surf->X
          -
          pt[0] * *patches[patches_hit[ray_dir][i]]->surf->W;
        poly2 =
          *patches[patches_hit[ray_dir][i]]->surf->Z
          -
          pt[2] * *patches[patches_hit[ray_dir][i]]->surf->W;
      }
      else  //  if (ray_dir == 4 || ray_dir == 5)
      {
        poly1 =
          *patches[patches_hit[ray_dir][i]]->surf->X
          -
          pt[0] * *patches[patches_hit[ray_dir][i]]->surf->W;
        poly2 =
          *patches[patches_hit[ray_dir][i]]->surf->Y
          -
          pt[1] * *patches[patches_hit[ray_dir][i]]->surf->W;
      }
      
      //  Find intersections of poly1 and poly2.
      
      patches[patches_hit[ray_dir][i]]->set_range();
      low_s  = patches[patches_hit[ray_dir][i]]->get_low_s();
      high_s = patches[patches_hit[ray_dir][i]]->get_high_s();
      low_t  = patches[patches_hit[ray_dir][i]]->get_low_t();
      high_t = patches[patches_hit[ray_dir][i]]->get_high_t();
      
      num_int_pts = num_int_pts_proto =
        get_pts(low_s, high_s, low_t, high_t, poly1, poly2, int_pts, 0, 0);
      
      //  Remove intersections outside of trim region.
      
      j = 0;
      
      while (j < num_int_pts)
      {
        if (!patches[patches_hit[ray_dir][i]]->contains(*int_pts[j]))
        {
          if (!--int_pts[j]->ref_count)
            delete int_pts[j];
          
          if (j < num_int_pts - 1)
//          {
            int_pts[j] = int_pts[num_int_pts - 1];
//            int_pts[j]->ref_count++;
//            
//            if (!--int_pts[num_int_pts - 1]->ref_count)
//              delete int_pts[num_int_pts - 1];
//          }
          else  //  if (j == num_int_pts - 1)
            j++;
          
          num_int_pts--;
        }
        else
          j++;
      }
      
      //  For each hit,
      //    see whether or not it is on the correct side of the ray.
      
      for (j = 0; j < num_int_pts; j++)
      {
        b         =
          patches[patches_hit[ray_dir][i]]->surf->
          get_range(int_pts[j]->get_low_s(), int_pts[j]->get_high_s(),
                    int_pts[j]->get_low_t(), int_pts[j]->get_high_t());
        separated = 1;
        
        if (((ray_dir == 0 || ray_dir == 1)
             &&
             ((b.low_infty[0] < 0 || !b.low_infty[0] && b.low[0] <= pt[0])
              &&
              (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] >= pt[0])))
            ||
            ((ray_dir == 2 || ray_dir == 3)
             &&
             ((b.low_infty[1] < 0 || !b.low_infty[1] && b.low[1] <= pt[1])
              &&
              (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] >= pt[1])))
            ||
            ((ray_dir == 4 || ray_dir == 5)
             &&
             ((b.low_infty[2] < 0 || !b.low_infty[2] && b.low[2] <= pt[2])
              &&
              (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] >= pt[2]))))
          separated = 0;
        
        while (!separated)
        {
          int_pts[j]->shrink(shrink_step, shrink_step);
          b         =
          patches[patches_hit[ray_dir][i]]->surf->
            get_range(int_pts[j]->get_low_s(), int_pts[j]->get_high_s(),
                      int_pts[j]->get_low_t(), int_pts[j]->get_high_t());
          separated = 1;
          
          if (((ray_dir == 0 || ray_dir == 1)
               &&
               ((b.low_infty[0] < 0
                 ||
                 !b.low_infty[0] && b.low[0] <= pt[0])
                &&
                (b.high_infty[0] > 0
                 ||
                 !b.high_infty[0] && b.high[0] >= pt[0])))
              ||
              ((ray_dir == 2 || ray_dir == 3)
               &&
               ((b.low_infty[1] < 0
                 ||
                 !b.low_infty[1] && b.low[1] <= pt[1])
                &&
                (b.high_infty[1] > 0
                 ||
                 !b.high_infty[1] && b.high[1] >= pt[1])))
              ||
              ((ray_dir == 4 || ray_dir == 5)
               &&
               ((b.low_infty[2] < 0
                 ||
                 !b.low_infty[2] && b.low[2] <= pt[2])
                &&
                (b.high_infty[2] > 0
                 ||
                 !b.high_infty[2] && b.high[2] >= pt[2]))))
            separated = 0;
        }
        
        if ((ray_dir == 0
             &&
             (b.low_infty[0] > 0 || !b.low_infty[0] && b.low[0] > pt[0]))
            ||
            (ray_dir == 1
             &&
             (b.high_infty[0] > 0 || !b.high_infty[0] && b.high[0] < pt[0]))
            ||
            (ray_dir == 2
             &&
             (b.low_infty[1] > 0 || !b.low_infty[1] && b.low[1] > pt[1]))
            ||
            (ray_dir == 3
             &&
             (b.high_infty[1] > 0 || !b.high_infty[1] && b.high[1] < pt[1]))
            ||
            (ray_dir == 4
             &&
             (b.low_infty[2] > 0 || !b.low_infty[2] && b.low[2] > pt[2]))
            ||
            (ray_dir == 5
             &&
             (b.high_infty[2] > 0 || !b.high_infty[2] && b.high[2] < pt[2])))
          total_num_hits++;
      }
      
      for (j = 0; j < num_int_pts; j++)
        if (!--int_pts[j]->ref_count)
          delete int_pts[j];
      
      if (num_int_pts_proto > 0)
        delete [] int_pts;
    }
    
    if (total_num_hits % 2)
      location = IN;
    else
      location = OUT;
  }
  
  for (i = 0; i < 6; i++)
    delete [] patches_hit[i];
  
  return location;
}

//K_SOLID gen_new_solid(K_PARTITION** partitions1,
//                      const unsigned long num_partitions1,
//                      K_PARTITION** partitions2,
//                      const unsigned long num_partitions2)
K_SOLID gen_new_solid(K_PARTITION** partitions1,
                      const unsigned long num_partitions1,
                      K_PARTITION** partitions2,
                      const unsigned long num_partitions2,
                      const char op)
{
  unsigned long i, j, k, l;
  K_PARTITION** new_partitions;
  K_PATCH**     new_patches;
  unsigned long num_new_patches;
  int           found;
  K_SOLID       r;
  
  num_new_patches = num_partitions1 + num_partitions2;
  new_partitions  = new K_PARTITION* [num_new_patches];
  new_patches     = new K_PATCH* [num_new_patches];
  
  for (i = 0; i < num_partitions1; i++)
  {
    new_partitions[i] = new K_PARTITION(*partitions1[i]);
    new_partitions[i]->ref_count++;
  }
  
  for (i = 0; i < num_partitions2; i++)
  {
    new_partitions[i + num_partitions1] = new K_PARTITION(*partitions2[i]);
    new_partitions[i + num_partitions1]->ref_count++;
  }
  
  for (i = 0; i < num_partitions1; i++)
    for (j = 0; j < partitions1[i]->num_trim_curves; j++)
      if (partitions1[i]->adj_partitions[j])
      {
        found = 0;
        k     = 0;
        
        while (!found && k < num_partitions2)
        {
          for (l = 0; !found && l < partitions2[k]->num_trim_curves; l++)
            if (partitions2[k]->trim_curves[l] ==
                partitions1[i]->trim_curves[j]->curve_in_other_dom)
              found = 1;
          
          if (!found)
            k++;
        }
        
//        cerr << " ksolid: gen_new_solid: 1-2: k = " << k << ", num_partitions2 = " << num_partitions2 << endl << flush;
        assert(k < num_partitions2);
        
        if (new_partitions[i]->adj_surfs[j]
            &&
            !--new_partitions[i]->adj_surfs[j]->ref_count)
          delete new_partitions[i]->adj_surfs[j];
        
        new_partitions[i]->adj_surfs[j] = partitions2[k]->from->surf;
        new_partitions[i]->adj_surfs[j]->ref_count++;
        
        new_partitions[i]->adj_curves[j] = k + num_partitions1;
      }
      else  //  if (!partitions1[i]->adj_partitions[j])
      {
        found = 0;
        k     = 0;
        
        while (!found && k < num_partitions1)
        {
          if (partitions1[k]->from ==
              partitions1[i]->from->adj_patches[partitions1[i]->adj_curves[j]])
            for (l = 0; !found && l < partitions1[k]->num_trim_curves; l++)
              if (partitions1[k]->trim_curves[l] ==
                  partitions1[i]->from->
                  trim_curves[partitions1[i]->adj_curves[j]]->
                  curve_in_other_dom)
                found = 1;
          
          if (!found)
            k++;
        }
        
//        cerr << " ksolid: gen_new_solid: 1-1: k = " << k << ", num_partitions1 = " << num_partitions1 << endl << flush;
        assert(k < num_partitions1);
        
        if (new_partitions[i]->adj_surfs[j]
            &&
            !--new_partitions[i]->adj_surfs[j]->ref_count)
          delete new_partitions[i]->adj_surfs[j];
        
        new_partitions[i]->adj_surfs[j] =
          partitions1[i]->from->adj_surfs[new_partitions[i]->adj_curves[j]];
        new_partitions[i]->adj_surfs[j]->ref_count++;
        
        new_partitions[i]->adj_curves[j] = k;
      }
  
  for (i = 0; i < num_partitions2; i++)
  {
    if (op == 'D')
      new_partitions[i + num_partitions1]->is_head =
        !new_partitions[i + num_partitions1]->is_head;
    
    for (j = 0; j < partitions2[i]->num_trim_curves; j++)
      if (partitions2[i]->adj_partitions[j])
      {
        found = 0;
        k     = 0;
        
        while (!found && k < num_partitions1)
        {
          for (l = 0; !found && l < partitions1[k]->num_trim_curves; l++)
            if (partitions1[k]->trim_curves[l] ==
                partitions2[i]->trim_curves[j]->curve_in_other_dom)
              found = 1;
          
          if (!found)
            k++;
        }
        
//        cerr << " ksolid: gen_new_solid: 2-1: k = " << k << ", num_partitions1 = " << num_partitions1 << endl << flush;
        assert(k < num_partitions1);
        
        if (new_partitions[i + num_partitions1]->adj_surfs[j]
            &&
            !--new_partitions[i + num_partitions1]->adj_surfs[j]->ref_count)
          delete new_partitions[i + num_partitions1]->adj_surfs[j];
        
        new_partitions[i + num_partitions1]->adj_surfs[j] =
          partitions1[k]->from->surf;
        new_partitions[i + num_partitions1]->adj_surfs[j]->ref_count++;
        
        new_partitions[i + num_partitions1]->adj_curves[j] = k;
      }
      else  //  if (!partitions2[i]->adj_partitions[j])
      {
        found = 0;
        k     = 0;
        
        while (!found && k < num_partitions2)
        {
          if (partitions2[k]->from ==
              partitions2[i]->from->adj_patches[partitions2[i]->adj_curves[j]])
            for (l = 0; !found && l < partitions2[k]->num_trim_curves; l++)
              if (partitions2[k]->trim_curves[l] ==
                  partitions2[i]->from->
                  trim_curves[partitions2[i]->adj_curves[j]]->
                  curve_in_other_dom)
                found = 1;
          
          if (!found)
            k++;
        }
        
//        cerr << " ksolid: gen_new_solid: 2-2: k = " << k << ", num_partitions2 = " << num_partitions2 << endl << flush;
        assert(k < num_partitions2);
        
        if (new_partitions[i + num_partitions1]->adj_surfs[j]
            &&
            !--new_partitions[i + num_partitions1]->adj_surfs[j]->ref_count)
          delete new_partitions[i + num_partitions1]->adj_surfs[j];
        
        new_partitions[i + num_partitions1]->adj_surfs[j] =
          partitions2[i]->from->
          adj_surfs[new_partitions[i + num_partitions1]->adj_curves[j]];
        new_partitions[i + num_partitions1]->adj_surfs[j]->ref_count++;
        
        new_partitions[i + num_partitions1]->adj_curves[j]
          = k + num_partitions1;
      }
  }
  
  for (i = 0; i < num_new_patches; i++)
  {
    new_patches[i] = new K_PATCH(*new_partitions[i]);
    new_patches[i]->ref_count++;
  }
  
  for (i = 0; i < num_new_patches; i++)
    for (j = 0; j < new_partitions[i]->num_trim_curves; j++)
    {
      new_patches[i]->adj_surfs[j]   = new_partitions[i]->adj_surfs[j];
      new_patches[i]->adj_surfs[j]->ref_count++;
      new_patches[i]->adj_patches[j] =
        new_patches[new_partitions[i]->adj_curves[j]];
      new_patches[i]->adj_patches[j]->ref_count++;
    }
  
  r = K_SOLID(new_patches, num_new_patches);
  
  for (i = 0; i < num_new_patches; i++)
    if (!--new_partitions[i]->ref_count)
      delete new_partitions[i];
  
  delete [] new_partitions;
  
  for (i = 0; i < num_new_patches; i++)
    if (!--new_patches[i]->ref_count)
      delete new_patches[i];
  
  delete [] new_patches;
  
  return r;
}

K_SOLID K_SOLID :: boolean(K_SOLID& s, const char op)
{
#ifdef _EXPERIMENT
  num_ksolid_boolean++;
#endif
  
  unsigned long      i, j, k;
  K_SOLID            r;
  K_BOX3D*           t_boxes;
  K_BOX3D*           s_boxes;
  unsigned long      num_closed_curves;
  K_PATCH**          new_patches;
  unsigned long      num_new_patches;
  K_PATCH**          patches_proto;
  unsigned long      num_patches_proto;
  K_PARTITION**      t_partitions;
  unsigned long      num_t_partitions;
  K_PARTITION**      s_partitions;
  unsigned long      num_s_partitions;
  K_PARTITION**      partitions_proto;
  unsigned long      num_partitions_proto;
  K_GRAPH*           t_partition_graph;
  long*              t_components;
  K_GRAPH*           t_component_graph;
  K_GRAPH*           s_partition_graph;
  long*              s_components;
  K_GRAPH*           s_component_graph;
  K_POINT2D          p;
//  bigrational        x, y, z;
  bigrational_vector v(3);
  int                t_classified, s_classified;
  int*               t_colors;
  int*               s_colors;
  int                t_new_color, s_new_color;
  K_PARTITION**      t_new_partitions;
  unsigned long      num_t_new_partitions;
  K_GRAPH*           t_SG;
  K_PARTITION**      s_new_partitions;
  unsigned long      num_s_new_partitions;
  K_GRAPH*           s_SG;
  
  if (num_patches > 0 && s.num_patches > 0)
  {
    t_boxes = new K_BOX3D [num_patches];
    s_boxes = new K_BOX3D [s.num_patches];
    
    for (i = 0; i < num_patches; i++)
      t_boxes[i] =
        patches[i]->surf->get_range(patches[i]->low_s, patches[i]->high_s,
                                    patches[i]->low_t, patches[i]->high_t);
    
    for (i = 0; i < s.num_patches; i++)
      s_boxes[i] =
        s.patches[i]->surf->get_range(s.patches[i]->low_s,
                                      s.patches[i]->high_s,
                                      s.patches[i]->low_t,
                                      s.patches[i]->high_t);
    
//    for (i = 0; i < num_patches; i++)
//    {
//      cerr << " ksolid: boolean: patches[" << i << "] = (" << patches[i]->low_s << ", " << patches[i]->high_s << ") x (" << patches[i]->low_t << ", " << patches[i]->high_t << ")" << endl << flush;
//      cerr << " ksolid: boolean: t_boxes[" << i << "] = " << t_boxes[i] << endl << flush;
//    }
    
//    for (i = 0; i < s.num_patches; i++)
//    {
//      cerr << " ksolid: boolean: s.patches[" << i << "] = (" << s.patches[i]->low_s << ", " << s.patches[i]->high_s << ") x (" << s.patches[i]->low_t << ", " << s.patches[i]->high_t << ")" << endl << flush;
//      cerr << " ksolid: boolean: s_boxes[" << i << "] = " << s_boxes[i] << endl << flush;
//    }
    
    for (i = 0; i < num_patches; i++)
      for (j = 0; j < s.num_patches; j++)
        if (t_boxes[i].overlap(s_boxes[j]))
        {
          cerr << " Intersecting patch " << i << " and " << j << "." << endl << flush;
          
          patches[i]->intersect(*s.patches[j]);
          
          cerr << "   Patch " << i << " on solid 1 has " << patches[i]->num_int_curves << " intersection curves, and " << endl << flush;
          cerr << "   patch " << j << " on solid 2 has " << s.patches[j]->num_int_curves << " intersection curves. " << endl << endl << flush;
        }
        else
          cerr << " Intersecting patch " << i << " and " << j << "; they do not overalp." << endl << endl << flush;
    
    delete [] t_boxes;
    
    for (i = 0; i < num_patches; i++)
      cerr << "   Patch " << i << " on solid 1 has " << patches[i]->num_int_curves << " intersection curves, and " << endl << flush;
    cerr << endl << flush;
    for (j = 0; j < s.num_patches; j++)
      cerr << "   patch " << j << " on solid 2 has " << s.patches[j]->num_int_curves << " intersection curves. " << endl << flush;
    cerr << endl << flush;
  }
  
  //  Split loops
  
  i = 0;
  
  while (i < num_patches)
  {
    num_closed_curves = patches[i]->merge_curves();
    cerr << " Patch " << i << " on solid 1 has " << num_closed_curves << " loops. " << endl << flush;
    
    if (num_closed_curves > 0)
    {
      num_new_patches = patches[i]->split_loops(new_patches);
//      cerr << " ksolid: boolean: num_new_patches = " << num_new_patches << endl << flush;
      num_patches_proto = num_patches + num_new_patches - 1;
      patches_proto     = new K_PATCH* [num_patches_proto];
      
      for (j = 0; j < i; j++)
      {
        patches_proto[j] = patches[j];
        patches_proto[j]->ref_count++;
      }
      
      for (k = 0; k < num_new_patches; k++)
      {
        patches_proto[i + k] = new_patches[k];
        patches_proto[i + k]->ref_count++;
      }
      
      for (j = i + 1; j < num_patches; j++)
      {
        patches_proto[j + num_new_patches - 1] = patches[j];
        patches_proto[j + num_new_patches - 1]->ref_count++;
      }
      
      for (j = 0; j < num_patches; j++)
        if (!--patches[j]->ref_count)
          delete patches[j];
      
      delete [] patches;  //  num_pathces > 0
      
      patches      = patches_proto;
      num_patches  = num_patches_proto;
      i           += num_new_patches - 1;
//      cerr << " ksolid: boolean: num_patches = " << num_patches << endl << flush;
      
      for (k = 0; k < num_new_patches; k++)
        if (!--new_patches[k]->ref_count)
          delete new_patches[k];
      
      delete [] new_patches;  //  num_new_patches >= 1
    }
    
    i++;
  }
  
  i = 0;
  
  while (i < s.num_patches)
  {
    num_closed_curves = s.patches[i]->merge_curves();
    cerr << " Patch " << i << " on solid 2 has " << num_closed_curves << " loops. " << endl << flush;
    
    if (num_closed_curves > 0)
    {
      num_new_patches = s.patches[i]->split_loops(new_patches);
//      cerr << " ksolid: boolean: num_new_patches = " << num_new_patches << endl << flush;
      num_patches_proto = s.num_patches + num_new_patches - 1;
      patches_proto     = new K_PATCH* [num_patches_proto];
      
      for (j = 0; j < i; j++)
      {
        patches_proto[j] = s.patches[j];
        patches_proto[j]->ref_count++;
      }
      
      for (k = 0; k < num_new_patches; k++)
      {
        patches_proto[i + k] = new_patches[k];
        patches_proto[i + k]->ref_count++;
      }
      
      for (j = i + 1; j < s.num_patches; j++)
      {
        patches_proto[j + num_new_patches - 1] = s.patches[j];
        patches_proto[j + num_new_patches - 1]->ref_count++;
      }
      
      for (j = 0; j < s.num_patches; j++)
        if (!--s.patches[j]->ref_count)
          delete s.patches[j];
      
      delete [] s.patches;  //  s.num_pathces > 0
      
      s.patches      = patches_proto;
      s.num_patches  = num_patches_proto;
      i             += num_new_patches - 1;
//      cerr << " ksolid: boolean: s.num_patches = " << s.num_patches << endl << flush;
      
      for (k = 0; k < num_new_patches; k++)
        if (!--new_patches[k]->ref_count)
          delete new_patches[k];
      
      delete [] new_patches;  //  num_new_patches >= 1
    }
    
    i++;
  }
  
  //  Generate partitions.
  
  cerr << " Generating partitions for patches on solid 1 " << endl << flush;
  
  num_t_partitions = 0;
  
  for (i = 0; i < num_patches; i++)
    num_t_partitions += patches[i]->num_merged + 1;
  
//  cerr << " ksolid: boolean: num_t_partitions = " << num_t_partitions << endl << flush;
  
  if (num_t_partitions > 0)
    t_partitions = new K_PARTITION* [num_t_partitions];
  else  //  if (num_t_partitions == 0)
    t_partitions = 0;
  
  num_t_partitions = 0;
  
  for (i = 0; i < num_patches; i++)
  {
    num_partitions_proto = gen_partitions(patches[i], partitions_proto);
    
    cerr << " Patch " << i << " forms " << num_partitions_proto << " partitions on solid 1 " << endl << flush;
    
    for (j = 0; j < num_partitions_proto; j++)
    {
      t_partitions[num_t_partitions] = partitions_proto[j];
      t_partitions[num_t_partitions]->ref_count++;
      num_t_partitions++;
    }
    
    for (j = 0; j < num_partitions_proto; j++)
      if (!--partitions_proto[j]->ref_count)
        delete partitions_proto[j];
    
    delete [] partitions_proto;
  }
  
  cerr << " Generating partitions for patches on solid 2 " << endl << flush;
  
  num_s_partitions = 0;
  
  for (i = 0; i < s.num_patches; i++)
    num_s_partitions += s.patches[i]->num_merged + 1;
  
//  cerr << " ksolid: boolean: num_s_partitions = " << num_s_partitions << endl << flush;
  
  if (num_s_partitions > 0)
    s_partitions = new K_PARTITION* [num_s_partitions];
  else  //  if (num_s_partitions == 0)
    s_partitions = 0;
  
  num_s_partitions = 0;
  
  for (i = 0; i < s.num_patches; i++)
  {
    num_partitions_proto = gen_partitions(s.patches[i], partitions_proto);
    
    cerr << " Patch " << i << " forms " << num_partitions_proto << " partitions on solid 2 " << endl << flush;
    
    for (j = 0; j < num_partitions_proto; j++)
    {
      s_partitions[num_s_partitions] = partitions_proto[j];
      s_partitions[num_s_partitions]->ref_count++;
      num_s_partitions++;
    }
    
    for (j = 0; j < num_partitions_proto; j++)
      if (!--partitions_proto[j]->ref_count)
        delete partitions_proto[j];
    
    delete [] partitions_proto;
  }
  
  //  Generate graphs.
  
  if (num_t_partitions > 0 && num_s_partitions > 0)
  {
    gen_adjacency(t_partitions, num_t_partitions,
                  t_partition_graph, t_components, t_component_graph);
    gen_adjacency(s_partitions, num_s_partitions,
                  s_partition_graph, s_components, s_component_graph);
    
    cerr << " Shooting 3D rays. " << endl << flush;
    
    p = t_partitions[0]->get_pt_in();
//    patches[0]->surf->param_to_coord(p.get_low_s(), p.get_low_t(), x, y, z);
//    t_classified = s.classify_pt(x, y, z);
    patches[0]->surf->param_to_coord(p.get_low(), v);
    t_classified = s.classify_pt(v);
//    cerr << " ksolid: boolean: p = " << p << ", ( x, y, z ) = ( " << x << ", " << y << ", " << z << " ), classified = " << t_classified << endl << flush;
    
    p = s_partitions[0]->get_pt_in();
//    s.patches[0]->surf->param_to_coord(p.get_low_s(), p.get_low_t(), x, y, z);
//    s_classified = classify_pt(x, y, z);
    s.patches[0]->surf->param_to_coord(p.get_low(), v);
    s_classified = classify_pt(v);
//    cerr << " ksolid: boolean: p = " << p << ", ( x, y, z ) = ( " << x << ", " << y << ", " << z << " ), classified = " << s_classified << endl << flush;
    
    t_component_graph->propagate_color(t_components[0], t_colors, t_classified);
    s_component_graph->propagate_color(s_components[0], s_colors, s_classified);
    
    cerr << " Selecting partitions of the resulting solid." << endl << flush;
    assert(op == 'U' || op == 'I' || op == 'D');
//    cerr << " ksolid: boolean: op = " << op << endl << flush;
    
    if (op == 'U')
    {
      t_new_color = OUT;
      s_new_color = OUT;
    }
    else if (op == 'I')
    {
      t_new_color = IN;
      s_new_color = IN;
    }
    else  //  if (op == 'D')
    {
      t_new_color = OUT;
      s_new_color = IN;
    }
    
    num_t_new_partitions = select_relevant_partitions(t_partitions,
                                                      num_t_partitions,
                                                      t_partition_graph,
                                                      t_components,
                                                      t_colors,
                                                      t_new_color,
                                                      t_new_partitions,
                                                      t_SG);
    num_s_new_partitions = select_relevant_partitions(s_partitions,
                                                      num_s_partitions,
                                                      s_partition_graph,
                                                      s_components,
                                                      s_colors,
                                                      s_new_color,
                                                      s_new_partitions,
                                                      s_SG);
    
    cerr << " " << num_t_new_partitions << " partitions are selected from solid 1." << endl << flush;
    cerr << " " << num_s_new_partitions << " partitions are selected from solid 2." << endl << flush;
  }
  
  //  Build a resulting solid.
  
  cerr << " Build the resulting solid." << endl << flush;
  
  if (num_t_partitions > 0 && num_s_partitions > 0)
//    r = gen_new_solid(t_new_partitions, num_t_new_partitions,
//                      s_new_partitions, num_s_new_partitions);
    r = gen_new_solid(t_new_partitions, num_t_new_partitions,
                      s_new_partitions, num_s_new_partitions,
                      op);
  else if (num_t_partitions > 0)
    r = *this;
  else if (num_s_partitions > 0)
    r = s;
  
  cerr << " The resulting solid has " << r.num_patches << " patches." << endl << flush;
  
  if (num_t_partitions > 0)
  {
    for (i = 0; i < num_t_partitions; i++)
      if (!--t_partitions[i]->ref_count)
        delete t_partitions[i];
    
    delete [] t_partitions;
  }
  
  if (num_s_partitions > 0)
  {
    for (i = 0; i < num_s_partitions; i++)
      if (!--s_partitions[i]->ref_count)
        delete s_partitions[i];
    
    delete [] s_partitions;
  }
  
  if (num_t_partitions > 0 && num_s_partitions > 0)
  {
    delete t_partition_graph;
    delete t_component_graph;
    delete s_partition_graph;
    delete s_component_graph;
    
    delete [] t_colors;
    delete [] s_colors;
    
    for (i = 0; i < num_t_new_partitions; i++)
      if (!--t_new_partitions[i]->ref_count)
        delete t_new_partitions[i];
    
    delete [] t_new_partitions;
    
    for (i = 0; i < num_s_new_partitions; i++)
      if (!--s_new_partitions[i]->ref_count)
        delete s_new_partitions[i];
    
    delete [] s_new_partitions;
    
    delete t_SG;
    delete s_SG;
  }
  
  return r;
}

int K_SOLID :: Bezier_output(ostream& out_fs) const
{
  unsigned long i;
  
  out_fs << num_patches << endl << flush;
  
  for (i = 0; i < num_patches; i++)
  {
    out_fs << endl << flush;
    cerr << " K_SOLID :: Bezier_output: patch " << i << endl << flush;
    patches[i]->Bezier_output(out_fs, 172, 172, 172);  //  color: grey
//    patches[i]->Bezier_output(out_fs, 0, 127, 127);  //  color: green
  }
  
  return 0;
}

K_SOLID read_solid(istream& in_fs, const bigrational& perturb_factor)
{
  unsigned long type;
  K_SOLID       s;
  
  in_fs >> type;
  cerr << " ksolid: read_solid: type = " << type << endl << flush;
  
  if (sgn(perturb_factor))
    cerr << " ksolid: read_solid: perturb_factor = " << perturb_factor << endl << flush;
  
  if (type == 1)
    s = read_box(in_fs, perturb_factor);
  else if (type == 11)
    s = read_BRLCAD_box(in_fs, perturb_factor);
  else if (type == 2)
    s = read_cyl(in_fs, perturb_factor);
  else if (type == 12)
    s = read_BRLCAD_cyl(in_fs, perturb_factor);
  else if (type == 3)
    s = read_ell(in_fs, perturb_factor);
  else if (type == 13)
    s = read_BRLCAD_ell(in_fs, perturb_factor);
  else if (type == 4)
    s = read_tor(in_fs, perturb_factor);
  else if (type == 14)
    s = read_BRLCAD_tor(in_fs, perturb_factor);
  else
  {
    cerr << " ksolid: read: wrong type: " << type << endl << flush;
    abort();
  }
  
  return s;
}

K_SOLID K_SOLID :: merge(const K_SOLID& s) const
{
  unsigned long i;
  unsigned long n;
  K_PATCH**     p;
  
  if ((n = num_patches + s.num_patches) > 0)
  {
    p = new K_PATCH* [n];
    
    for (i = 0; i < num_patches; i++)
    {
      p[i] = patches[i];
      p[i]->ref_count++;
    }
    
    for (i = 0; i < s.num_patches; i++)
    {
      p[i + num_patches] = s.patches[i];
      p[i + num_patches]->ref_count++;
    }
  }
  else  //  if (!n)
    p = 0;
  
  return K_SOLID(p, n);
}

#define DIM_BRLCAD_MATRIX 4

bigrational_matrix read_BRLCAD_matrix(istream& in_fs)
{
  unsigned long      i, j;
  float              f;
  bigrational_matrix T(DIM_BRLCAD_MATRIX, DIM_BRLCAD_MATRIX);
  
  cerr << " ksolid: read_BRLCAD_matrix: T = " << endl << flush;
  
  in_fs >> f;
  
  if (!in_fs.eof())
  {
    for (i = 0; i < DIM_BRLCAD_MATRIX - 1; i++)
    {
      for (j = 0; j < DIM_BRLCAD_MATRIX; j++)
      {
        T(i, j) = as_bigrational(f);
        cerr << f << " => " << T(i, j) << ", " << flush;
        in_fs >> f;
      }
      
      cerr << endl << flush;
    }
    
    for (j = 0; j < DIM_BRLCAD_MATRIX - 1; j++)
    {
      T(DIM_BRLCAD_MATRIX - 1, j) = as_bigrational(f);
      cerr << f << " => " << T(DIM_BRLCAD_MATRIX - 1, j) << ", " << flush;
      in_fs >> f;
    }
    
    T(DIM_BRLCAD_MATRIX - 1, DIM_BRLCAD_MATRIX - 1) = as_bigrational(f);
    cerr << f << " => " << T(DIM_BRLCAD_MATRIX - 1, DIM_BRLCAD_MATRIX - 1) << ", " << endl << flush;
  }
  else  //  if (in_fs.eof())
    for (i = 0; i < DIM_BRLCAD_MATRIX; i++)
    {
      for (j = 0; j < DIM_BRLCAD_MATRIX; j++)
      {
        T(i, j) = i == j ? 1 : 0;
        cerr << T(i, j) << ", " << flush;
      }
      
      cerr << endl << flush;
    }
  
  return T;
}

#define MAX_NUM_SURFS 256

K_SOLID K_SOLID :: transform(const bigrational_matrix& T) const
{
  K_SOLID t(*this);
  
  if (num_patches > 0 && !T.is_identity())
  {
    unsigned long      i, j, k;
    K_SURF**           surfs;
    unsigned long      num_surfs;
    int                transformed;
    K_SURF*            s;
    K_RATPOLY*         X_org;
    K_RATPOLY*         Y_org;
    K_RATPOLY*         Z_org;
    K_RATPOLY*         W_org;
    bigrational_matrix inv_T;
    K_RATPOLY*         Impl_org;
    long               d[3];
    
    surfs = new K_SURF* [MAX_NUM_SURFS];
    
//    for (i = 0; i < MAX_NUM_SURFS; i++)
//      surfs[i] = 0;
    
    num_surfs = 0;
    
    for (i = 0; i < t.num_patches; i++)
    {
      //  See whether or not t.patches[i]->surf has already been transformed.
      
      transformed = 0;
      s           = t.patches[i]->surf;
      j           = 0;
      
      while (!transformed && j < num_surfs)
        if (surfs[i] == s)
          transformed = 1;
        else
          j++;
      
      if (!transformed)
      {
        //  Transform parametric forms.
        
        if (s->mon_ok)
        {
          X_org = s->X;
//          X_org->ref_count++;
//          
//          if (!--s->X->ref_count)
//            delete s->X;
          
          Y_org = s->Y;
//          Y_org->ref_count++;
//          
//          if (!--s->Y->ref_count)
//            delete s->Y;
          
          Z_org = s->Z;
//          Z_org->ref_count++;
//          
//          if (!--s->Z->ref_count)
//            delete s->Z;
          
          W_org = s->W;
//          W_org->ref_count++;
//          
//          if (!--s->W->ref_count)
//            delete s->W;
          
          s->X = new K_RATPOLY(T(0, 0) * *X_org + T(0, 1) * *Y_org +
                               T(0, 2) * *Z_org + T(0, 3) * *W_org);
          s->X->ref_count++;
          
          s->Y = new K_RATPOLY(T(1, 0) * *X_org + T(1, 1) * *Y_org +
                               T(1, 2) * *Z_org + T(1, 3) * *W_org);
          s->Y->ref_count++;
          
          s->Z = new K_RATPOLY(T(2, 0) * *X_org + T(2, 1) * *Y_org +
                               T(2, 2) * *Z_org + T(2, 3) * *W_org);
          s->Z->ref_count++;
          
          s->W = new K_RATPOLY(T(3, 0) * *X_org + T(3, 1) * *Y_org +
                               T(3, 2) * *Z_org + T(3, 3) * *W_org);
          s->W->ref_count++;
          
          if (!--X_org->ref_count)
            delete X_org;
          
          if (!--Y_org->ref_count)
            delete Y_org;
          
          if (!--Z_org->ref_count)
            delete Z_org;
          
          if (!--W_org->ref_count)
            delete W_org;
        }
        
        //  Transform implicit forms.
        
        if (s->Impl_ok)
        {
          Impl_org = s->Impl;
//          Impl_org->ref_count++;
//          
//          if (!--s->Impl->ref_count)
//            delete s->Impl;
          
          inv_T      = T.inverse();
          s->Impl    = new K_RATPOLY(Impl_org->transform_Impl(inv_T));
          s->Impl->ref_count++;
          
          if (!--Impl_org->ref_count)
            delete Impl_org;
        }
        
        //  Add s to the list of transformed surfaces.
        
        assert(num_surfs < MAX_NUM_SURFS);
        surfs[num_surfs] = s;
        surfs[num_surfs]->ref_count++;
        num_surfs++;
      }
      
      //  Transform adj_surfs of trim_curves.
      
      for (j = 0; j < t.patches[i]->num_trim_curves; j++)
      {
        //  See whether or not
        //    t.patches[i]->adj_surfs[j] has already been transformed.
        
        transformed = 0;
        s           = t.patches[i]->adj_surfs[j];
        k           = 0;
        
        while (!transformed && k < num_surfs)
          if (surfs[k] == s)
            transformed = 1;
          else
            k++;
        
        if (!transformed)
        {
          //  Transform parametric forms.
          
          if (s->mon_ok)
          {
            X_org = s->X;
//            X_org->ref_count++;
//            
//            if (!--s->X->ref_count)
//              delete s->X;
            
            Y_org = s->Y;
//            Y_org->ref_count++;
//            
//            if (!--s->Y->ref_count)
//              delete s->Y;
            
            Z_org = s->Z;
//            Z_org->ref_count++;
//            
//            if (!--s->Z->ref_count)
//              delete s->Z;
            
            W_org = s->W;
//            W_org->ref_count++;
//            
//            if (!--s->W->ref_count)
//              delete s->W;
            
            s->X = new K_RATPOLY(T(0, 0) * *X_org + T(0, 1) * *Y_org +
                                 T(0, 2) * *Z_org + T(0, 3) * *W_org);
            s->X->ref_count++;
            
            s->Y = new K_RATPOLY(T(1, 0) * *X_org + T(1, 1) * *Y_org +
                                 T(1, 2) * *Z_org + T(1, 3) * *W_org);
            s->Y->ref_count++;
            
            s->Z = new K_RATPOLY(T(2, 0) * *X_org + T(2, 1) * *Y_org +
                                 T(2, 2) * *Z_org + T(2, 3) * *W_org);
            s->Z->ref_count++;
            
            s->W = new K_RATPOLY(T(3, 0) * *X_org + T(3, 1) * *Y_org +
                                 T(3, 2) * *Z_org + T(3, 3) * *W_org);
            s->W->ref_count++;
            
            if (!--X_org->ref_count)
              delete X_org;
            
            if (!--Y_org->ref_count)
              delete Y_org;
            
            if (!--Z_org->ref_count)
              delete Z_org;
            
            if (!--W_org->ref_count)
              delete W_org;
          }
          
          //  Transform implicit forms.
          
          if (s->Impl_ok)
          {
            Impl_org = s->Impl;
//            Impl_org->ref_count++;
//            
//            if (!--s->Impl->ref_count)
//              delete s->Impl;
            
//            inv_T      = T.inverse();
            s->Impl    = new K_RATPOLY(Impl_org->transform_Impl(inv_T));
            s->Impl->ref_count++;
            
            if (!--Impl_org->ref_count)
              delete Impl_org;
          }
          
          //  Add s to the list of transformed surfaces.
          
          assert(num_surfs < MAX_NUM_SURFS);
          surfs[num_surfs] = s;
          surfs[num_surfs]->ref_count++;
          num_surfs++;
        }
      }
    }
    
    for (i = 0; i < num_surfs; i++)
      if (!--surfs[i]->ref_count)
        delete surfs[i];
    
    delete [] surfs;
  }
  
  return t;
}

#define MAX_LEN_FILE_NAME 128
#define MAX_LEN_PATH_NAME 1024

K_SOLID read_CSG(const char* solid_info_dir_proto, const char* out_file_proto,
                 const bigrational& perturb_factor)
{
  ifstream           in_fs;
  ofstream           out_fs;
  char               solid_info[MAX_LEN_PATH_NAME];
  char               solid_info_subdir1[MAX_LEN_FILE_NAME];
  char               solid_info_subdir2[MAX_LEN_FILE_NAME];
  char               solid_info_subdir[MAX_LEN_FILE_NAME];
  char               solid_info_dir1[MAX_LEN_PATH_NAME];
  char               solid_info_dir2[MAX_LEN_PATH_NAME];
  char               solid_info_dir[MAX_LEN_PATH_NAME];
  unsigned long      type;
  unsigned long      op_num;
  char               op;
  bigrational_matrix T;
  char               out_file[MAX_LEN_PATH_NAME];
  K_SOLID            s1_proto, s1, s2_proto, s2, s3;
  bigrational        perturb_factor1, perturb_factor2;
  
  strcpy(solid_info, solid_info_dir_proto);
  strcat(solid_info, "/solidinfo");
  in_fs.open(solid_info);
  in_fs >> type;
  cerr << " type = " << type << endl << flush;
  cerr << " solid_info = " << solid_info << endl << flush;
  
  if (sgn(perturb_factor))
    cerr << " ksolid: read_CSG: perturb_factor = " << perturb_factor << endl << flush;
  
  if (type == 1)
    s3 = read_box(in_fs, perturb_factor);
  else if (type == 11)
    s3 = read_BRLCAD_box(in_fs, perturb_factor);
  else if (type == 2)
    s3 = read_cyl(in_fs, perturb_factor);
  else if (type == 12)
    s3 = read_BRLCAD_cyl(in_fs, perturb_factor);
  else if (type == 3)
    s3 = read_ell(in_fs, perturb_factor);
  else if (type == 13)
    s3 = read_BRLCAD_ell(in_fs, perturb_factor);
  else if (type == 4)
    s3 = read_tor(in_fs, perturb_factor);
  else if (type == 14)
    s3 = read_BRLCAD_tor(in_fs, perturb_factor);
  else if (type == 6)
  //  binary operation
  {
    in_fs >> solid_info_subdir1;
    strcpy(solid_info_dir1, solid_info_dir_proto);
    strcat(solid_info_dir1, "/");
    strcat(solid_info_dir1, solid_info_subdir1);
    cerr << " solid_info_dir1 = " << solid_info_dir1 << endl << flush;
    
    in_fs >> op_num;
    assert(op_num == 1 || op_num == 2 || op_num == 3);
    
    if (op_num == 1)
      op = 'U';
    else if (op_num == 2)
      op = 'I';
    else if (op_num == 3)
      op = 'D';
    
    cerr << " op_num = " << op_num << ", op = " << op << endl << flush;
    
    in_fs >> solid_info_subdir2;
    strcpy(solid_info_dir2, solid_info_dir_proto);
    strcat(solid_info_dir2, "/");
    strcat(solid_info_dir2, solid_info_subdir2);
    cerr << " solid_info_dir2 = " << solid_info_dir2 << endl << flush;
    
    if (op_num == 1)
    {
      perturb_factor1 = perturb_factor;
      perturb_factor2 = perturb_factor;
    }
    else if (op_num == 2)
    {
      perturb_factor1 = - perturb_factor;
      perturb_factor2 = - perturb_factor;
    }
    else  //  if (op_num == 3)
    {
      perturb_factor1 = perturb_factor;
      perturb_factor2 = - perturb_factor;
    }
    
    s1 = read_CSG(solid_info_dir1, out_file_proto, perturb_factor1);
//    s1 = read_CSG(solid_info_dir1, 0, perturb_factor1);
    
    s2_proto = read_CSG(solid_info_dir2, out_file_proto, perturb_factor2);
//    s2_proto = read_CSG(solid_info_dir2, 0, perturb_factor2);
    
    T  = read_BRLCAD_matrix(in_fs);
    s2 = s2_proto.transform(T);
    
    s3 = s1.boolean(s2, op);
  }
  else if (type == 7)
  //  join
  {
    while (in_fs >> solid_info_subdir)
    {
      strcpy(solid_info_dir, solid_info_dir_proto);
      strcat(solid_info_dir, "/");
      strcat(solid_info_dir, solid_info_subdir);
      cerr << " solid_info_dir = " << solid_info_dir << endl << flush;
      
      s1_proto = read_CSG(solid_info_dir, out_file_proto, perturb_factor);
//      s1_proto = read_CSG(solid_info_dir, 0, perturb_factor);
      
      T  = read_BRLCAD_matrix(in_fs);
      s1 = s1_proto.transform(T);
      
      s3 = s3.merge(s1);
    }
  }
  else
  {
    cerr << " ksolid: read_CSG: wrong type: " << type << endl << flush;
    abort();
  }
  
  in_fs.close();
  
//  if (out_file_proto)
//  {
//    strcpy(out_file, solid_info_dir_proto);
//    strcat(out_file, "/");
//    strcat(out_file, out_file_proto);
//    out_fs.open(out_file);
//    s3.Bezier_output(out_fs);
//    out_fs.close();
//  }
  
  return s3;
}

