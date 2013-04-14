//  file:    genell.cc
//  update:  03/22/03

//#include <config.h>

#include <genell.h>

#include <bigrational_vector.h>
#include <fpconversion.h>
#include <genbox.h>

//  int get_impl_ell(const bigrational_vector pts[4], K_RATPOLY*& impl)
//    computes the implicit formula for the ellipsoid
//      whose center is pts[0] and axes are pts[1], pts[2] and pts[3].

int get_impl_ell(const bigrational_vector pts[4], K_RATPOLY*& impl)
{
  long d[3];
  long p[3];
  
  d[0] = d[1] = d[2] = 1;
  
  K_RATPOLY poly1(3, d);
  K_RATPOLY poly2(3, d);
  K_RATPOLY poly3(3, d);
  
  p[0] = 1;
  p[1] = 0;
  p[2] = 0;
  poly1.get_coeff(p) = pts[1][0];
  poly2.get_coeff(p) = pts[2][0];
  poly3.get_coeff(p) = pts[3][0];
  
  p[0] = 0;
  p[1] = 1;
  p[2] = 0;
  poly1.get_coeff(p) = pts[1][1];
  poly2.get_coeff(p) = pts[2][1];
  poly3.get_coeff(p) = pts[3][1];
  
  p[0] = 0;
  p[1] = 0;
  p[2] = 1;
  poly1.get_coeff(p) = pts[1][2];
  poly2.get_coeff(p) = pts[2][2];
  poly3.get_coeff(p) = pts[3][2];
  
  p[0] = 0;
  p[1] = 0;
  p[2] = 0;
  poly1.get_coeff(p) = - in_prod(pts[1], pts[0]);
  poly2.get_coeff(p) = - in_prod(pts[2], pts[0]);
  poly3.get_coeff(p) = - in_prod(pts[3], pts[0]);
  
  bigrational p11, p22, p33;
  
  p11 = in_prod(pts[1], pts[1]);
  p22 = in_prod(pts[2], pts[2]);
  p33 = in_prod(pts[3], pts[3]);
  
  d[0] = d[1] = d[2] = 0;
  
  K_RATPOLY poly4(3, d);
  
  p[0] = p[1] = p[2] = 0;
  poly4.get_coeff(p) = - p11 * p22 * p33;
  
  impl  = new K_RATPOLY;
  *impl = (p22 * p33) / p11 * (poly1 * poly1) +
          (p33 * p11) / p22 * (poly2 * poly2) +
          (p11 * p22) / p33 * (poly3 * poly3) + poly4;
  
  return 0;
}

//  int get_impl_sphere(const bigrational_vector& center,
//                      const bigrational r_sq,
//                      K_RATPOLY*& impl)
//    computes the implicit formula for the sphere
//      centered at center and radius sqrt(r_sq).

int get_impl_sphere(const bigrational_vector& center, const bigrational r_sq,
                    K_RATPOLY*& impl)
{
  assert(center.get_dim() == 3);
  
  long d[3];
  long p[3];
  
  d[0] = d[1] = d[2] = 2;
  impl = new K_RATPOLY(3, d);
  
  p[0] = 2;
  p[1] = 0;
  p[2] = 0;
  impl->get_coeff(p) = 1;
  
  p[0] = 0;
  p[1] = 2;
  p[2] = 0;
  impl->get_coeff(p) = 1;
  
  p[0] = 0;
  p[1] = 0;
  p[2] = 2;
  impl->get_coeff(p) = 1;
  
  p[0] = 1;
  p[1] = 0;
  p[2] = 0;
  impl->get_coeff(p) = - 2 * center[0];
  
  p[0] = 0;
  p[1] = 1;
  p[2] = 0;
  impl->get_coeff(p) = - 2 * center[1];
  
  p[0] = 0;
  p[1] = 0;
  p[2] = 1;
  impl->get_coeff(p) = - 2 * center[2];
  
  p[0] = 0;
  p[1] = 0;
  p[2] = 0;
  impl->get_coeff(p) = in_prod(center, center);
  
  return 0;
}

//  int get_param_ell(const bigrational_vector& center,
//                    const bigrational_vector& a1,
//                    const bigrational_vector& a2,
//                    const bigrational_vector& a3,
//                    K_RATPOLY*& X,
//                    K_RATPOLY*& Y,
//                    K_RATPOLY*& Z,
//                    K_RATPOLY*& W)
//    computes the parametric formula (*X/*W, *Y/*W, *Z/*W) for the ellipsoid
//      whose center is "center" and axes are a1, a2 and a3.

int get_param_ell(const bigrational_vector& center,
                  const bigrational_vector& a1,
                  const bigrational_vector& a2,
                  const bigrational_vector& a3,
                  K_RATPOLY*& X, K_RATPOLY*& Y, K_RATPOLY*& Z, K_RATPOLY*& W)
{
  assert(center.get_dim() == 3);
  assert(a1.get_dim() == 3);
  assert(a2.get_dim() == 3);
  assert(a3.get_dim() == 3);
  
  long d[2];
  long p[2];
  
  d[0] = d[1] = 2;
  X    = new K_RATPOLY(2, d);
  Y    = new K_RATPOLY(2, d);
  Z    = new K_RATPOLY(2, d);
  W    = new K_RATPOLY(2, d);
  
  p[0] = 2;
  p[1] = 0;
  X->get_coeff(p) = center[0] - a1[0];
  Y->get_coeff(p) = center[1] - a1[1];
  Z->get_coeff(p) = center[2] - a1[2];
  
  p[0] = 0;
  p[1] = 2;
  X->get_coeff(p) = center[0] - a1[0];
  Y->get_coeff(p) = center[1] - a1[1];
  Z->get_coeff(p) = center[2] - a1[2];
  
  p[0] = 1;
  p[1] = 0;
  X->get_coeff(p) = 2 * a2[0];
  Y->get_coeff(p) = 2 * a2[1];
  Z->get_coeff(p) = 2 * a2[2];
  
  p[0] = 0;
  p[1] = 1;
  X->get_coeff(p) = 2 * a3[0];
  Y->get_coeff(p) = 2 * a3[1];
  Z->get_coeff(p) = 2 * a3[2];
  
  p[0] = 0;
  p[1] = 0;
  X->get_coeff(p) = center[0] + a1[0];
  Y->get_coeff(p) = center[1] + a1[1];
  Z->get_coeff(p) = center[2] + a1[2];
  
  p[0] = 2;
  p[1] = 0;
  W->get_coeff(p) = 1;
  
  p[0] = 0;
  p[1] = 2;
  W->get_coeff(p) = 1;
  
  p[0] = 0;
  p[1] = 0;
  W->get_coeff(p) = 1;
  
  return 0;
}


//  unsigned long get_patches_ell(const bigrational_vector pts[4],
//                                K_PATCH**& patches)
//    computes 8 patches for the ellipsoid
//      whose center is pts[0] and axes are pts[1], pts[2] and pts[3].

unsigned long get_patches_ell(const bigrational_vector pts[4],
                              K_PATCH**& patches)
{
  const unsigned long num_patches = 8;
  
  unsigned long i, j;
  int           non_std;
  bigrational   p12, p13, p30, p11, p22, p33;
  K_RATPOLY*    impl;
  
  //  1.  See if a given ellipsoid is non_standard i.e., it is a sphere.
  
  non_std = 0;
  p12     = in_prod(pts[1], pts[2]);
  
  if (sgn(p12))
    non_std = 1;
  else  //  if (!sgn(p12))
  {
    p13 = in_prod(pts[1], pts[3]);
    
    if (sgn(p13))
      non_std = 1;
    else  //  if (!sgn(p13))
    {
      p30 = in_prod(pts[3], pts[0]);
      
      if (!sgn(p30))
        non_std = 1;
      else  //  if (!sgn(p30))
      {
        p11 = in_prod(pts[1], pts[1]);
        p22 = in_prod(pts[2], pts[2]);
        
        if (p11 != p22)
          non_std = 1;
        else  //  if (p11 == p22)
        {
          p33 = in_prod(pts[3], pts[3]);
          
          if (p11 != p33)
            non_std = 1;
        }
      }
    }
  }
  
  //  !non_std iff !p12 && !p13 && !p30 && p11 == p22 == p33
  
  //  2.  Construct patches.
  
  //  2-1.  Compute the implicit formula for surf.
  
  if (non_std)
    get_impl_ell(pts, impl);
  else  //  if (!non_std)
    get_impl_sphere(pts[0], p11, impl);
  
  bigrational_vector minus_pts[4];
  K_RATPOLY*         X;
  K_RATPOLY*         Y;
  K_RATPOLY*         Z;
  K_RATPOLY*         W;
  K_SURF*            surf;
  long               d[3];
  long               p[3];
  K_RATPOLY*         poly;
  K_SEGMENT**        segment;
  unsigned long      num_trim_curves;
  K_CURVE**          trim_curves;
  
  for (i = 1; i < 4; i++)
    minus_pts[i] = - pts[i];
  
  for (i = 0; i < num_patches; i++)
  {
    //  2-2.  Compute the parametric formulae for surf.
    
    if (i == 0)
      get_param_ell(pts[0], pts[1], pts[2], pts[3], X, Y, Z, W);
    else if (i == 1)
      get_param_ell(pts[0], pts[1], minus_pts[2], pts[3], X, Y, Z, W);
    else if (i == 2)
      get_param_ell(pts[0], minus_pts[1], minus_pts[2], pts[3], X, Y, Z, W);
    else if (i == 3)
      get_param_ell(pts[0], minus_pts[1], pts[2], pts[3], X, Y, Z, W);
    else if (i == 4)
      get_param_ell(pts[0], pts[1], pts[2], minus_pts[3], X, Y, Z, W);
    else if (i == 5)
      get_param_ell(pts[0], pts[1], minus_pts[2], minus_pts[3], X, Y, Z, W);
    else if (i == 6)
      get_param_ell(pts[0], minus_pts[1], minus_pts[2], minus_pts[3],
                    X, Y, Z, W);
    else  //  if (i == 7)
      get_param_ell(pts[0], minus_pts[1], pts[2], minus_pts[3], X, Y, Z, W);
    
    //  2-3.  Set surf.
    
    surf = new K_SURF(impl, X, Y, Z, W);
    
    //  2-4.  Set trim_curves.
    
    trim_curves = new K_CURVE* [num_trim_curves = 3];
    
    if (i == 0 || i == 2 || i == 5 || i == 7)
    {
      //  T = 0
      
      d[0] = 0;
      d[1] = 1;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 0;
      p[1] = 1;
      poly->get_coeff(p) = 1;
      
      segment    = new K_SEGMENT* [1];
      segment[0] = new K_SEGMENT(new K_POINT2D(0, 0), new K_POINT2D(1, 0));
      segment[0]->ref_count++;
      
      trim_curves[0] = new K_CURVE(poly, segment, 1);
      trim_curves[0]->ref_count++;
      
      if (!--segment[0]->ref_count)
        delete segment[0];
      
      delete [] segment;
      
      //  S^2 + T^2 - 1 = 0
      
      d[0] = 2;
      d[1] = 2;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 2;
      p[1] = 0;
      poly->get_coeff(p) = 1;
      
      p[0] = 0;
      p[1] = 2;
      poly->get_coeff(p) = 1;
      
      p[0] = 0;
      p[1] = 0;
      poly->get_coeff(p) = - 1;
      
      segment    = new K_SEGMENT* [1];
      segment[0] = new K_SEGMENT(new K_POINT2D(1, 0), new K_POINT2D(0, 1));
      segment[0]->ref_count++;
      
      trim_curves[1] = new K_CURVE(poly, segment, 1);
      trim_curves[1]->ref_count++;
      
      if (!--segment[0]->ref_count)
        delete segment[0];
      
      delete [] segment;
      
      //  S = 0
      
      d[0] = 1;
      d[1] = 0;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 1;
      p[1] = 0;
      poly->get_coeff(p) = 1;
      
      segment    = new K_SEGMENT* [1];
      segment[0] = new K_SEGMENT(new K_POINT2D(0, 1), new K_POINT2D(0, 0));
      segment[0]->ref_count++;
      
      trim_curves[2] = new K_CURVE(poly, segment, 1);
      trim_curves[2]->ref_count++;
      
      if (!--segment[0]->ref_count)
        delete segment[0];
      
      delete [] segment;
    }
    else  //  if (i == 1 || i == 3 || i == 4 || i == 6)
    {
      //  S = 0
      
      d[0] = 1;
      d[1] = 0;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 1;
      p[1] = 0;
      poly->get_coeff(p) = 1;
      
      segment    = new K_SEGMENT* [1];
      segment[0] = new K_SEGMENT(new K_POINT2D(0, 0), new K_POINT2D(0, 1));
      segment[0]->ref_count++;
      
      trim_curves[0] = new K_CURVE(poly, segment, 1);
      trim_curves[0]->ref_count++;
      
      if (!--segment[0]->ref_count)
        delete segment[0];
      
      delete [] segment;
      
      //  S^2 + T^2 - 1 = 0
      
      d[0] = 2;
      d[1] = 2;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 2;
      p[1] = 0;
      poly->get_coeff(p) = 1;
      
      p[0] = 0;
      p[1] = 2;
      poly->get_coeff(p) = 1;
      
      p[0] = 0;
      p[1] = 0;
      poly->get_coeff(p) = - 1;
      
      segment    = new K_SEGMENT* [1];
      segment[0] = new K_SEGMENT(new K_POINT2D(0, 1), new K_POINT2D(1, 0));
      segment[0]->ref_count++;
      
      trim_curves[1] = new K_CURVE(poly, segment, 1);
      trim_curves[1]->ref_count++;
      
      if (!--segment[0]->ref_count)
        delete segment[0];
      
      delete [] segment;
      
      //  T = 0
      
      d[0] = 0;
      d[1] = 1;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 0;
      p[1] = 1;
      poly->get_coeff(p) = 1;
      
      segment    = new K_SEGMENT* [1];
      segment[0] = new K_SEGMENT(new K_POINT2D(1, 0), new K_POINT2D(0, 0));
      segment[0]->ref_count++;
      
      trim_curves[2] = new K_CURVE(poly, segment, 1);
      trim_curves[2]->ref_count++;
      
      if (!--segment[0]->ref_count)
        delete segment[0];
      
      delete [] segment;
    }
    
    //  2-5.  Construct a patch.
    
    patches[i] = new K_PATCH(surf, trim_curves, num_trim_curves);
    
    for (j = 0; j < num_trim_curves; j++)
      if (!--trim_curves[j]->ref_count)
        delete trim_curves[j];
    
    delete [] trim_curves;
  }
  
  return num_patches;
}

//  int perturb_ell(const bigrational_vector pts[4],
//                  bigrational_vector       perturbed_pts[4],
//                  const bigrational&       factor)
//    perturbs an ellipsoid specified by 4 vectors pts to
//             the ellipsoid specified by 4 vectors perturbed_pts
//      outward by "factor" if factor > 0 and
//      inward  by |factor| if factor < 0.

int perturb_ell(const bigrational_vector pts[4],
                bigrational_vector       perturbed_pts[4],
                const bigrational&       factor)
{
  assert(abs(factor) < 1);
  
  const unsigned long num_pts = 4;
  
  unsigned long       i, j;
  bigrational_vector  center(3);
  bigrational_vector* from_center;
  bigrational         squared_dist_proto, squared_dist;
  
  center = pts[0];
  
  from_center = new bigrational_vector [num_pts - 1];  //  num_pts - 1 == 3
  
  from_center[0] = pts[1];
  squared_dist   = in_prod(from_center[0], from_center[0]);
//  j              = 0;
  
  for (i = 1; i < num_pts; i++)
  {
    from_center[i - 1] = pts[i];
//    squared_dist_proto = in_prod(from_center[i - 1], from_center[i - 1]);
//    
//    if (squared_dist < squared_dist_proto)
//    {
//      squared_dist = squared_dist_proto;
//      j            = i - 1;
//    }
  }
  
//  center = center + scalar_mul(factor * squared_dist, from_center[j]);
  
  perturbed_pts[0] = center;
  
  for (i = 1; i < num_pts; i++)
    perturbed_pts[i] = center + scalar_mul(1 + factor, from_center[i - 1]);
  
  delete [] from_center;  //  num_pts - 1 == 3 => from_center != 0
  
  return 0;
}

//  K_SOLID gen_ell(const bigrational_vector pts[4])
//    returns an ellipsoid specified by 4 vectors pts.

K_SOLID gen_ell(const bigrational_vector pts[4])
{
  const unsigned long num_patches = 8;
  
  unsigned long i;
  K_PATCH**     patches;
  K_SOLID       s;
  
  //  pathces
  
  //  1.  Set patches.
  
  patches = new K_PATCH* [num_patches];
  get_patches_ell(pts, patches);
  
  //  2.  Compute surf_split1, surf_split2 and surf_split3.
  
  bigrational_vector a1(3);
  bigrational_vector a2(3);
  bigrational_vector a3(3);
  bigrational_vector b1(3);
  bigrational_vector b2(3);
  bigrational_vector b3(3);
  bigrational_vector four_pts[4];
  K_RATPOLY*   impl_bilin;
  K_SURF*      surf_split1;
  K_SURF*      surf_split2;
  K_SURF*      surf_split3;
  
  a1 = pts[0] + pts[1];
  a2 = pts[0] + pts[2];
  a3 = pts[0] + pts[3];
  b1 = pts[0] - pts[1];
  b2 = pts[0] - pts[2];
  b3 = pts[0] - pts[3];
  
  four_pts[0] = a1;
  four_pts[1] = a2;
  four_pts[2] = b1;
  four_pts[3] = b2;
  get_impl_plane_bilin(four_pts, impl_bilin);
  surf_split1 = new K_SURF(*impl_bilin);
  delete impl_bilin;
  
  four_pts[0] = a2;
  four_pts[1] = a3;
  four_pts[2] = b2;
  four_pts[3] = b3;
  get_impl_plane_bilin(four_pts, impl_bilin);
  surf_split2 = new K_SURF(*impl_bilin);
  delete impl_bilin;
  
  four_pts[0] = a1;
  four_pts[1] = a3;
  four_pts[2] = b1;
  four_pts[3] = b3;
  get_impl_plane_bilin(four_pts, impl_bilin);
  surf_split3 = new K_SURF(*impl_bilin);
  delete impl_bilin;
  
  //  3.  For each patch, set its adj_surfs and adj_patches.
  
  //  patches[0]:
  
  patches[0]->adj_surfs[0] = surf_split1;
  patches[0]->adj_surfs[0]->ref_count++;
  patches[0]->adj_surfs[1] = surf_split2;
  patches[0]->adj_surfs[1]->ref_count++;
  patches[0]->adj_surfs[2] = surf_split3;
  patches[0]->adj_surfs[2]->ref_count++;
  
  patches[0]->adj_patches[0] = patches[4];
  patches[0]->adj_patches[0]->ref_count++;
  patches[0]->adj_patches[1] = patches[3];
  patches[0]->adj_patches[1]->ref_count++;
  patches[0]->adj_patches[2] = patches[1];
  patches[0]->adj_patches[2]->ref_count++;
  
  //  patches[1]:
  
  patches[1]->adj_surfs[0] = surf_split3;
  patches[1]->adj_surfs[0]->ref_count++;
  patches[1]->adj_surfs[1] = surf_split2;
  patches[1]->adj_surfs[1]->ref_count++;
  patches[1]->adj_surfs[2] = surf_split1;
  patches[1]->adj_surfs[2]->ref_count++;
  
  patches[1]->adj_patches[0] = patches[0];
  patches[1]->adj_patches[0]->ref_count++;
  patches[1]->adj_patches[1] = patches[2];
  patches[1]->adj_patches[1]->ref_count++;
  patches[1]->adj_patches[2] = patches[5];
  patches[1]->adj_patches[2]->ref_count++;
  
  //  patches[2]:
  
  patches[2]->adj_surfs[0] = surf_split1;
  patches[2]->adj_surfs[0]->ref_count++;
  patches[2]->adj_surfs[1] = surf_split2;
  patches[2]->adj_surfs[1]->ref_count++;
  patches[2]->adj_surfs[2] = surf_split3;
  patches[2]->adj_surfs[2]->ref_count++;
  
  patches[2]->adj_patches[0] = patches[6];
  patches[2]->adj_patches[0]->ref_count++;
  patches[2]->adj_patches[1] = patches[1];
  patches[2]->adj_patches[1]->ref_count++;
  patches[2]->adj_patches[2] = patches[3];
  patches[2]->adj_patches[2]->ref_count++;
  
  //  patches[3]:
  
  patches[3]->adj_surfs[0] = surf_split3;
  patches[3]->adj_surfs[0]->ref_count++;
  patches[3]->adj_surfs[1] = surf_split2;
  patches[3]->adj_surfs[1]->ref_count++;
  patches[3]->adj_surfs[2] = surf_split1;
  patches[3]->adj_surfs[2]->ref_count++;
  
  patches[3]->adj_patches[0] = patches[2];
  patches[3]->adj_patches[0]->ref_count++;
  patches[3]->adj_patches[1] = patches[0];
  patches[3]->adj_patches[1]->ref_count++;
  patches[3]->adj_patches[2] = patches[7];
  patches[3]->adj_patches[2]->ref_count++;
  
  //  patches[4]:
  
  patches[4]->adj_surfs[0] = surf_split3;
  patches[4]->adj_surfs[0]->ref_count++;
  patches[4]->adj_surfs[1] = surf_split2;
  patches[4]->adj_surfs[1]->ref_count++;
  patches[4]->adj_surfs[2] = surf_split1;
  patches[4]->adj_surfs[2]->ref_count++;
  
  patches[4]->adj_patches[0] = patches[5];
  patches[4]->adj_patches[0]->ref_count++;
  patches[4]->adj_patches[1] = patches[7];
  patches[4]->adj_patches[1]->ref_count++;
  patches[4]->adj_patches[2] = patches[0];
  patches[4]->adj_patches[2]->ref_count++;
  
  //  patches[5]:
  
  patches[5]->adj_surfs[0] = surf_split1;
  patches[5]->adj_surfs[0]->ref_count++;
  patches[5]->adj_surfs[1] = surf_split2;
  patches[5]->adj_surfs[1]->ref_count++;
  patches[5]->adj_surfs[2] = surf_split3;
  patches[5]->adj_surfs[2]->ref_count++;
  
  patches[5]->adj_patches[0] = patches[1];
  patches[5]->adj_patches[0]->ref_count++;
  patches[5]->adj_patches[1] = patches[6];
  patches[5]->adj_patches[1]->ref_count++;
  patches[5]->adj_patches[2] = patches[4];
  patches[5]->adj_patches[2]->ref_count++;
  
  //  patches[6]:
  
  patches[6]->adj_surfs[0] = surf_split3;
  patches[6]->adj_surfs[0]->ref_count++;
  patches[6]->adj_surfs[1] = surf_split2;
  patches[6]->adj_surfs[1]->ref_count++;
  patches[6]->adj_surfs[2] = surf_split1;
  patches[6]->adj_surfs[2]->ref_count++;
  
  patches[6]->adj_patches[0] = patches[7];
  patches[6]->adj_patches[0]->ref_count++;
  patches[6]->adj_patches[1] = patches[5];
  patches[6]->adj_patches[1]->ref_count++;
  patches[6]->adj_patches[2] = patches[2];
  patches[6]->adj_patches[2]->ref_count++;
  
  //  patches[7]:
  
  patches[7]->adj_surfs[0] = surf_split1;
  patches[7]->adj_surfs[0]->ref_count++;
  patches[7]->adj_surfs[1] = surf_split2;
  patches[7]->adj_surfs[1]->ref_count++;
  patches[7]->adj_surfs[2] = surf_split3;
  patches[7]->adj_surfs[2]->ref_count++;
  
  patches[7]->adj_patches[0] = patches[3];
  patches[7]->adj_patches[0]->ref_count++;
  patches[7]->adj_patches[1] = patches[4];
  patches[7]->adj_patches[1]->ref_count++;
  patches[7]->adj_patches[2] = patches[6];
  patches[7]->adj_patches[2]->ref_count++;
  
  //  4.  Associate edges.
  
  patches[0]->trim_curves[0]->assoc(patches[4]->trim_curves[2], - 1);
  patches[0]->trim_curves[1]->assoc(patches[3]->trim_curves[1], - 1);
  patches[0]->trim_curves[2]->assoc(patches[1]->trim_curves[0], - 1);
  
  patches[1]->trim_curves[2]->assoc(patches[5]->trim_curves[0], - 1);
  patches[1]->trim_curves[1]->assoc(patches[2]->trim_curves[1], - 1);
  
  patches[2]->trim_curves[0]->assoc(patches[6]->trim_curves[2], - 1);
  patches[2]->trim_curves[2]->assoc(patches[3]->trim_curves[0], - 1);
  
  patches[3]->trim_curves[2]->assoc(patches[7]->trim_curves[0], - 1);
  
  patches[4]->trim_curves[1]->assoc(patches[7]->trim_curves[1], - 1);
  patches[4]->trim_curves[0]->assoc(patches[5]->trim_curves[2], - 1);
  
  patches[5]->trim_curves[1]->assoc(patches[6]->trim_curves[1], - 1);
  
  patches[6]->trim_curves[0]->assoc(patches[7]->trim_curves[2], - 1);
  
  for (i = 0; i < num_patches; i++)
  {
    patches[i]->surf->X->reduce_deg();
    patches[i]->surf->Y->reduce_deg();
    patches[i]->surf->Z->reduce_deg();
    patches[i]->surf->W->reduce_deg();
  }
  
  s = K_SOLID(patches, num_patches);
  
  return s;
}

K_SOLID read_ell(istream& in_fs, const bigrational& perturb_factor)
{
  const unsigned long num_pts = 4;
  
  unsigned long      i, j;
//  unsigned long      type;
  bigrational_vector pts[4];
  bigrational_vector perturbed_pts[4];
  K_SOLID            s;
  
//  in_fs >> type;
//  cerr << " genell: read_perturb_ell: type = " << type << endl << flush;
//  assert(type == 3);
  cerr << endl << flush;
  
  //  1.  Read vectors from is_fs.
  
  for (i = 0; i < num_pts; i++)
  {
    pts[i] = bigrational_vector(3);
    
    cerr << " genell: read_ell: pts[" << i << "] = ( " << flush;
    for (j = 0; j < 3; j++)
    {
      in_fs >> pts[i][j];
      
      cerr << pts[i][j];
      if (j < 2)
        cerr << ", ";
    }
    cerr << " )" << endl << flush;
  }
  cerr << endl << flush;
  
  //  2.  Perturb pts if necessarily and generate an ellipsoid.
  
  if (!sgn(perturb_factor))
    s = gen_ell(pts);
  else  //  if (sgn(perturb_factor))
  {
    perturb_ell(pts, perturbed_pts, perturb_factor);
    
    for (i = 0; i < num_pts; i++)
      cerr << " genell: read_ell: perturb_pts[" << i << "] =" << perturbed_pts[i] << endl << flush;
    cerr << endl << flush;
    
    s = gen_ell(perturbed_pts);
  }
  
  return s;
}

K_SOLID read_BRLCAD_ell(istream& in_fs, const bigrational& perturb_factor)
{
  const unsigned long num_pts = 4;
  
  unsigned long      i, j;
//  unsigned long      type;
  float              f_pts[4][3];
  bigrational_vector pts[4];
  bigrational        o;
  bigrational_vector perturbed_pts[4];
  K_SOLID            s;
  
//  in_fs >> type;
//  cerr << " genell: read_perturb_BRLCAD_ell: type = " << type << endl << flush;
//  assert(type == 13);
  cerr << endl << flush;
  
  for (i = 0; i < num_pts; i++)
  {
    cerr << " genell: read_BRALCAD_ell: f_pts[" << i << "] = ( " << flush;
    for (j = 0; j < 3; j++)
    {
      in_fs >> f_pts[i][j];
      
      cerr << f_pts[i][j];
      if (j < 2)
        cerr << ", ";
    }
    cerr << " )" << endl << flush;
  }
  cerr << endl << flush;
  
  for (i = 0; i < num_pts; i++)
  {
    pts[i] = bigrational_vector(3);
    
    for (j = 0; j < 3; j++)
      pts[i][j] = as_bigrational(f_pts[i][j]);
  }
  
  o = (pts[1][1] * pts[2][2] - pts[1][2] * pts[2][1]) * pts[3][0]
    + (pts[1][2] * pts[2][0] - pts[1][0] * pts[2][2]) * pts[3][1]
    + (pts[1][0] * pts[2][1] - pts[1][1] * pts[2][0]) * pts[3][2];
  
  if (sgn(o) < 0)
    pts[3] = scalar_mul(- 1, pts[3]);
  
  for (i = 0; i < num_pts; i++)
    cerr << " genell: read_BRALCAD_ell: pts[" << i << "] =" << pts[i] << endl << flush;
  cerr << endl << flush;
  
  //  2.  Perturb pts if necessarily and generate an ellipsoid.
  
  if (!sgn(perturb_factor))
    s = gen_ell(pts);
  else  //  if (sgn(perturb_factor))
  {
    perturb_ell(pts, perturbed_pts, perturb_factor);
    
    for (i = 0; i < num_pts; i++)
      cerr << " genell: read_BRALCAD_ell: perturbed_pts[" << i << "] =" << perturbed_pts[i] << endl << flush;
    cerr << endl << flush;
    
    s = gen_ell(perturbed_pts);
  }
  
  return s;
}

