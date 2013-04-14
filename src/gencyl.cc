//  file:    gencyl.cc
//  update:  01/22/03

//#include <config.h>

#include <gencyl.h>

#include <fpconversion.h>
#include <genbox.h>

//  int get_impl_plane(const bigrational_vector& pts0,
//                     const bigrational_vector& pts1,
//                     const bigrational_vector& pts2,
//                     K_RATPOLY*& I)
//    computes the implicit formula for the plane
//      that passes 3 points pts0, pts1 and pts2.

int get_impl_plane(const bigrational_vector& pts0,
                   const bigrational_vector& pts1,
                   const bigrational_vector& pts2,
                   K_RATPOLY*& I)
{
  assert(pts0.get_dim() == 3);
  assert(pts1.get_dim() == 3);
  assert(pts2.get_dim() == 3);
  
  bigrational* plane;
  long         d[3];
  long         p[3];
  
  get_plane_coeffs(pts0, pts1, pts2, plane);
  d[0] = d[1] = d[2] = 1;
  
  I = new K_RATPOLY(3, d);
  
  p[0] = p[1] = p[2] = 0;
  I->get_coeff(p) = plane[3];
  p[0] = 1;
  p[1] = p[2] = 0;
  I->get_coeff(p) = plane[0];
  p[0] = 0;
  p[1] = 1;
  p[2] = 0;
  I->get_coeff(p) = plane[1];
  p[0] = p[1] = 0;
  p[2] = 1;
  I->get_coeff(p) = plane[2];
  
  delete [] plane;
  
  return 0;
}

//  int get_cyl_cap(const bigrational_vector& center,
//                  const bigrational_vector& vecA,
//                  const bigrational_vector& vecB,
//                  const int is_ccw,
//                  K_PATCH*& cap)
//    computes the ellipse (the cap of a cylinder)
//      whose center is "center" and semi-major/minor axes are vecA/B.
//    trim_curves are oriented counterclockwise (the bottom of the cylinder)
//      if is_ccw == 1 and
//    clockwise (the top of the cylinder) otherwise.

int get_cyl_cap(const bigrational_vector& center,
                const bigrational_vector& vecA,
                const bigrational_vector& vecB,
                const int is_ccw,
                K_PATCH*& cap)
{
  assert(center.get_dim() == 3);
  assert(vecA.get_dim() == 3);
  assert(vecB.get_dim() == 3);
  
  unsigned long      i;
  K_RATPOLY*         I;
  bigrational_vector pts[4];
  K_RATPOLY*         X;
  K_RATPOLY*         Y;
  K_RATPOLY*         Z;
  K_RATPOLY*         W;
  K_SURF*            surf;
  K_RATPOLY*         poly;
  K_SEGMENT**        segment;
  long               d[2];
  long               p[2];
  K_CURVE**          trim_curves;
  unsigned long      num_trim_curves;
  
  //  1.  Set surf.
  
  //  1-1.  Compute the implicit formula for the cap.
  
  get_impl_plane(center, vecA, vecB, I);
  
  //  1-2.  Compute the parametric formulae for the cap.
  
  pts[0] = center - (vecB - center) - (vecA - center);
  pts[1] = center + (vecB - center) - (vecA - center);
  pts[2] = center + (vecB - center) + (vecA - center);
  pts[3] = center - (vecB - center) + (vecA - center);
  
  get_param_plane(pts[0], pts[1], pts[3], X, Y, Z, W);
  
  //  1-3.  Set surf.
  
  surf = new K_SURF(I, X, Y, Z, W);
  
  //  2.  Set trim_curves.
  
  trim_curves = new K_CURVE* [num_trim_curves = 4];
  
  //  All trim_curves have the same polynomial s^2 - s + t^2 - t + 1/4.
  
  d[0] = d[1] = 2;
  poly = new K_RATPOLY(2, d);
  
  p[0] = 2;
  p[1] = 0;
  poly->get_coeff(p) = 1;
  p[0] = 1;
  p[1] = 0;
  poly->get_coeff(p) = - 1;
  p[0] = 0;
  p[1] = 2;
  poly->get_coeff(p) = 1;
  p[0] = 0;
  p[1] = 1;
  poly->get_coeff(p) = - 1;
  p[0] = 0;
  p[1] = 0;
  poly->get_coeff(p) = bigrational(1, 4);
  
  for (i = 0; i < 4; i++)
  {
    segment = new K_SEGMENT* [1];
    
    if (!is_ccw)
      if (i == 0)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(1, bigrational(1, 2)),
                        new K_POINT2D(bigrational(1, 2), 0));
      else if (i == 1)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(bigrational(1, 2), 0),
                        new K_POINT2D(0, bigrational(1, 2)));
      else if (i == 2)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(0, bigrational(1, 2)),
                        new K_POINT2D(bigrational(1, 2), 1));
      else  //  if (i == 3)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(bigrational(1, 2), 1),
                        new K_POINT2D(1, bigrational(1, 2)));
    else  //  if (is_ccw)
      if (i == 0)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(bigrational(1, 2), 0),
                        new K_POINT2D(1, bigrational(1, 2)));
      else if (i == 1)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(1, bigrational(1, 2)),
                        new K_POINT2D(bigrational(1, 2), 1));
      else if (i == 2)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(bigrational(1, 2), 1),
                        new K_POINT2D(0, bigrational(1, 2)));
      else  //  if (i == 3)
        segment[0] =
          new K_SEGMENT(new K_POINT2D(0, bigrational(1, 2)),
                        new K_POINT2D(bigrational(1, 2), 0));
    
    segment[0]->ref_count++;
    
    trim_curves[i] = new K_CURVE(poly, segment, 1);
    trim_curves[i]->ref_count++;
    
    if (!--segment[0]->ref_count)
      delete segment[0];
    
    delete [] segment;
  }
  
  //  3.  Construct the patch.
  
  cap = new K_PATCH(surf, trim_curves, num_trim_curves);
  
  for (i = 0; i < num_trim_curves; i++)
    if (!--trim_curves[i]->ref_count)
      delete trim_curves[i];
  
  delete [] trim_curves;
  
  return 0;
}

//  int get_param_cone(const bigrational_vector& base,
//                     const bigrational_vector& centerline,
//                     const bigrational_vector& base_vecA,
//                     const bigrational_vector& base_vecB,
//                     const bigrational_vector& top_vecA,
//                     const bigrational_vector& top_vecB,
//                     K_RATPOLY*& X,
//                     K_RATPOLY*& Y,
//                     K_RATPOLY*& Z,
//                     K_RATPOLY*& W)
//    computes the parametric formula (*X/*W, *Y/*W, *Z/*W) for the cone
//      whose bottom ellipse is specified by
//              center base, semi-major/minor axes are base_vecA/B and
//            top ellipse is specified by
//              center base+centerline, semi-major/minor axes are top_vecA/B.

int get_param_cone(const bigrational_vector& base,
                   const bigrational_vector& centerline,
                   const bigrational_vector& base_vecA,
                   const bigrational_vector& base_vecB,
                   const bigrational_vector& top_vecA,
                   const bigrational_vector& top_vecB,
                   K_RATPOLY*& X, K_RATPOLY*& Y, K_RATPOLY*& Z, K_RATPOLY*& W)
{
  assert(base.get_dim() == 3);
  assert(centerline.get_dim() == 3);
  assert(base_vecA.get_dim() == 3);
  assert(base_vecB.get_dim() == 3);
  assert(top_vecA.get_dim() == 3);
  assert(top_vecB.get_dim() == 3);
  
  long d[2];
  long p[2];
  
  d[0] = 2;
  d[1] = 1;
  X    = new K_RATPOLY(2, d);
  Y    = new K_RATPOLY(2, d);
  Z    = new K_RATPOLY(2, d);
  d[0] = 2;
  d[1] = 0;
  W    = new K_RATPOLY(2, d);
  
  p[0] = 2;
  p[1] = 1;
  X->get_coeff(p) = centerline[0] + base_vecA[0] - top_vecA[0];
  p[0] = 2;
  p[1] = 0;
  X->get_coeff(p) = base[0] - base_vecA[0];
  p[0] = 1;
  p[1] = 1;
  X->get_coeff(p) = 2 * (top_vecB[0] - base_vecB[0]);
  p[0] = 1;
  p[1] = 0;
  X->get_coeff(p) = 2 * base_vecB[0];
  p[0] = 0;
  p[1] = 1;
  X->get_coeff(p) = centerline[0] - base_vecA[0] + top_vecA[0];
  p[0] = 0;
  p[1] = 0;
  X->get_coeff(p) = base[0] + base_vecA[0];
  
  p[0] = 2;
  p[1] = 1;
  Y->get_coeff(p) = centerline[1] + base_vecA[1] - top_vecA[1];
  p[0] = 2;
  p[1] = 0;
  Y->get_coeff(p) = base[1] - base_vecA[1];
  p[0] = 1;
  p[1] = 1;
  Y->get_coeff(p) = 2 * (top_vecB[1] - base_vecB[1]);
  p[0] = 1;
  p[1] = 0;
  Y->get_coeff(p) = 2 * base_vecB[1];
  p[0] = 0;
  p[1] = 1;
  Y->get_coeff(p) = centerline[1] - base_vecA[1] + top_vecA[1];
  p[0] = 0;
  p[1] = 0;
  Y->get_coeff(p) = base[1] + base_vecA[1];
  
  p[0] = 2;
  p[1] = 1;
  Z->get_coeff(p) = centerline[2] + base_vecA[2] - top_vecA[2];
  p[0] = 2;
  p[1] = 0;
  Z->get_coeff(p) = base[2] - base_vecA[2];
  p[0] = 1;
  p[1] = 1;
  Z->get_coeff(p) = 2 * (top_vecB[2] - base_vecB[2]);
  p[0] = 1;
  p[1] = 0;
  Z->get_coeff(p) = 2 * base_vecB[2];
  p[0] = 0;
  p[1] = 1;
  Z->get_coeff(p) = centerline[2] - base_vecA[2] + top_vecA[2];
  p[0] = 0;
  p[1] = 0;
  Z->get_coeff(p) = base[2] + base_vecA[2];
  
  p[0] = 2;
  p[1] = 0;
  W->get_coeff(p) = 1;
  p[0] = 0;
  p[1] = 0;
  W->get_coeff(p) = 1;
  
  return 0;
}

//  int get_cyl_side(K_RATPOLY* const I,
//                   K_RATPOLY* const X,
//                   K_RATPOLY* const Y,
//                   K_RATPOLY* const Z,
//                   K_RATPOLY* const W,
//                   K_PATCH*& side)
//    computes the cone (the side of a cylinder)
//      whose implicit formula is *I and
//            parametric formula is (*X/*W, *Y/*W, *Z/*W).

int get_cyl_side(K_RATPOLY* const I,
                 K_RATPOLY* const X,
                 K_RATPOLY* const Y,
                 K_RATPOLY* const Z,
                 K_RATPOLY* const W,
                 K_PATCH*& side)
{
  unsigned long i;
  K_SURF*       surf;
  K_CURVE**     trim_curves;
  unsigned long num_trim_curves;
  K_RATPOLY*    poly;
  K_SEGMENT**   segment;
  long          d[2];
  long          p[2];
  
  //  1.  Set surf.
  
  surf = new K_SURF(I, X, Y, Z, W);
  
  //  2.  Set trim_curves.
  
  trim_curves = new K_CURVE* [num_trim_curves = 4];
  
  for (i = 0; i < 4; i++)
  {
    segment = new K_SEGMENT* [1];
    
    if (i == 0)
    {
      d[0] = 0;
      d[1] = 1;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 0;
      p[1] = 1;
      poly->get_coeff(p) = 1;
      
      segment[0] = new K_SEGMENT(new K_POINT2D(0, 0), new K_POINT2D(1, 0), 1);
    }
    else if (i == 1)
    {
      d[0] = 1;
      d[1] = 0;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 1;
      p[1] = 0;
      poly->get_coeff(p) = 1;
      p[0] = 0;
      p[1] = 0;
      poly->get_coeff(p) = - 1;
      
      segment[0] = new K_SEGMENT(new K_POINT2D(1, 0), new K_POINT2D(1, 1), 1);
    }
    else if (i == 2)
    {
      d[0] = 0;
      d[1] = 1;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 0;
      p[1] = 1;
      poly->get_coeff(p) = 1;
      p[0] = 0;
      p[1] = 0;
      poly->get_coeff(p) = - 1;
      
      segment[0] = new K_SEGMENT(new K_POINT2D(1, 1), new K_POINT2D(0, 1), 1);
    }
    else  //  if (i == 3)
    {
      d[0] = 1;
      d[1] = 0;
      poly = new K_RATPOLY(2, d);
      
      p[0] = 1;
      p[1] = 0;
      poly->get_coeff(p) = 1;
      
      segment[0] = new K_SEGMENT(new K_POINT2D(0, 1), new K_POINT2D(0, 0), 1);
    }
    
    segment[0]->ref_count++;
    
    trim_curves[i] = new K_CURVE(poly, segment, 1);
    trim_curves[i]->ref_count++;
    
    if (!--segment[0]->ref_count)
      delete segment[0];
    
    delete [] segment;
  }
  
  //  3.  Construct the patch.
  
   side = new K_PATCH(surf, trim_curves, num_trim_curves);
  
  for (i = 0; i < num_trim_curves; i++)
    if (!--trim_curves[i]->ref_count)
      delete trim_curves[i];
  
  delete [] trim_curves;
  
  return 0;
}

//  int perturb_cyl(const bigrational_vector& base,
//                  const bigrational_vector& centerline,
//                  const bigrational_vector& base_vecA,
//                  const bigrational_vector& base_vecB,
//                  const bigrational_vector& top_vecA,
//                  const bigrational_vector& top_vecB,
//                  bigrational_vector&       perturbed_base,
//                  bigrational_vector&       perturbed_centerline,
//                  bigrational_vector&       perturbed_base_vecA,
//                  bigrational_vector&       perturbed_base_vecB,
//                  bigrational_vector&       perturbed_top_vecA,
//                  bigrational_vector&       perturbed_top_vecB,
//                  const bigrational&        factor)
//    perturbs a cylinder specified by vectors (base, ...) to
//             the cylinder specified by vectors (perturbed_base, ...)
//      outward by "factor" if factor > 0 and
//      inward  by |factor| if factor < 0.

int perturb_cyl(const bigrational_vector& base,
                const bigrational_vector& centerline,
                const bigrational_vector& base_vecA,
                const bigrational_vector& base_vecB,
                const bigrational_vector& top_vecA,
                const bigrational_vector& top_vecB,
                bigrational_vector&       perturbed_base,
                bigrational_vector&       perturbed_centerline,
                bigrational_vector&       perturbed_base_vecA,
                bigrational_vector&       perturbed_base_vecB,
                bigrational_vector&       perturbed_top_vecA,
                bigrational_vector&       perturbed_top_vecB,
                const bigrational&        factor)
{
  assert(abs(factor) < 1);
  
  const unsigned long num_pts = 6;
  
  unsigned long       i, j;
  bigrational_vector  center(3);
  bigrational_vector* from_center;
  bigrational         squared_dist_proto, squared_dist;
  
  center = base + scalar_mul(bigrational(1, 2), centerline);
  
  from_center = new bigrational_vector [num_pts];  //  num_pts == 6
  
  from_center[0] = base - center;
//  squared_dist   = in_prod(from_center[0], from_center[0]);
//  j              = 0;
  
  for (i = 1; i < num_pts; i++)
  {
    if (i == 1)
      from_center[i] = base + centerline - center;
    else if (i == 2)
      from_center[i] = base + base_vecA - center;
    else if (i == 3)
      from_center[i] = base + base_vecB - center;
    else if (i == 4)
      from_center[i] = base + centerline + top_vecA - center;
    else  //  if (i == 5)
      from_center[i] = base + centerline + top_vecB - center;
//    
//    squared_dist_proto = in_prod(from_center[i], from_center[i]);
//    
//    if (squared_dist < squared_dist_proto)
//    {
//      squared_dist = squared_dist_proto;
//      j            = i;
//    }
  }
  
//  center = center + scalar_mul(factor * factor * squared_dist, from_center[j]);
  
  perturbed_base       = center + scalar_mul(1 + factor, from_center[0]);
  perturbed_centerline =
    center + scalar_mul(1 + factor, from_center[1]) - perturbed_base;
  perturbed_base_vecA  =
    center + scalar_mul(1 + factor, from_center[2]) - perturbed_base;
  perturbed_base_vecB  =
    center + scalar_mul(1 + factor, from_center[3]) - perturbed_base;
  perturbed_top_vecA   =
    center + scalar_mul(1 + factor, from_center[4]) - perturbed_base - perturbed_centerline;
  perturbed_top_vecB   =
    center + scalar_mul(1 + factor, from_center[5]) - perturbed_base - perturbed_centerline;
  
  delete [] from_center;  //  num_pts == 6 => from_center != 0
  
  return 0;
}

//  K_SOLID gen_cyl(const bigrational_vector& base,
//                  const bigrational_vector& centerline,
//                  const bigrational_vector& base_vecA,
//                  const bigrational_vector& base_vecB,
//                  const bigrational_vector& top_vecA,
//                  const bigrational_vector& top_vecB)
//    returns a cylinder specified by vectors (base, ...).

K_SOLID gen_cyl(const bigrational_vector& base,
                const bigrational_vector& centerline,
                const bigrational_vector& base_vecA,
                const bigrational_vector& base_vecB,
                const bigrational_vector& top_vecA,
                const bigrational_vector& top_vecB)
{
  assert(base.get_dim() == 3);
  assert(centerline.get_dim() == 3);
  assert(base_vecA.get_dim() == 3);
  assert(base_vecB.get_dim() == 3);
  assert(top_vecA.get_dim() == 3);
  assert(top_vecB.get_dim() == 3);
  
  const unsigned long num_patches = 6;
  
  unsigned long      i;
  bigrational_vector center(3);
  bigrational_vector vecA(3);
  bigrational_vector vecB(3);
  bigrational_vector minus_base_vecA(3);
  bigrational_vector minus_base_vecB(3);
  bigrational_vector minus_top_vecA(3);
  bigrational_vector minus_top_vecB(3);
  K_RATPOLY*         X;
  K_RATPOLY*         Y;
  K_RATPOLY*         Z;
  K_RATPOLY*         W;
  K_RATPOLY*         I;
  bigrational_vector four_pts[4];
  K_RATPOLY*         IB;
  K_SURF*            surf_split1;
  K_SURF*            surf_split2;
  K_PATCH*           patches[6];
  K_SOLID            s;
  
  //  1.  Set patches.
  
  //  pathces[0] -- bottom circle
  
  center = base;
  vecA   = base + base_vecA;
  vecB   = base + base_vecB;
  
  get_cyl_cap(center, vecA, vecB, 1, patches[0]);
  
  //  patches[1], ..., patches[4] -- 4 sides
  
  minus_base_vecA = - base_vecA;
  minus_base_vecB = - base_vecB;
  minus_top_vecA  = - top_vecA;
  minus_top_vecB  = - top_vecB;
  
  get_param_cone(base, centerline,
                 base_vecA, base_vecB,
                 top_vecA, top_vecB,
                 X, Y, Z, W);
  assert(implicitize(*X, *Y, *Z, *W, 2, I));
  get_cyl_side(I, X, Y, Z, W, patches[1]);
  
  get_param_cone(base, centerline,
                 minus_base_vecB, base_vecA,
                 minus_top_vecB,  top_vecA,
                 X, Y, Z, W);
  get_cyl_side(I, X, Y, Z, W, patches[2]);
  
  get_param_cone(base, centerline,
                 minus_base_vecA, minus_base_vecB,
                 minus_top_vecA,  minus_top_vecB,
                 X, Y, Z, W);
  get_cyl_side(I, X, Y, Z, W, patches[3]);
  
  get_param_cone(base, centerline,
                 base_vecB, minus_base_vecA,
                 top_vecB,  minus_top_vecA,
                 X, Y, Z, W);
  get_cyl_side(I, X, Y, Z, W, patches[4]);
  
  //  pathces[5] -- top circle
  
  center = base + centerline;
  vecA   = base + centerline + top_vecA;
  vecB   = base + centerline + top_vecB;
  
  get_cyl_cap(center, vecA, vecB, 0, patches[5]);
  
  //  2.  Compute surf_split1 and surf_split2.
  
  four_pts[0] = base + base_vecA;
  four_pts[1] = base - base_vecA;
  four_pts[2] = base + centerline + top_vecA;
  four_pts[3] = base + centerline - top_vecA;
  get_impl_plane_bilin(four_pts, IB);
  surf_split1 = new K_SURF(*IB);
  
  four_pts[0] = base + base_vecB;
  four_pts[1] = base - base_vecB;
  four_pts[2] = base + centerline + top_vecB;
  four_pts[3] = base + centerline - top_vecB;
  get_impl_plane_bilin(four_pts, IB);
  surf_split2 = new K_SURF(*IB);
  
  //  3.  For each patch, set its adj_surfs and adj_patches.
  
  //  patches[0]:
  
  patches[0]->adj_surfs[0] = patches[4]->surf;
  patches[0]->adj_surfs[0]->ref_count++;
  patches[0]->adj_surfs[1] = patches[1]->surf;
  patches[0]->adj_surfs[1]->ref_count++;
  patches[0]->adj_surfs[2] = patches[2]->surf;
  patches[0]->adj_surfs[2]->ref_count++;
  patches[0]->adj_surfs[3] = patches[3]->surf;
  patches[0]->adj_surfs[3]->ref_count++;
  
  patches[0]->adj_patches[0] = patches[4];
  patches[0]->adj_patches[0]->ref_count++;
  patches[0]->adj_patches[1] = patches[1];
  patches[0]->adj_patches[1]->ref_count++;
  patches[0]->adj_patches[2] = patches[2];
  patches[0]->adj_patches[2]->ref_count++;
  patches[0]->adj_patches[3] = patches[3];
  patches[0]->adj_patches[3]->ref_count++;
  
  //  patches[1]:
  
  patches[1]->adj_surfs[0] = patches[0]->surf;
  patches[1]->adj_surfs[0]->ref_count++;
  patches[1]->adj_surfs[1] = surf_split2;
  patches[1]->adj_surfs[1]->ref_count++;
  patches[1]->adj_surfs[2] = patches[5]->surf;
  patches[1]->adj_surfs[2]->ref_count++;
  patches[1]->adj_surfs[3] = surf_split1;
  patches[1]->adj_surfs[3]->ref_count++;
  
  patches[1]->adj_patches[0] = patches[0];
  patches[1]->adj_patches[0]->ref_count++;
  patches[1]->adj_patches[1] = patches[4];
  patches[1]->adj_patches[1]->ref_count++;
  patches[1]->adj_patches[2] = patches[5];
  patches[1]->adj_patches[2]->ref_count++;
  patches[1]->adj_patches[3] = patches[2];
  patches[1]->adj_patches[3]->ref_count++;
  
  //  patches[2]:
  
  patches[2]->adj_surfs[0] = patches[0]->surf;
  patches[2]->adj_surfs[0]->ref_count++;
  patches[2]->adj_surfs[1] = surf_split1;
  patches[2]->adj_surfs[1]->ref_count++;
  patches[2]->adj_surfs[2] = patches[5]->surf;
  patches[2]->adj_surfs[2]->ref_count++;
  patches[2]->adj_surfs[3] = surf_split2;
  patches[2]->adj_surfs[3]->ref_count++;
  
  patches[2]->adj_patches[0] = patches[0];
  patches[2]->adj_patches[0]->ref_count++;
  patches[2]->adj_patches[1] = patches[1];
  patches[2]->adj_patches[1]->ref_count++;
  patches[2]->adj_patches[2] = patches[5];
  patches[2]->adj_patches[2]->ref_count++;
  patches[2]->adj_patches[3] = patches[3];
  patches[2]->adj_patches[3]->ref_count++;
  
  //  patches[3]:
  
  patches[3]->adj_surfs[0] = patches[0]->surf;
  patches[3]->adj_surfs[0]->ref_count++;
  patches[3]->adj_surfs[1] = surf_split2;
  patches[3]->adj_surfs[1]->ref_count++;
  patches[3]->adj_surfs[2] = patches[5]->surf;
  patches[3]->adj_surfs[2]->ref_count++;
  patches[3]->adj_surfs[3] = surf_split1;
  patches[3]->adj_surfs[3]->ref_count++;
  
  patches[3]->adj_patches[0] = patches[0];
  patches[3]->adj_patches[0]->ref_count++;
  patches[3]->adj_patches[1] = patches[2];
  patches[3]->adj_patches[1]->ref_count++;
  patches[3]->adj_patches[2] = patches[5];
  patches[3]->adj_patches[2]->ref_count++;
  patches[3]->adj_patches[3] = patches[4];
  patches[3]->adj_patches[3]->ref_count++;
  
  //  patches[4]:
  
  patches[4]->adj_surfs[0] = patches[0]->surf;
  patches[4]->adj_surfs[0]->ref_count++;
  patches[4]->adj_surfs[1] = surf_split1;
  patches[4]->adj_surfs[1]->ref_count++;
  patches[4]->adj_surfs[2] = patches[5]->surf;
  patches[4]->adj_surfs[2]->ref_count++;
  patches[4]->adj_surfs[3] = surf_split2;
  patches[4]->adj_surfs[3]->ref_count++;
  
  patches[4]->adj_patches[0] = patches[0];
  patches[4]->adj_patches[0]->ref_count++;
  patches[4]->adj_patches[1] = patches[3];
  patches[4]->adj_patches[1]->ref_count++;
  patches[4]->adj_patches[2] = patches[5];
  patches[4]->adj_patches[2]->ref_count++;
  patches[4]->adj_patches[3] = patches[1];
  patches[4]->adj_patches[3]->ref_count++;
  
  //  patches[5]:
  
  patches[5]->adj_surfs[0] = patches[4]->surf;
  patches[5]->adj_surfs[0]->ref_count++;
  patches[5]->adj_surfs[1] = patches[3]->surf;
  patches[5]->adj_surfs[1]->ref_count++;
  patches[5]->adj_surfs[2] = patches[2]->surf;
  patches[5]->adj_surfs[2]->ref_count++;
  patches[5]->adj_surfs[3] = patches[1]->surf;
  patches[5]->adj_surfs[3]->ref_count++;
  
  patches[5]->adj_patches[0] = patches[4];
  patches[5]->adj_patches[0]->ref_count++;
  patches[5]->adj_patches[1] = patches[3];
  patches[5]->adj_patches[1]->ref_count++;
  patches[5]->adj_patches[2] = patches[2];
  patches[5]->adj_patches[2]->ref_count++;
  patches[5]->adj_patches[3] = patches[1];
  patches[5]->adj_patches[3]->ref_count++;
  
  //  4.  Associate edges.
  
  patches[0]->trim_curves[0]->assoc(patches[4]->trim_curves[0], - 1);
  patches[0]->trim_curves[1]->assoc(patches[1]->trim_curves[0], - 1);
  patches[0]->trim_curves[2]->assoc(patches[2]->trim_curves[0], - 1);
  patches[0]->trim_curves[3]->assoc(patches[3]->trim_curves[0], - 1);
  
  patches[5]->trim_curves[0]->assoc(patches[4]->trim_curves[2], - 1);
  patches[5]->trim_curves[1]->assoc(patches[3]->trim_curves[2], - 1);
  patches[5]->trim_curves[2]->assoc(patches[2]->trim_curves[2], - 1);
  patches[5]->trim_curves[3]->assoc(patches[1]->trim_curves[2], - 1);
  
  patches[1]->trim_curves[3]->assoc(patches[2]->trim_curves[1], - 1);
  patches[2]->trim_curves[3]->assoc(patches[3]->trim_curves[1], - 1);
  patches[3]->trim_curves[3]->assoc(patches[4]->trim_curves[1], - 1);
  patches[4]->trim_curves[3]->assoc(patches[1]->trim_curves[1], - 1);
  
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

K_SOLID read_cyl(istream& in_fs, const bigrational& perturb_factor)
{
  unsigned long      i;
//  unsigned long      type;
  bigrational_vector base(3);
  bigrational_vector centerline(3);
  bigrational_vector base_vecA(3);
  bigrational_vector base_vecB(3);
  bigrational_vector top_vecA(3);
  bigrational_vector top_vecB(3);
  bigrational_vector perturbed_base(3);
  bigrational_vector perturbed_centerline(3);
  bigrational_vector perturbed_base_vecA(3);
  bigrational_vector perturbed_base_vecB(3);
  bigrational_vector perturbed_top_vecA(3);
  bigrational_vector perturbed_top_vecB(3);
  K_SOLID            s;
  
//  in_fs >> type;
//  cerr << " gencyl: read_cyl: type = " << type << endl << flush;
//  assert(type == 2);
  cerr << endl << flush;
   
  //  1.  Read vectors from is_fs.
  
  cerr << " gencyl: read_cyl: base       = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> base[i];
    
    cerr << base[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  
  cerr << " gencyl: read_cyl: centerline = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> centerline[i];
    
    cerr << centerline[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  
  cerr << " gencyl: read_cyl: base_vecA  = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> base_vecA[i];
    
    cerr << base_vecA[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  
  cerr << " gencyl: read_cyl: base_vecB  = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> base_vecB[i];
    
    cerr << base_vecB[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  
  cerr << " gencyl: read_cyl: top_vecA   = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> top_vecA[i];
    
    cerr << top_vecA[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  
  cerr << " gencyl: read_cyl: top_vecB   = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> top_vecB[i];
    
    cerr << top_vecB[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  cerr << endl << flush;
  
  //  2.  Perturb vectors if necessarily and generate a cylinder.
  
  if (!sgn(perturb_factor))
    s = gen_cyl(base, centerline, base_vecA, base_vecB, top_vecA, top_vecB);
  else  //  if (sgn(perturb_factor))
  {
    perturb_cyl(base, centerline,
                base_vecA, base_vecB,
                top_vecA, top_vecB,
                perturbed_base, perturbed_centerline,
                perturbed_base_vecA, perturbed_base_vecB,
                perturbed_top_vecA, perturbed_top_vecB,
                perturb_factor);
    
    cerr << " gencyl: read_cyl: perturbed_base       = " << perturbed_base << endl << flush;
    cerr << " gencyl: read_cyl: perturbed_centerline = " << perturbed_centerline << endl << flush;
    cerr << " gencyl: read_cyl: perturbed_base_vecA  = " << perturbed_base_vecA << endl << flush;
    cerr << " gencyl: read_cyl: perturbed_base_vecB  = " << perturbed_base_vecB << endl << flush;
    cerr << " gencyl: read_cyl: perturbed_top_vecA   = " << perturbed_top_vecA << endl << flush;
    cerr << " gencyl: read_cyl: perturbed_top_vecB   = " << perturbed_top_vecB << endl << flush;
    cerr << endl << flush;
    
    s = gen_cyl(perturbed_base, perturbed_centerline,
                perturbed_base_vecA, perturbed_base_vecB,
                perturbed_top_vecA, perturbed_top_vecB);
  }
  
  return s;
}

K_SOLID read_BRLCAD_cyl(istream& in_fs, const bigrational& perturb_factor)
{
  unsigned long      i;
//  unsigned long      type;
  float              f_base[3];
  float              f_centerline[3];
  float              f_base_vecA[3];
  float              f_base_vecB[3];
  float              size_f_base_vecA, size_f_base_vecB;
  float              f_top_vecA[3];
  float              f_top_vecB[3];
  float              size_f_top_vecA, size_f_top_vecB;
  float              tol;
  int                is_cylinder;
  bigrational_vector base(3);
  bigrational_vector centerline(3);
  bigrational_vector base_vecA(3);
  bigrational_vector base_vecB(3);
  bigrational        size_base_vecA, size_base_vecB;
  bigrational_vector top_vecA(3);
  bigrational_vector top_vecB(3);
  bigrational        size_top_vecA, size_top_vecB;
  bigrational        o;
  bigrational_vector perturbed_base(3);
  bigrational_vector perturbed_centerline(3);
  bigrational_vector perturbed_base_vecA(3);
  bigrational_vector perturbed_base_vecB(3);
  bigrational_vector perturbed_top_vecA(3);
  bigrational_vector perturbed_top_vecB(3);
  K_SOLID            s;
  
//  in_fs >> type;
//  cerr << " gencyl: read_BRLCAD_cyl: type = " << type << endl << flush;
//  assert(type == 12);
  cerr << endl << flush;
  
  //  1.  Read vectors from is_fs.
  
  cerr << " gencyl: read_BRLCAD_cyl: f_base       = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> f_base[i];
    
    cerr << f_base[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  
  cerr << " gencyl: read_BRLCAD_cyl: f_centerline = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> f_centerline[i];
    
    cerr << f_centerline[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  
  for (i = 0; i < 3; i++)
  {
    base[i]       = as_bigrational(f_base[i]);
    centerline[i] = as_bigrational(f_centerline[i]);
  }
  
  size_f_base_vecA = 0.0;
  
  cerr << " gencyl: read_BRLCAD_cyl: f_base_vecA  = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> f_base_vecA[i];
    
    cerr << f_base_vecA[i];
    if (i < 2)
      cerr << ", ";
    
    size_f_base_vecA += f_base_vecA[i] * f_base_vecA[i];
  }
  cerr << " )" << endl << flush;
  
  size_f_base_vecA = sqrt(size_f_base_vecA);
  cerr << endl << " gencyl: read_BRLCAD_cyl: size_f_base_vecA = " << size_f_base_vecA << endl << flush;
  
  size_f_base_vecB = 0.0;
  
  cerr << " gencyl: read_BRLCAD_cyl: f_base_vecB  = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> f_base_vecB[i];
    
    cerr << f_base_vecB[i];
    if (i < 2)
      cerr << ", ";
    
    size_f_base_vecB += f_base_vecB[i] * f_base_vecB[i];
  }
  cerr << " )" << endl << flush;
  
  size_f_base_vecB = sqrt(size_f_base_vecB);
  cerr << endl << " gencyl: read_BRLCAD_cyl: size_f_base_vecB = " << size_f_base_vecB << endl << flush;
  
  size_f_top_vecA = 0.0;
  
  cerr << " gencyl: read_BRLCAD_cyl: f_top_vecA   = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> f_top_vecA[i];
    
    cerr << f_top_vecA[i];
    if (i < 2)
      cerr << ", ";
    
    size_f_top_vecA += f_top_vecA[i] * f_top_vecA[i];
  }
  cerr << " )" << endl << flush;
  
  size_f_top_vecA = sqrt(size_f_top_vecA);
  cerr << endl << " gencyl: read_BRLCAD_cyl: size_f_top_vecA = " << size_f_top_vecA << endl << flush;
  
  size_f_top_vecB = 0.0;
  
  cerr << " gencyl: read_BRLCAD_cyl: f_top_vecB   = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> f_top_vecB[i];
    
    cerr << f_top_vecB[i];
    if (i < 2)
      cerr << ", ";
    
    size_f_top_vecB += f_top_vecB[i] * f_top_vecB[i];
  }
  cerr << " )" << endl << flush;
  
  size_f_top_vecB = sqrt(size_f_top_vecB);
  cerr << endl << " gencyl: read_BRLCAD_cyl: size_f_top_vecB = " << size_f_top_vecB << endl << flush;
  cerr << endl << flush;
  
  tol         = 0.00001;
  is_cylinder = 1;
  
  if (size_f_top_vecA / size_f_base_vecA - size_f_top_vecB / size_f_base_vecB
      > tol)
    is_cylinder = 0;
  
  for (i = 0; i < 3; i++)
  {
    if (f_base_vecA[i] / size_f_base_vecA - f_top_vecA[i] / size_f_top_vecA
        > tol)
      is_cylinder = 0;
    
    if (f_base_vecB[i] / size_f_base_vecB - f_top_vecB[i] / size_f_top_vecB
        > tol)
      is_cylinder = 0;
  }
  
  cerr << " gencyl: read_BRLCAD_cyl: is_cylinder = " << is_cylinder << endl << flush;
  cerr << endl << flush;
  
  if (is_cylinder)
  {
    size_base_vecA = as_bigrational(size_f_base_vecA);
    size_base_vecB = as_bigrational(size_f_base_vecB);
    size_top_vecA  = as_bigrational(size_f_top_vecA);
    size_top_vecB  = size_top_vecA * size_base_vecB / size_base_vecA;
    
    for (i = 0; i < 3; i++)
    {
      base_vecA[i] =
        as_bigrational(f_base_vecA[i] / size_f_base_vecA) * size_base_vecA;
      base_vecB[i] =
        as_bigrational(f_base_vecB[i] / size_f_base_vecB) * size_base_vecB;
      top_vecA[i]  =
        as_bigrational(f_base_vecA[i] / size_f_base_vecA) * size_top_vecA;
      top_vecB[i]  =
        as_bigrational(f_base_vecB[i] / size_f_base_vecB) * size_top_vecB;
    }
  }
  else  //  if (!is_cylinder)
    for (i = 0; i < 3; i++)
    {
      base_vecA[i]  = as_bigrational(f_base_vecA[i]);
      base_vecB[i]  = as_bigrational(f_base_vecB[i]);
      top_vecA[i]   = as_bigrational(f_top_vecA[i]);
      top_vecB[i]   = as_bigrational(f_top_vecB[i]);
    }
  
  o = (base_vecA[1] * base_vecB[2] - base_vecA[2] * base_vecB[1])
    * centerline[0]
    + (base_vecA[2] * base_vecB[0] - base_vecA[0] * base_vecB[2])
    * centerline[1]
    + (base_vecA[0] * base_vecB[1] - base_vecA[1] * base_vecB[0])
    * centerline[2];
  
  if (sgn(o) < 0)
  {
    base_vecA = scalar_mul(- 1, base_vecA);
    top_vecA  = scalar_mul(- 1, top_vecA);
  }
  
  cerr << " gencyl: read_BRLCAD_cyl: base       = " << base << endl << flush;
  cerr << " gencyl: read_BRLCAD_cyl: centerline = " << centerline << endl << flush;
  cerr << " gencyl: read_BRLCAD_cyl: base_vecA  = " << base_vecA << endl << flush;
  cerr << " gencyl: read_BRLCAD_cyl: base_vecB  = " << base_vecB << endl << flush;
  cerr << " gencyl: read_BRLCAD_cyl: top_vecA   = " << top_vecA << endl << flush;
  cerr << " gencyl: read_BRLCAD_cyl: top_vecB   = " << top_vecB << endl << flush;
  cerr << endl << flush;
  
  //  2.  Perturb vectors if necessarily and generate a cylinder.
  
  if (!sgn(perturb_factor))
    s = gen_cyl(base, centerline, base_vecA, base_vecB, top_vecA, top_vecB);
  else  //  if (sgn(perturb_factor))
  {
    perturb_cyl(base, centerline,
                base_vecA, base_vecB,
                top_vecA, top_vecB,
                perturbed_base, perturbed_centerline,
                perturbed_base_vecA, perturbed_base_vecB,
                perturbed_top_vecA, perturbed_top_vecB,
                perturb_factor);
    
    cerr << " gencyl: read_BRLCAD_cyl: perturbed_base       = " << perturbed_base << endl << flush;
    cerr << " gencyl: read_BRLCAD_cyl: perturbed_centerline = " << perturbed_centerline << endl << flush;
    cerr << " gencyl: read_BRLCAD_cyl: perturbed_base_vecA  = " << perturbed_base_vecA << endl << flush;
    cerr << " gencyl: read_BRLCAD_cyl: perturbed_base_vecB  = " << perturbed_base_vecB << endl << flush;
    cerr << " gencyl: read_BRLCAD_cyl: perturbed_top_vecA   = " << perturbed_top_vecA << endl << flush;
    cerr << " gencyl: read_BRLCAD_cyl: perturbed_top_vecB   = " << perturbed_top_vecB << endl << flush;
    cerr << endl << flush;
    
    s = gen_cyl(perturbed_base, perturbed_centerline,
                perturbed_base_vecA, perturbed_base_vecB,
                perturbed_top_vecA, perturbed_top_vecB);
  }
  
  return s;
}

