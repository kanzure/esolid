#include <config.h>

#include <gentor.h>

#include <fpconversion.h>
#include <genbox.h>

//  int get_param_tor(const bigrational_vector& center,
//                    const bigrational_vector& a0,
//                    const bigrational_vector& a1,
//                    const bigrational_vector& a2,
//                    const bigrational&        r0,
//                    const bigrational&        r1,
//                    K_RATPOLY*&               X,
//                    K_RATPOLY*&               Y,
//                    K_RATPOLY*&               Z,
//                    K_RATPOLY*&               W)
//    computes the parametric formula (*X/*W, *Y/*W, *Z/*W) for the torus
//      specified by (center, ...).

int get_param_tor(const bigrational_vector& center,
                  const bigrational_vector& a0,
                  const bigrational_vector& a1,
                  const bigrational_vector& a2,
                  const bigrational&        r0,
                  const bigrational&        r1,
                  K_RATPOLY*&               X,
                  K_RATPOLY*&               Y,
                  K_RATPOLY*&               Z,
                  K_RATPOLY*&               W)
{
  assert(center.get_dim() == 3);
  assert(a0.get_dim() == 3);
  assert(a1.get_dim() == 3);
  assert(a2.get_dim() == 3);

  long d[2];
  long p[2];

  d[0] = d[1] = 2;

  K_RATPOLY X0(2, d);
  K_RATPOLY Y0(2, d);
  K_RATPOLY Z0(2, d);
  K_RATPOLY W0(2, d);

  p[0] = 2;
  p[1] = 2;
  X0.get_coeff(p) = center[0] - r0 * a2[0] - r1 * a0[0];
  Y0.get_coeff(p) = center[1] - r0 * a2[1] - r1 * a0[1];
  Z0.get_coeff(p) = center[2] - r0 * a2[2] - r1 * a0[2];

  p[0] = 2;
  p[1] = 1;
  X0.get_coeff(p) = 2 * r1 * a2[0];
  Y0.get_coeff(p) = 2 * r1 * a2[1];
  Z0.get_coeff(p) = 2 * r1 * a2[2];

  p[0] = 2;
  p[1] = 0;
  X0.get_coeff(p) = center[0] - r0 * a2[0] + r1 * a0[0];
  Y0.get_coeff(p) = center[1] - r0 * a2[1] + r1 * a0[1];
  Z0.get_coeff(p) = center[2] - r0 * a2[2] + r1 * a0[2];

  p[0] = 1;
  p[1] = 2;
  X0.get_coeff(p) = - 2 * r0 * a1[0];
  Y0.get_coeff(p) = - 2 * r0 * a1[1];
  Z0.get_coeff(p) = - 2 * r0 * a1[2];

  p[0] = 1;
  p[1] = 1;
  X0.get_coeff(p) = 4 * r1 * a1[0];
  Y0.get_coeff(p) = 4 * r1 * a1[1];
  Z0.get_coeff(p) = 4 * r1 * a1[2];

  p[0] = 1;
  p[1] = 0;
  X0.get_coeff(p) = - 2 * r0 * a1[0];
  Y0.get_coeff(p) = - 2 * r0 * a1[1];
  Z0.get_coeff(p) = - 2 * r0 * a1[2];

  p[0] = 0;
  p[1] = 2;
  X0.get_coeff(p) = center[0] + r0 * a2[0] - r1 * a0[0];
  Y0.get_coeff(p) = center[1] + r0 * a2[1] - r1 * a0[1];
  Z0.get_coeff(p) = center[2] + r0 * a2[2] - r1 * a0[2];

  p[0] = 0;
  p[1] = 1;
  X0.get_coeff(p) = - 2 * r1 * a2[0];
  Y0.get_coeff(p) = - 2 * r1 * a2[1];
  Z0.get_coeff(p) = - 2 * r1 * a2[2];

  p[0] = 0;
  p[1] = 0;
  X0.get_coeff(p) = center[0] + r0 * a2[0] + r1 * a0[0];
  Y0.get_coeff(p) = center[1] + r0 * a2[1] + r1 * a0[1];
  Z0.get_coeff(p) = center[2] + r0 * a2[2] + r1 * a0[2];

  p[0] = 2;
  p[1] = 2;
  W0.get_coeff(p) = 1;

  p[0] = 2;
  p[1] = 0;
  W0.get_coeff(p) = 1;

  p[0] = 0;
  p[1] = 2;
  W0.get_coeff(p) = 1;

  p[0] = 0;
  p[1] = 0;
  W0.get_coeff(p) = 1;

  //  Make some substitutions to scale.
  //  Patches lie in [- 1, 1] x [- 1, 1]. Want them to lie in [0, 1] x [0, 1].
  //  Substitute 2 S - 1 for S and 2 T - 1 for T.

  K_RATPOLY X1;
  K_RATPOLY Y1;
  K_RATPOLY Z1;
  K_RATPOLY W1;

  X1 = X0.add_var(1);
  Y1 = Y0.add_var(1);
  Z1 = Z0.add_var(1);
  W1 = W0.add_var(1);

  d[0] = 1;
  d[1] = 0;

  K_RATPOLY sub_s(2, d);

  p[0] = 1;
  p[1] = 0;
  sub_s.get_coeff(p) = 2;

  p[0] = 0;
  p[1] = 0;
  sub_s.get_coeff(p) = - 1;

  K_RATPOLY X2;
  K_RATPOLY Y2;
  K_RATPOLY Z2;
  K_RATPOLY W2;

  X2 = X1.subst_expr(0, sub_s);
  Y2 = Y1.subst_expr(0, sub_s);
  Z2 = Z1.subst_expr(0, sub_s);
  W2 = W1.subst_expr(0, sub_s);

  K_RATPOLY X3;
  K_RATPOLY Y3;
  K_RATPOLY Z3;
  K_RATPOLY W3;

  X3 = X2.add_var(2);
  Y3 = Y2.add_var(2);
  Z3 = Z2.add_var(2);
  W3 = W2.add_var(2);

  d[0] = 0;
  d[1] = 1;

  K_RATPOLY sub_t(2, d);

  p[0] = 0;
  p[1] = 1;
  sub_t.get_coeff(p) = 2;

  p[0] = 0;
  p[1] = 0;
  sub_t.get_coeff(p) = - 1;

  X = new K_RATPOLY(X3.subst_expr(1, sub_t));
  Y = new K_RATPOLY(Y3.subst_expr(1, sub_t));
  Z = new K_RATPOLY(Z3.subst_expr(1, sub_t));
  W = new K_RATPOLY(W3.subst_expr(1, sub_t));

  return 0;
}

//  int get_patch_tor(K_RATPOLY* const impl,
//                    K_RATPOLY* const X,
//                    K_RATPOLY* const Y,
//                    K_RATPOLY* const Z,
//                    K_RATPOLY* const W,
//                    K_PATCH*& patch)
//    computes the patch for the torus
//      whose implicit and parametric formulae are
//              *I and (*X/*W, *Y/*W, *Z/*W), resp.

int get_patch_tor(K_RATPOLY* const impl,
                  K_RATPOLY* const X,
                  K_RATPOLY* const Y,
                  K_RATPOLY* const Z,
                  K_RATPOLY* const W,
                  K_PATCH*& patch)
{
  unsigned long i, j;
  K_SURF*       surf;
  long          d[2];
  long          p[2];
  K_RATPOLY*    poly;
  K_SEGMENT**   segment;
  unsigned long num_trim_curves;
  K_CURVE**     trim_curves;

  //  1.  Set surf.

  surf = new K_SURF(impl, X, Y, Z, W);

  //  2.  Set trim_curves.

  trim_curves = new K_CURVE* [num_trim_curves = 4];

  for (i = 0; i < num_trim_curves; i++)
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

      segment[0] = new K_SEGMENT(new K_POINT2D(0, 0), new K_POINT2D(1, 0));
      segment[0]->ref_count++;
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

      segment[0] = new K_SEGMENT(new K_POINT2D(1, 0), new K_POINT2D(1, 1));
      segment[0]->ref_count++;
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

      segment[0] = new K_SEGMENT(new K_POINT2D(1, 1), new K_POINT2D(0, 1));
      segment[0]->ref_count++;
    }
    else  //  if (i == 3)
    {
      d[0] = 1;
      d[1] = 0;
      poly = new K_RATPOLY(2, d);

      p[0] = 1;
      p[1] = 0;
      poly->get_coeff(p) = 1;

      segment[0] = new K_SEGMENT(new K_POINT2D(0, 1), new K_POINT2D(0, 0));
      segment[0]->ref_count++;
    }

    trim_curves[i] = new K_CURVE(poly, segment, 1);
    trim_curves[i]->ref_count++;

    if (!--segment[0]->ref_count)
      delete segment[0];

    delete [] segment;
  }

  //  3.  Construct the patch.

  patch = new K_PATCH(surf, trim_curves, num_trim_curves);

  for (i = 0; i < num_trim_curves; i++)
    if (!--trim_curves[i]->ref_count)
      delete trim_curves[i];

  delete [] trim_curves;

  return 0;
}

//  int perturb_tor(const bigrational_vector& center,
//                  const bigrational_vector& a0,
//                  const bigrational_vector& a1,
//                  const bigrational_vector& a2,
//                  const bigrational&        r0,
//                  const bigrational&        r1,
//                  bigrational_vector&       perturbed_center,
//                  bigrational_vector&       perturbed_a0,
//                  bigrational_vector&       perturbed_a1,
//                  bigrational_vector&       perturbed_a2,
//                  bigrational_vector&       perturbed_r0,
//                  bigrational_vector&       perturbed_r1,
//                  const bigrational&        factor)
//    perturbs a cylinder specified by vectors (center, ...) to
//             the cylinder specified by vectors (perturbed_center, ...)
//      outward by "factor" if factor > 0 and
//      inward  by |factor| if factor < 0.

int perturb_tor(const bigrational_vector& center,
                const bigrational_vector& a0,
                const bigrational_vector& a1,
                const bigrational_vector& a2,
                const bigrational&        r0,
                const bigrational&        r1,
                bigrational_vector&       perturbed_center,
                bigrational_vector&       perturbed_a0,
                bigrational_vector&       perturbed_a1,
                bigrational_vector&       perturbed_a2,
                bigrational&              perturbed_r0,
                bigrational&              perturbed_r1,
                const bigrational&        factor)
{
  assert(abs(factor) < 1);

  const unsigned long num_pts = 4;

  unsigned long       i, j;
  bigrational_vector* from_center;
  bigrational         squared_dist_proto, squared_dist;

  from_center = new bigrational_vector [num_pts - 1];  //  num_pts - 1 == 3

  from_center[0] = a0 - center;
//  squared_dist   = in_prod(from_center[0], from_center[0]);
//  j              = 0;

  for (i = 1; i < num_pts - 1; i++)
  {
    if (i == 1)
      from_center[i] = a1 - center;
    else  //  if (i == 2)
      from_center[i] = a2 - center;
//
//    squared_dist_proto = in_prod(from_center[i], from_center[i]);
//
//    if (squared_dist < squared_dist_proto)
//    {
//      squared_dist = squared_dist_proto;
//      j            = i;
//    }
  }

  perturbed_center = center;
//  perturbed_center =
//    center + scalar_mul(factor * factor * squared_dist, from_center[j]);
  perturbed_a0 = perturbed_center + scalar_mul(1 + factor, from_center[0]);
  perturbed_a1 = perturbed_center + scalar_mul(1 + factor, from_center[1]);
  perturbed_a2 = perturbed_center + scalar_mul(1 + factor, from_center[2]);
  perturbed_r0 = (1 + factor) * r0;
  perturbed_r1 = (1 + factor) * r1;

  delete [] from_center;  //  num_pts - 1 == 3 => from_center != 0

  return 0;
}

//  K_SOLID gen_tor(const bigrational_vector& center,
//                  const bigrational_vector& a0,
//                  const bigrational_vector& a1,
//                  const bigrational_vector& a2,
//                  const bigrational& r0,
//                  const bigrational& r1)
//    returns a torus specified by (center, ...).

#ifndef ALL_CCW
K_SOLID gen_tor(const bigrational_vector& center,
                const bigrational_vector& a0,
                const bigrational_vector& a1,
                const bigrational_vector& a2,
                const bigrational& r0,
                const bigrational& r1)
{
  assert(center.get_dim() == 3);
  assert(a0.get_dim() == 3);
  assert(a1.get_dim() == 3);
  assert(a2.get_dim() == 3);

  const unsigned long num_patches = 4;

  unsigned long      i;
  K_PATCH**          patches;
  bigrational_vector minus_a0;
  bigrational_vector minus_a1;
  bigrational_vector minus_a2;
  K_RATPOLY*         X;
  K_RATPOLY*         Y;
  K_RATPOLY*         Z;
  K_RATPOLY*         W;
  K_RATPOLY*         impl;
  bigrational_vector center_plus_a0;
  bigrational_vector center_minus_a0;
  bigrational_vector center_plus_a1;
  bigrational_vector center_minus_a1;
  bigrational_vector center_plus_a2;
  bigrational_vector center_minus_a2;
  bigrational_vector four_pts[4];
  K_RATPOLY*         impl_bilin;
  K_SURF*            surf_split1;
  K_SURF*            surf_split2;
  K_SOLID            s;

  //  1.  Set patches.

  patches = new K_PATCH* [num_patches];

  minus_a0 = - a0;
  minus_a1 = - a1;
  minus_a2 = - a2;

  //  patches[0]:  top right

  get_param_tor(center, a0, a1, a2, r0, r1, X, Y, Z, W);
  assert(implicitize(*X, *Y, *Z, *W, 4, impl));
  get_patch_tor(impl, X, Y, Z, W, patches[0]);

  //  patches[1]:  top left

  get_param_tor(center, a0, minus_a1, minus_a2, r0, r1, X, Y, Z, W);
  get_patch_tor(impl, X, Y, Z, W, patches[1]);

  //  patches[2]:  bottom right

  get_param_tor(center, minus_a0, minus_a1, a2, r0, r1, X, Y, Z, W);
  get_patch_tor(impl, X, Y, Z, W, patches[2]);

  //  patches[3]:  bottom left

  get_param_tor(center, minus_a0, a1, minus_a2, r0, r1, X, Y, Z, W);
  get_patch_tor(impl, X, Y, Z, W, patches[3]);

  //  2.  Compute surf_split1 and surf_split2.

  center_plus_a0  = center + a0;
  center_minus_a0 = center - a0;
  center_plus_a1  = center + a1;
  center_minus_a1 = center - a1;
  center_plus_a2  = center + a2;
  center_minus_a2 = center - a2;

  four_pts[0] = center_plus_a1;
  four_pts[1] = center_minus_a1;
  four_pts[2] = center_plus_a2;
  four_pts[3] = center_minus_a2;
  get_impl_plane_bilin(four_pts, impl_bilin);
  surf_split1 = new K_SURF(*impl_bilin);
  delete impl_bilin;

  four_pts[0] = center_plus_a1;
  four_pts[1] = center_minus_a1;
  four_pts[2] = center_plus_a0;
  four_pts[3] = center_minus_a0;
  get_impl_plane_bilin(four_pts, impl_bilin);
  surf_split2 = new K_SURF(*impl_bilin);
  delete impl_bilin;

  //  3.  For each patch, set its adj_surfs and adj_patches.

  //  patches[0]:  top right

  patches[0]->adj_surfs[0] = surf_split1;
  patches[0]->adj_surfs[0]->ref_count++;
  patches[0]->adj_surfs[1] = surf_split2;
  patches[0]->adj_surfs[1]->ref_count++;
  patches[0]->adj_surfs[2] = surf_split1;
  patches[0]->adj_surfs[2]->ref_count++;
  patches[0]->adj_surfs[3] = surf_split2;
  patches[0]->adj_surfs[3]->ref_count++;

  patches[0]->adj_patches[0] = patches[2];
  patches[0]->adj_patches[0]->ref_count++;
  patches[0]->adj_patches[1] = patches[1];
  patches[0]->adj_patches[1]->ref_count++;
  patches[0]->adj_patches[2] = patches[2];
  patches[0]->adj_patches[2]->ref_count++;
  patches[0]->adj_patches[3] = patches[1];
  patches[0]->adj_patches[3]->ref_count++;

  //  patches[1]:  top left

  patches[1]->adj_surfs[0] = surf_split1;
  patches[1]->adj_surfs[0]->ref_count++;
  patches[1]->adj_surfs[1] = surf_split2;
  patches[1]->adj_surfs[1]->ref_count++;
  patches[1]->adj_surfs[2] = surf_split1;
  patches[1]->adj_surfs[2]->ref_count++;
  patches[1]->adj_surfs[3] = surf_split2;
  patches[1]->adj_surfs[3]->ref_count++;

  patches[1]->adj_patches[0] = patches[3];
  patches[1]->adj_patches[0]->ref_count++;
  patches[1]->adj_patches[1] = patches[0];
  patches[1]->adj_patches[1]->ref_count++;
  patches[1]->adj_patches[2] = patches[3];
  patches[1]->adj_patches[2]->ref_count++;
  patches[1]->adj_patches[3] = patches[0];
  patches[1]->adj_patches[3]->ref_count++;

  //  patches[2]:  bottom right

  patches[2]->adj_surfs[0] = surf_split1;
  patches[2]->adj_surfs[0]->ref_count++;
  patches[2]->adj_surfs[1] = surf_split2;
  patches[2]->adj_surfs[1]->ref_count++;
  patches[2]->adj_surfs[2] = surf_split1;
  patches[2]->adj_surfs[2]->ref_count++;
  patches[2]->adj_surfs[3] = surf_split2;
  patches[2]->adj_surfs[3]->ref_count++;

  patches[2]->adj_patches[0] = patches[0];
  patches[2]->adj_patches[0]->ref_count++;
  patches[2]->adj_patches[1] = patches[3];
  patches[2]->adj_patches[1]->ref_count++;
  patches[2]->adj_patches[2] = patches[0];
  patches[2]->adj_patches[2]->ref_count++;
  patches[2]->adj_patches[3] = patches[3];
  patches[2]->adj_patches[3]->ref_count++;

  //  patches[3]:  bottom left

  patches[3]->adj_surfs[0] = surf_split1;
  patches[3]->adj_surfs[0]->ref_count++;
  patches[3]->adj_surfs[1] = surf_split2;
  patches[3]->adj_surfs[1]->ref_count++;
  patches[3]->adj_surfs[2] = surf_split1;
  patches[3]->adj_surfs[2]->ref_count++;
  patches[3]->adj_surfs[3] = surf_split2;
  patches[3]->adj_surfs[3]->ref_count++;

  patches[3]->adj_patches[0] = patches[1];
  patches[3]->adj_patches[0]->ref_count++;
  patches[3]->adj_patches[1] = patches[2];
  patches[3]->adj_patches[1]->ref_count++;
  patches[3]->adj_patches[2] = patches[1];
  patches[3]->adj_patches[2]->ref_count++;
  patches[3]->adj_patches[3] = patches[2];
  patches[3]->adj_patches[3]->ref_count++;

  //  4.  Associate edges.

  patches[0]->trim_curves[0]->assoc(patches[2]->trim_curves[0], - 1);
  patches[0]->trim_curves[1]->assoc(patches[1]->trim_curves[1], - 1);
  patches[0]->trim_curves[2]->assoc(patches[2]->trim_curves[2], - 1);
  patches[0]->trim_curves[3]->assoc(patches[1]->trim_curves[3], - 1);

  patches[3]->trim_curves[0]->assoc(patches[1]->trim_curves[0], - 1);
  patches[3]->trim_curves[1]->assoc(patches[2]->trim_curves[1], - 1);
  patches[3]->trim_curves[2]->assoc(patches[1]->trim_curves[2], - 1);
  patches[3]->trim_curves[3]->assoc(patches[2]->trim_curves[3], - 1);

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
#else
K_SOLID gen_tor(const bigrational_vector& center,
                const bigrational_vector& a0,
                const bigrational_vector& a1,
                const bigrational_vector& a2,
                const bigrational& r0,
                const bigrational& r1)
{
  assert(center.get_dim() == 3);
  assert(a0.get_dim() == 3);
  assert(a1.get_dim() == 3);
  assert(a2.get_dim() == 3);

  const unsigned long num_patches = 4;

  unsigned long      i;
  K_PATCH**          patches;
  bigrational_vector minus_a0;
  bigrational_vector minus_a2;
  K_RATPOLY*         X;
  K_RATPOLY*         Y;
  K_RATPOLY*         Z;
  K_RATPOLY*         W;
  K_RATPOLY*         impl;
  bigrational_vector center_plus_a0;
  bigrational_vector center_minus_a0;
  bigrational_vector center_plus_a1;
  bigrational_vector center_minus_a1;
  bigrational_vector center_plus_a2;
  bigrational_vector center_minus_a2;
  bigrational_vector four_pts[4];
  K_RATPOLY*         impl_bilin;
  K_SURF*            surf_split1;
  K_SURF*            surf_split2;
  K_SOLID            s;

  //  1.  Set patches.

  patches = new K_PATCH* [num_patches];

  minus_a0 = - a0;
  minus_a2 = - a2;

  //  patches[0]:  top right

  get_param_tor(center, a0, a1, a2, r0, r1, X, Y, Z, W);
  assert(implicitize(*X, *Y, *Z, *W, 4, impl));
  get_patch_tor(impl, X, Y, Z, W, patches[0]);

  //  patches[1]:  top left

  get_param_tor(center, a0, a1, minus_a2, r0, r1, X, Y, Z, W);
  get_patch_tor(impl, X, Y, Z, W, patches[1]);

  //  patches[2]:  bottom right

  get_param_tor(center, minus_a0, a1, a2, r0, r1, X, Y, Z, W);
  get_patch_tor(impl, X, Y, Z, W, patches[2]);

  //  patches[3]:  bottom left

  get_param_tor(center, minus_a0, a1, minus_a2, r0, r1, X, Y, Z, W);
  get_patch_tor(impl, X, Y, Z, W, patches[3]);

  //  2.  Compute surf_split1 and surf_split2.

  center_plus_a0  = center + a0;
  center_minus_a0 = center - a0;
  center_plus_a1  = center + a1;
  center_minus_a1 = center - a1;
  center_plus_a2  = center + a2;
  center_minus_a2 = center - a2;

  four_pts[0] = center_plus_a1;
  four_pts[1] = center_minus_a1;
  four_pts[2] = center_plus_a2;
  four_pts[3] = center_minus_a2;
  get_impl_plane_bilin(four_pts, impl_bilin);
  surf_split1 = new K_SURF(*impl_bilin);
  delete impl_bilin;

  four_pts[0] = center_plus_a1;
  four_pts[1] = center_minus_a1;
  four_pts[2] = center_plus_a0;
  four_pts[3] = center_minus_a0;
  get_impl_plane_bilin(four_pts, impl_bilin);
  surf_split2 = new K_SURF(*impl_bilin);
  delete impl_bilin;

  //  3.  For each patch, set its adj_surfs and adj_patches.

  //  patches[0]:  top right

  patches[0]->adj_surfs[0] = surf_split1;
  patches[0]->adj_surfs[0]->ref_count++;
  patches[0]->adj_surfs[1] = surf_split2;
  patches[0]->adj_surfs[1]->ref_count++;
  patches[0]->adj_surfs[2] = surf_split1;
  patches[0]->adj_surfs[2]->ref_count++;
  patches[0]->adj_surfs[3] = surf_split2;
  patches[0]->adj_surfs[3]->ref_count++;

  patches[0]->adj_patches[0] = patches[2];
  patches[0]->adj_patches[0]->ref_count++;
  patches[0]->adj_patches[1] = patches[1];
  patches[0]->adj_patches[1]->ref_count++;
  patches[0]->adj_patches[2] = patches[2];
  patches[0]->adj_patches[2]->ref_count++;
  patches[0]->adj_patches[3] = patches[1];
  patches[0]->adj_patches[3]->ref_count++;

  //  patches[1]:  top left

  patches[1]->adj_surfs[0] = surf_split1;
  patches[1]->adj_surfs[0]->ref_count++;
  patches[1]->adj_surfs[1] = surf_split2;
  patches[1]->adj_surfs[1]->ref_count++;
  patches[1]->adj_surfs[2] = surf_split1;
  patches[1]->adj_surfs[2]->ref_count++;
  patches[1]->adj_surfs[3] = surf_split2;
  patches[1]->adj_surfs[3]->ref_count++;

  patches[1]->adj_patches[0] = patches[3];
  patches[1]->adj_patches[0]->ref_count++;
  patches[1]->adj_patches[1] = patches[0];
  patches[1]->adj_patches[1]->ref_count++;
  patches[1]->adj_patches[2] = patches[3];
  patches[1]->adj_patches[2]->ref_count++;
  patches[1]->adj_patches[3] = patches[0];
  patches[1]->adj_patches[3]->ref_count++;

  //  patches[2]:  bottom right

  patches[2]->adj_surfs[0] = surf_split1;
  patches[2]->adj_surfs[0]->ref_count++;
  patches[2]->adj_surfs[1] = surf_split2;
  patches[2]->adj_surfs[1]->ref_count++;
  patches[2]->adj_surfs[2] = surf_split1;
  patches[2]->adj_surfs[2]->ref_count++;
  patches[2]->adj_surfs[3] = surf_split2;
  patches[2]->adj_surfs[3]->ref_count++;

  patches[2]->adj_patches[0] = patches[0];
  patches[2]->adj_patches[0]->ref_count++;
  patches[2]->adj_patches[1] = patches[3];
  patches[2]->adj_patches[1]->ref_count++;
  patches[2]->adj_patches[2] = patches[0];
  patches[2]->adj_patches[2]->ref_count++;
  patches[2]->adj_patches[3] = patches[3];
  patches[2]->adj_patches[3]->ref_count++;

  //  patches[3]:  bottom left

  patches[3]->adj_surfs[0] = surf_split1;
  patches[3]->adj_surfs[0]->ref_count++;
  patches[3]->adj_surfs[1] = surf_split2;
  patches[3]->adj_surfs[1]->ref_count++;
  patches[3]->adj_surfs[2] = surf_split1;
  patches[3]->adj_surfs[2]->ref_count++;
  patches[3]->adj_surfs[3] = surf_split2;
  patches[3]->adj_surfs[3]->ref_count++;

  patches[3]->adj_patches[0] = patches[1];
  patches[3]->adj_patches[0]->ref_count++;
  patches[3]->adj_patches[1] = patches[2];
  patches[3]->adj_patches[1]->ref_count++;
  patches[3]->adj_patches[2] = patches[1];
  patches[3]->adj_patches[2]->ref_count++;
  patches[3]->adj_patches[3] = patches[2];
  patches[3]->adj_patches[3]->ref_count++;

  //  4.  Associate edges.

  patches[0]->trim_curves[0]->assoc(patches[2]->trim_curves[0], 1);
  patches[0]->trim_curves[1]->assoc(patches[1]->trim_curves[1], 1);
  patches[0]->trim_curves[2]->assoc(patches[2]->trim_curves[2], 1);
  patches[0]->trim_curves[3]->assoc(patches[1]->trim_curves[3], 1);

  patches[3]->trim_curves[0]->assoc(patches[1]->trim_curves[0], 1);
  patches[3]->trim_curves[1]->assoc(patches[2]->trim_curves[1], 1);
  patches[3]->trim_curves[2]->assoc(patches[1]->trim_curves[2], 1);
  patches[3]->trim_curves[3]->assoc(patches[2]->trim_curves[3], 1);

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
#endif

K_SOLID read_tor(istream& in_fs, const bigrational& perturb_factor)
{
  unsigned long      i;
  unsigned long      type;
  bigrational_vector center(3);
  bigrational_vector a0(3);
  bigrational_vector a1(3);
  bigrational_vector a2(3);
  bigrational        r0;
  bigrational        r1;
  bigrational_vector perturbed_center(3);
  bigrational_vector perturbed_a0(3);
  bigrational_vector perturbed_a1(3);
  bigrational_vector perturbed_a2(3);
  bigrational        perturbed_r0;
  bigrational        perturbed_r1;
  K_SOLID            s;

//  in_fs >> type;
//  cerr << " gentor: read_tor: type = " << type << endl << flush;
//  assert(type == 4);
  cerr << endl << flush;

  //  1.  Read vectors from is_fs.

  cerr << " gentor: read_tor: center = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> center[i];

    cerr << center[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;

  cerr << " gentor: read_tor: a0     = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> a0[i];

    cerr << a0[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;

  cerr << " gentor: read_tor: a1     = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> a1[i];

    cerr << a1[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;

  cerr << " gentor: read_tor: a2     = ( " << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> a2[i];

    cerr << a2[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;

  cerr << " gentor: read_tor: r0     = " << flush;
  in_fs >> r0;
  cerr << r0 << endl << flush;

  cerr << " gentor: read_tor: r1     = " << flush;
  in_fs >> r1;
  cerr << r1 << endl << flush;
  cerr << endl << flush;

  //  2.  Perturb pts if necessarily and generate a torus.

  if (!sgn(perturb_factor))
    s = gen_tor(center, a0, a1, a2, r0, r1);
  else  //  if (sgn(perturb_factor))
  {
    perturb_tor(center,
                a0, a1, a2,
                r0, r1,
                perturbed_center,
                perturbed_a0, perturbed_a1, perturbed_a2,
                perturbed_r0, perturbed_r1,
                perturb_factor);

    cerr << " gencyl: read_tor: perturbed_center = " << perturbed_center << endl << flush;
    cerr << " gencyl: read_tor: perturbed_a0     = " << perturbed_a0 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_a1     = " << perturbed_a1 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_a2     = " << perturbed_a2 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_r0     = " << perturbed_r0 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_r1     = " << perturbed_r1 << endl << flush;
    cerr << endl << flush;

    s = gen_tor(perturbed_center,
                perturbed_a0, perturbed_a1, perturbed_a2,
                perturbed_r0, perturbed_r1);
  }

  return s;
}

K_SOLID read_BRLCAD_tor(istream& in_fs, const bigrational& perturb_factor)
{
  unsigned long      i, j;
//  unsigned long      type;
  float              f_center[3];
  float              s_vec[3];
  float              d_vec[3];
  float              n_vec[3];
  float              s_r;
  float              d_r;
  bigrational_vector center(3);
  bigrational_vector a0(3);
  bigrational_vector a1(3);
  bigrational_vector a2(3);
  bigrational        r0;
  bigrational        r1;
  bigrational_vector perturbed_center(3);
  bigrational_vector perturbed_a0(3);
  bigrational_vector perturbed_a1(3);
  bigrational_vector perturbed_a2(3);
  bigrational        perturbed_r0;
  bigrational        perturbed_r1;
  K_SOLID            s;

//  in_fs >> type;
//  cerr << " gentor: read_BRLCAD_tor: type = " << type << endl << flush;
//  assert(type == 14);
  cerr << endl << flush;

  //  1.  Read vectors from is_fs.

  cerr << " gentor: read_BRLCAD_tor: f_center = (" << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> f_center[i];

    cerr << f_center[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;

  cerr << " gentor: read_BRLCAD_tor: s_vec    = (" << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> s_vec[i];

    cerr << s_vec[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;

  cerr << " gentor: read_BRLCAD_tor: d_vec    = (" << flush;
  for (i = 0; i < 3; i++)
  {
    in_fs >> d_vec[i];

    cerr << d_vec[i];
    if (i < 2)
      cerr << ", ";
  }
  cerr << " )" << endl << flush;
  cerr << endl << flush;

  //  2.  Compute axes and radii.

  s_r = sqrt(s_vec[0] * s_vec[0] + s_vec[1] * s_vec[1] + s_vec[2] * s_vec[2]);
  d_r = sqrt(d_vec[0] * d_vec[0] + d_vec[1] * d_vec[1] + d_vec[2] * d_vec[2]);

  for (i = 0; i < 3; i++)
  {
    s_vec[i] /= s_r;
    d_vec[i] /= d_r;
  }

  n_vec[0] = s_vec[1] * d_vec[2] - s_vec[2] * d_vec[1];
  n_vec[1] = s_vec[2] * d_vec[0] - s_vec[0] * d_vec[2];
  n_vec[2] = s_vec[0] * d_vec[1] - s_vec[1] * d_vec[0];

  for (i = 0; i < 3; i++)
  {
    center[i] = as_bigrational(f_center[i]);
    a0[i]     = as_bigrational(s_vec[i]);
    a1[i]     = as_bigrational(d_vec[i]);
    a2[i]     = as_bigrational(n_vec[i]);
  }

  r0 = as_bigrational(d_r);
  r1 = as_bigrational(s_r);

  cerr << " gencyl: read_BRLCAD_tor: center = " << center << endl << flush;
  cerr << " gencyl: read_BRLCAD_tor: a0     = " << a0 << endl << flush;
  cerr << " gencyl: read_BRLCAD_tor: a1     = " << a1 << endl << flush;
  cerr << " gencyl: read_BRLCAD_tor: a2     = " << a2 << endl << flush;
  cerr << " gencyl: read_BRLCAD_tor: r0     = " << r0 << endl << flush;
  cerr << " gencyl: read_BRLCAD_tor: r1     = " << r1 << endl << flush;
  cerr << endl << flush;

  //  2.  Perturb vectors if necessarily and generate a torus.

  if (!sgn(perturb_factor))
    s = gen_tor(center, a0, a1, a2, r0, r1);
  else  //  if (sgn(perturb_factor))
  {
    perturb_tor(center,
                a0, a1, a2,
                r0, r1,
                perturbed_center,
                perturbed_a0, perturbed_a1, perturbed_a2,
                perturbed_r0, perturbed_r1,
                perturb_factor);

    cerr << " gencyl: read_tor: perturbed_center = " << perturbed_center << endl << flush;
    cerr << " gencyl: read_tor: perturbed_a0     = " << perturbed_a0 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_a1     = " << perturbed_a1 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_a2     = " << perturbed_a2 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_r0     = " << perturbed_r0 << endl << flush;
    cerr << " gencyl: read_tor: perturbed_r1     = " << perturbed_r1 << endl << flush;
    cerr << endl << flush;

    s = gen_tor(perturbed_center,
                perturbed_a0, perturbed_a1, perturbed_a2,
                perturbed_r0, perturbed_r1);
  }

  return s;
}

