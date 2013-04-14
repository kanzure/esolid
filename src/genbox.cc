//#include <config.h>

#include <genbox.h>

#include <bigrational_matrix.h>
#include <fpconversion.h>

unsigned long solve(const bigrational_matrix& M,
                    const bigrational_vector& rhs,
                    bigrational_matrix& sol)
{
  long               i, j, k, l, p;
  unsigned long      num_row, num_col;
  unsigned long*     free_var;
  unsigned long      rank, nullity;
  bigrational_vector T;
  bigrational        t;
  unsigned long      num_sol;

  num_row = M.get_num_row();
  num_col = M.get_num_col();
  assert(num_row == rhs.get_dim());

  if (num_col > 0)
  {
    free_var = new unsigned long [num_col];

    for (i = 0; i < num_col; i++)
      free_var[i] = num_col;
  }
  else
    free_var = 0;

  rank = nullity = 0;
  T    = bigrational_vector(num_col);

  //  Gaussian elimination.

  for (i = j = 0; i < num_row; i++)
  {
    //  Find a pivot.

    p = num_row;

    while (p == num_row && j < num_col)
    {
      for (p = i; p < num_row && !sgn(M(p, j)); p++)
        ;

      if (p == num_row)  //  if there is no pivot in column j
        free_var[j++] = nullity++;
    }

    if (p < num_row)  //  if M(p, j) is a pivot
    {
      rank = i + 1;

      if (i < p)
      {
        //  Swap row i and row p.

        for (l = 0; l < num_col; l++)
        {
          T[l]    = M(i, l);
          M(i, l) = M(p, l);
          M(p, l) = T[l];
        }

        t      = rhs[i];
        rhs[i] = rhs[p];
        rhs[p] = t;
      }

      //  Eliminate the rest of column j.

      for (k = i + 1; k < num_row; k++)
      {
        for (l = j + 1; l < num_col; l++)
          M(k, l) -= M(i, l) * M(k, j) / M(i, j);

        rhs[k] -= rhs[i] * M(k, j) / M(i, j);
        M(k, j) = 0;
      }

      j++;
    }
  }

  //  The remaining columns are free.

  for (; j < num_col; j++)
    free_var[j] = nullity++;

  //  Assert homomorphism theorem.

//  cerr << " genbox: solve: rank = " << rank << ", nullity = " << nullity << endl << flush;
  assert(rank + nullity == num_col);

  //  See if there exists a solution.

  for (i = rank; i < num_row && !sgn(rhs[i]); i++)
    ;

  if (i < num_row)
    num_sol = 0;
  else  //  if (i == num_row)
    num_sol = nullity + 1;

  if (num_sol > 0)
  {
    //  Back substitution.

    sol = bigrational_matrix(num_col, num_sol);

    for (l = 0; l < num_sol; l++)
    {
      i = rank - 1;

      for (k = num_col - 1; k >= 0; k--)
        if (free_var[k] == num_col)
        {
          sol(k, l) = rhs[i];

          for (j = k + 1; j < num_col; j++)
            sol(k, l) -= M(i, j) * sol(j, l);

          sol(k, l) /= M(i, k);
          i--;
        }
        else  //  if (free_var[k] < num_col)
          if (l == free_var[k] + 1)
            sol(k, l) = 1;
          else  //  if (l != free_var[k] + 1)
            sol(k, l) = 0;

      assert(i == - 1);
    }
  }

//  cerr << " genbox: solve: num_sol = " << num_sol << endl << flush;

  if (num_col > 0)
    delete [] free_var;

  return num_sol;
}

//  int intersect_planes(const bigrational   plane0[4],
//                       const bigrational   plane1[4],
//                       const bigrational   plane2[4],
//                       bigrational_vector& pt)
//    computes an intersection pt of 3 planes plane0, plane1 and plane2.
//    returns 1 if they intersect and
//            0 otherwise.

int intersect_planes(const bigrational   plane0[4],
                     const bigrational   plane1[4],
                     const bigrational   plane2[4],
                     bigrational_vector& pt)
{
  unsigned long      i;
  bigrational_matrix M(3, 3);
  bigrational_vector rhs(3);
  bigrational_matrix sol;
  int                intersect;

  for (i = 0; i < 3; i++)
  {
    M(0, i) = plane0[i];
    M(1, i) = plane1[i];
    M(2, i) = plane2[i];
  }

  rhs[0] = - plane0[3];
  rhs[1] = - plane1[3];
  rhs[2] = - plane2[3];

  pt = bigrational_vector(3);

  if (solve(M, rhs, sol) > 0)
  {
    for (i = 0; i < 3; i++)
      pt[i] = sol(i, 0);

    intersect = 1;
  }
  else  //  if there is no solution
  {
    for (i = 0; i < 3; i++)
      pt[i] = 0;

    intersect = 0;
  }

  return intersect;
}

//  unsigned long get_impl_coeffs_bilin(const bigrational param_coeffs[12],
//                                      bigrational*&     impl_coeffs)
//    computes the coefficients of the implicit formula for
//    the bilinear surface
//      X(s, t) = param_coeffs[0] +  param_coeffs[1] s
//                                +  param_coeffs[2] t +  param_coeffs[3] st,
//      Y(s, t) = param_coeffs[4] +  param_coeffs[5] s
//                                +  param_coeffs[6] t +  param_coeffs[7] st,
//      Z(s, t) = param_coeffs[8] +  param_coeffs[9] s
//                                + param_coeffs[10] t + param_coeffs[11] st.
//    returns the number of the coefficients of the implicit formula.

unsigned long get_impl_coeffs_bilin(const bigrational param_coeffs[12],
                                    bigrational*&     impl_coeffs)
{
  const unsigned long num_impl_coeffs = 10;

  bigrational x0, x1, x2, x3;
  bigrational y0, y1, y2, y3;
  bigrational z0, z1, z2, z3;
  bigrational cx00, cx01, cx02, cx03, cx11, cx12, cx13, cx22, cx23, cx33;
  bigrational cy00, cy01, cy02, cy03, cy11, cy12, cy13, cy22, cy23, cy33;
  bigrational cz00, cz01, cz02, cz03, cz11, cz12, cz13, cz22, cz23, cz33;

  //  X(s, t) = x0 + x1 s + x2 t + x3 st,
  //  Y(s, t) = y0 + y1 s + y2 t + y3 st,
  //  Z(s, t) = z0 + z1 s + z2 t + z3 st,

  x0 = param_coeffs[0];
  x1 = param_coeffs[1];
  x2 = param_coeffs[2];
  x3 = param_coeffs[3];
  y0 = param_coeffs[4];
  y1 = param_coeffs[5];
  y2 = param_coeffs[6];
  y3 = param_coeffs[7];
  z0 = param_coeffs[8];
  z1 = param_coeffs[9];
  z2 = param_coeffs[10];
  z3 = param_coeffs[11];

  cx01 = x0 * x1;
  cx02 = x0 * x2;
  cx03 = x0 * x3;
  cx12 = x1 * x2;
  cx13 = x1 * x3;
  cx23 = x2 * x3;
  cy01 = y0 * y1;
  cy02 = y0 * y2;
  cy03 = y0 * y3;
  cy12 = y1 * y2;
  cy13 = y1 * y3;
  cy23 = y2 * y3;
  cz01 = z0 * z1;
  cz02 = z0 * z2;
  cz03 = z0 * z3;
  cz12 = z1 * z2;
  cz13 = z1 * z3;
  cz23 = z2 * z3;

  cx00 = x0 * x0;
  cx11 = x1 * x1;
  cx22 = x2 * x2;
  cx33 = x3 * x3;
  cy00 = y0 * y0;
  cy11 = y1 * y1;
  cy22 = y2 * y2;
  cy33 = y3 * y3;
  cz00 = z0 * z0;
  cz11 = z1 * z1;
  cz22 = z2 * z2;
  cz33 = z3 * z3;

  impl_coeffs = new bigrational [num_impl_coeffs];

  // constant term

  impl_coeffs[0] =
      cx00 * (cy12 * cz33 - cy13 * cz23
              + cy33 * cz12 - cy23 * cz13)
    + cx11 * (cy02 * cz23 - cy03 * cz22
              + cy23 * cz02 - cy22 * cz03)
    + cx22 * (cy01 * cz13 - cy03 * cz11
              + cy13 * cz01 - cy11 * cz03)
    + cx33 * (cy00 * cz12 - cy01 * cz02
              + cy12 * cz00 - cy02 * cz01)
    + cx01 * (- cy02 * cz33 + cy03 * cz23 - cy12 * cz23 + cy13 * cz22
              - cy33 * cz02 + cy23 * cz03 - cy23 * cz12 + cy22 * cz13)
    + cx02 * (- cy01 * cz33 + cy03 * cz13 - cy12 * cz13 + cy23 * cz11
              - cy33 * cz01 + cy13 * cz03 - cy13 * cz12 + cy11 * cz23)
    + cx13 * (- cy00 * cz23 + cy01 * cz22 - cy02 * cz12 + cy02 * cz03
              - cy23 * cz00 + cy22 * cz01 - cy12 * cz02 + cy03 * cz02)
    + cx23 * (- cy00 * cz13 + cy02 * cz11 - cy01 * cz12 + cy01 * cz03
              - cy13 * cz00 + cy11 * cz02 - cy12 * cz01 + cy03 * cz01)
    + 2 * cx03 * (- cy03 * cz12 - cy12 * cz03 + cy12 * cz12)
    + cx03 * (cy01 * cz23 + cy02 * cz13 - cy11 * cz22
              + cy23 * cz01 + cy13 * cz02 - cy22 * cz11)
    + 2 * cx12 * (- cy03 * cz03 + cy03 * cz12 + cy12 * cz03)
    + cx12 * (cy00 * cz33 - cy01 * cz23 - cy02 * cz13
              + cy33 * cz00 - cy23 * cz01 - cy13 * cz02);

  // X coefficient

  impl_coeffs[1] =
    2 * x0 * (- cy33 * cz12 + cy23 * cz13
              - cy12 * cz33 + cy13 * cz23)
    + x1 * (cy33 * cz02 + cy23 * cz12 - cy23 * cz03 - cy22 * cz13
            + cy02 * cz33 + cy12 * cz23 - cy03 * cz23 - cy13 * cz22)
    + x2 * (cy33 * cz01 + cy13 * cz12- cy13 * cz03 - cy11 * cz23
            + cy01 * cz33 + cy12 * cz13 - cy03 * cz13 - cy23 * cz11)
    + 2 * x3 * (- cy12 * cz12 + cy12 * cz03 + cy03 * cz12)
    + x3 * (- cy23 * cz01 + cy22 * cz11 - cy13 * cz02
            - cy01 * cz23 + cy11 * cz22 - cy02 * cz13);

  //  Y coefficient

  impl_coeffs[2] =
    2 * y0 * (- cx33 * cz12 + cx23 * cz13
              - cx12 * cz33 + cx13 * cz23)
    + y1 * (cx33 * cz02 + cx23 * cz12 - cx23 * cz03 - cx22 * cz13
            + cx02 * cz33 + cx12 * cz23 - cx03 * cz23 - cx13 * cz22)
    + y2 * (cx33 * cz01 + cx13 * cz12- cx13 * cz03 - cx11 * cz23
            + cx01 * cz33 + cx12 * cz13 - cx03 * cz13 - cx23 * cz11)
    + 2 * y3 * (- cx12 * cz12 + cx12 * cz03 + cx03 * cz12)
    + y3 * (- cx23 * cz01 + cx22 * cz11 - cx13 * cz02
            - cx01 * cz23 + cx11 * cz22 - cx02 * cz13);

  //  Z coefficient

  impl_coeffs[3] =
    2 * z0 * (- cx33 * cy12 + cx23 * cy13
              - cx12 * cy33 + cx13 * cy23)
    + z1 * (cx33 * cy02 + cx23 * cy12 - cx23 * cy03 - cx22 * cy13
            + cx02 * cy33 + cx12 * cy23 - cx03 * cy23 - cx13 * cy22)
    + z2 * (cx33 * cy01 + cx13 * cy12- cx13 * cy03 - cx11 * cy23
            + cx01 * cy33 + cx12 * cy13 - cx03 * cy13 - cx23 * cy11)
    + 2 * z3 * (- cx12 * cy12 + cx12 * cy03 + cx03 * cy12)
    + z3 * (- cx23 * cy01 + cx22 * cy11 - cx13 * cy02
            - cx01 * cy23 + cx11 * cy22 - cx02 * cy13);

  //  XY coefficient

  impl_coeffs[4] =
    - 2 * cz12 * x3 * y3
    + cz13 * (x2 * y3 + x3 * y2)
    + cz23 * (x1 * y3 + x3 * y1)
    - cz33 * (x1 * y2 + x2 * y1);

  //  YZ coefficient

  impl_coeffs[5] =
    - 2 * cx12 * y3 * z3
    + cx13 * (z2 * y3 + z3 * y2)
    + cx23 * (z1 * y3 + z3 * y1)
    - cx33 * (z1 * y2 + z2 * y1);

  //  XZ coefficient

  impl_coeffs[6] =
    - 2 * cy12 * x3 * z3
    + cy13 * (x2 * z3 + x3 * z2)
    + cy23 * (x1 * z3 + x3 * z1)
    - cy33 * (x1 * z2 + x2 * z1);

  //  X^2 coefficient

  impl_coeffs[7] = cy33 * cz12 - cy23 * cz13 + cy12 * cz33 - cy13 * cz23;

  //  Y^2 coefficient

  impl_coeffs[8] = cx33 * cz12 - cx23 * cz13 + cx12 * cz33 - cx13 * cz23;

  //  Z^2 coefficient

  impl_coeffs[9] = cx33 * cy12 - cx23 * cy13 + cx12 * cy33 - cx13 * cy23;

  return num_impl_coeffs;
}

//  int get_impl_plane_bilin(const bigrational_vector pts[4], K_RATPOLY*& impl)
//    computes the implicit formula for the plane of the bilinear surface
//      that passes 4 points pts.
//    returns 1 if 4 points pts are coplanar and
//            0 otherwise, i.e., 4 points pts are on some bilinear surface.

int get_impl_plane_bilin(const bigrational_vector pts[4], K_RATPOLY*& impl)
{
  bigrational* plane;
  bigrational  s;
  long         d[3];
  long         p[3];
  int          is_plane;

  get_plane_coeffs(pts[0], pts[1], pts[2], plane);
  s = plane[0] * pts[3][0] + plane[1] * pts[3][1] + plane[2] * pts[3][2]
      + plane[3];

  if (!sgn(s))
  {
    d[0] = d[1] = d[2] = 1;
    impl = new K_RATPOLY(3, d);

    p[0] = p[1] = p[2] = 0;
    impl->get_coeff(p) = plane[3];
    p[0] = 1;
    p[1] = p[2] = 0;
    impl->get_coeff(p) = plane[0];
    p[0] = 0;
    p[1] = 1;
    p[2] = 0;
    impl->get_coeff(p) = plane[1];
    p[0] = p[1] = 0;
    p[2] = 1;
    impl->get_coeff(p) = plane[2];

    impl->reduce_deg();
    impl->reduce_num_coeffs();

    is_plane = 1;
  }
  else  //  if (sgn(s))
  {
    bigrational  param_coeffs[12];
    bigrational* impl_coeffs;

    param_coeffs[0]  = pts[0][0];
    param_coeffs[1]  = pts[1][0] - pts[0][0];
    param_coeffs[2]  = pts[3][0] - pts[0][0];
    param_coeffs[3]  = pts[2][0] - pts[3][0] - param_coeffs[1];
    param_coeffs[4]  = pts[0][1];
    param_coeffs[5]  = pts[1][1] - pts[0][1];
    param_coeffs[6]  = pts[3][1] - pts[0][1];
    param_coeffs[7]  = pts[2][1] - pts[3][1] - param_coeffs[5];
    param_coeffs[8]  = pts[0][2];
    param_coeffs[9]  = pts[1][2] - pts[0][2];
    param_coeffs[10] = pts[3][2] - pts[0][2];
    param_coeffs[11] = pts[2][2] - pts[3][2] - param_coeffs[9];

    assert(get_impl_coeffs_bilin(param_coeffs, impl_coeffs) == 10);

    d[0] = d[1] = d[2] = 2;
    impl = new K_RATPOLY(3, d);

    p[0] = p[1] = p[2] = 0;
    impl->get_coeff(p) = impl_coeffs[0];
    p[0] = 1;
    p[1] = p[2] = 0;
    impl->get_coeff(p) = impl_coeffs[1];
    p[0] = 0;
    p[1] = 1;
    p[2] = 0;
    impl->get_coeff(p) = impl_coeffs[2];
    p[0] = p[1] = 0;
    p[2] = 1;
    impl->get_coeff(p) = impl_coeffs[3];
    p[0] = p[1] = 1;
    p[2] = 0;
    impl->get_coeff(p) = impl_coeffs[4];
    p[0] = 0;
    p[1] = p[2] = 1;
    impl->get_coeff(p) = impl_coeffs[5];
    p[0] = 1;
    p[1] = 0;
    p[2] = 1;
    impl->get_coeff(p) = impl_coeffs[6];
    p[0] = 2;
    p[1] = p[2] = 0;
    impl->get_coeff(p) = impl_coeffs[7];
    p[0] = 0;
    p[1] = 2;
    p[2] = 0;
    impl->get_coeff(p) = impl_coeffs[8];
    p[0] = p[1] = 0;
    p[2] = 2;
    impl->get_coeff(p) = impl_coeffs[9];

    delete [] impl_coeffs;

    impl->reduce_deg();
    impl->reduce_num_coeffs();

    is_plane = 0;
  }

  return is_plane;
}

//  int get_param_plane(const bigrational_vector& x,
//                      const bigrational_vector& y,
//                      const bigrational_vector& z,
//                      K_RATPOLY*& X,
//                      K_RATPOLY*& Y,
//                      K_RATPOLY*& Z,
//                      K_RATPOLY*& W)
//    computes the parametric formula (*X/*W, *Y/*W, *Z/*W) for the plane
//      that passes 3 points x, y and z.

int get_param_plane(const bigrational_vector& x,
                    const bigrational_vector& y,
                    const bigrational_vector& z,
                    K_RATPOLY*& X,
                    K_RATPOLY*& Y,
                    K_RATPOLY*& Z,
                    K_RATPOLY*& W)
{
  assert(x.get_dim() == 3);
  assert(y.get_dim() == 3);
  assert(z.get_dim() == 3);

  long d[2];
  long p[2];

  d[0] = d[1] = 1;

  X = new K_RATPOLY(2, d);
  Y = new K_RATPOLY(2, d);
  Z = new K_RATPOLY(2, d);
  W = new K_RATPOLY(2, d);

  p[0] = p[1] = 0;
  X->get_coeff(p) = x[0];
  p[0] = 1;
  p[1] = 0;
  X->get_coeff(p) = y[0] - x[0];
  p[0] = 0;
  p[1] = 1;
  X->get_coeff(p) = z[0] - x[0];

  p[0] = p[1] = 0;
  Y->get_coeff(p) = x[1];
  p[0] = 1;
  p[1] = 0;
  Y->get_coeff(p) = y[1] - x[1];
  p[0] = 0;
  p[1] = 1;
  Y->get_coeff(p) = z[1] - x[1];

  p[0] = p[1] = 0;
  Z->get_coeff(p) = x[2];
  p[0] = 1;
  p[1] = 0;
  Z->get_coeff(p) = y[2] - x[2];
  p[0] = 0;
  p[1] = 1;
  Z->get_coeff(p) = z[2] - x[2];

  p[0] = p[1] = 0;
  W->get_coeff(p) = 1;

  return 0;
}

//  int get_param_bilin(const bigrational_vector pts[4],
//                      K_RATPOLY*& X,
//                      K_RATPOLY*& Y,
//                      K_RATPOLY*& Z,
//                      K_RATPOLY*& W)
//    computes the parametric formula
//      (*X/*W, *Y/*W, *Z/*W)
//    for the bilinear surface that passes 4 points pts.

int get_param_bilin(const bigrational_vector pts[4],
                    K_RATPOLY*& X,
                    K_RATPOLY*& Y,
                    K_RATPOLY*& Z,
                    K_RATPOLY*& W)
{
  long d[2];
  long p[2];

  d[0] = d[1] = 1;

  X = new K_RATPOLY(2, d);
  Y = new K_RATPOLY(2, d);
  Z = new K_RATPOLY(2, d);
  W = new K_RATPOLY(2, d);

  p[0] = p[1] = 0;
  X->get_coeff(p) = pts[0][0];
  p[0] = 1;
  p[1] = 0;
  X->get_coeff(p) = pts[1][0] - pts[0][0];
  p[0] = 0;
  p[1] = 1;
  X->get_coeff(p) = pts[3][0] - pts[0][0];
  p[0] = p[1] = 1;
  X->get_coeff(p) = pts[2][0] + pts[0][0] - pts[3][0] - pts[1][0];

  p[0] = p[1] = 0;
  Y->get_coeff(p) = pts[0][1];
  p[0] = 1;
  p[1] = 0;
  Y->get_coeff(p) = pts[1][1] - pts[0][1];
  p[0] = 0;
  p[1] = 1;
  Y->get_coeff(p) = pts[3][1] - pts[0][1];
  p[0] = p[1] = 1;
  Y->get_coeff(p) = pts[2][1] + pts[0][1] - pts[3][1] - pts[1][1];

  p[0] = p[1] = 0;
  Z->get_coeff(p) = pts[0][2];
  p[0] = 1;
  p[1] = 0;
  Z->get_coeff(p) = pts[1][2] - pts[0][2];
  p[0] = 0;
  p[1] = 1;
  Z->get_coeff(p) = pts[3][2] - pts[0][2];
  p[0] = p[1] = 1;
  Z->get_coeff(p) = pts[2][2] + pts[0][2] - pts[3][2] - pts[1][2];

  p[0] = p[1] = 0;
  W->get_coeff(p) = 1;

  return 0;
}

//  int get_patch4(const bigrational_vector pts[4], K_PATCH*& patch)
//    computes the planar or bilinear patch that passes 4 points pts.

int get_patch4(const bigrational_vector pts[4], K_PATCH*& patch)
{
  unsigned long      i, j;
  int                is_plane;
  K_RATPOLY*         impl;
  K_RATPOLY*         X;
  K_RATPOLY*         Y;
  K_RATPOLY*         Z;
  K_RATPOLY*         W;
  K_SURF*            surf;
  bigrational_matrix M(3, 3);
  bigrational_vector v(2);
  bigrational        t;
  long               d[2];
  long               p[2];
  K_RATPOLY          P;
  K_POINT2D*         start;
  K_POINT2D*         end;
  K_SEGMENT**        segment;
  K_CURVE**          trim_curves;
  unsigned long      num_trim_curves;

  if (is_plane = get_impl_plane_bilin(pts, impl))
    get_param_plane(pts[0], pts[1], pts[3], X, Y, Z, W);
  else  //  if (!is_plane)
    get_param_bilin(pts, X, Y, Z, W);

  surf        = new K_SURF(impl, X, Y, Z, W);
  trim_curves = new K_CURVE* [num_trim_curves = 4];

  if (is_plane)
  {
    //  A planar patch.

    p[0] = 1; p[1] = 0;
    M(0, 0) = X->get_coeff(p);
    p[0] = 0; p[1] = 1;
    M(0, 1) = X->get_coeff(p);
    p[0] = 0; p[1] = 0;
    M(0, 2) = X->get_coeff(p) - pts[2][0];

    p[0] = 1; p[1] = 0;
    M(1, 0) = Y->get_coeff(p);
    p[0] = 0; p[1] = 1;
    M(1, 1) = Y->get_coeff(p);
    p[0] = 0; p[1] = 0;
    M(1, 2) = Y->get_coeff(p) - pts[2][1];

    p[0] = 1; p[1] = 0;
    M(2, 0) = Z->get_coeff(p);
    p[0] = 0; p[1] = 1;
    M(2, 1) = Z->get_coeff(p);
    p[0] = 0; p[1] = 0;
    M(2, 2) = Z->get_coeff(p) - pts[2][2];

    assert(sgn(M(0, 0)) || sgn(M(1, 0)) || sgn(M(2, 0)));

    if (sgn(M(0, 0)))
      i = 0;
    else if (sgn(M(1, 0)))
      i = 1;
    else  //  if (sgn(M(2, 0)))
      i = 2;

    if (i > 0)
      for (j = 0; j < 3; j++)
      {
        t       = M(0, j);
        M(0, j) = M(i, j);
        M(i, j) = t;
      }

    for (i = 1; i < 3; i++)
    {
      t = M(i, 0) / M(0, 0);

      for (j = 0; j < 3; j++)
        M(i, j) -= t * M(0, j);
    }

    if (sgn(M(1, 1)))
      v[1] = - M(1, 2) / M(1, 1);
    else
      v[1] = - M(2, 2) / M(2, 1);

    v[0] = - (M(0, 2) + M(0, 1) * v[1]) / M(0, 0);

    //  Set trim_curves.

    for (i = 0; i < num_trim_curves; i++)
    {
      if (i == 0)  //  trim_curves[0]: t = 0 from (0, 0) to (1, 0)
      {
        d[0] = 0;
        d[1] = 1;
        P    = K_RATPOLY(2, d);

        p[0] = 0;
        p[1] = 1;
        P.get_coeff(p) = 1;

        P.reduce_deg();

        start = new K_POINT2D(0, 0);
        end   = new K_POINT2D(1, 0);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }
      else if (i == 1)
      //  trim_curves[1]: - 1 + s + (1 - v[0]) / v[1] t = 0 from (1, 0) to v
      {
        d[0] = d[1] = 1;
        P    = K_RATPOLY(2, d);

        p[0] = p[1] = 0;
        P.get_coeff(p) = - 1;
        p[0] = 1;
        p[1] = 0;
        P.get_coeff(p) = 1;
        p[0] = 0;
        p[1] = 1;
        P.get_coeff(p) = (1 - v[0]) / v[1];

        P.reduce_deg();

        start = new K_POINT2D(1, 0);
        end   = new K_POINT2D(v[0], v[1]);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }
      else if (i == 2)
      //  trim_curves[2]: - 1 + (1 - v[1]) / v[0] s + t from v to (0, 1)
      {
        d[0] = d[1] = 1;
        P    = K_RATPOLY(2, d);

        p[0] = p[1] = 0;
        P.get_coeff(p) = - 1;
        p[0] = 1;
        p[1] = 0;
        P.get_coeff(p) = (1 - v[1]) / v[0];
        p[0] = 0;
        p[1] = 1;
        P.get_coeff(p) = 1;

        P.reduce_deg();

        start = new K_POINT2D(v[0], v[1]);
        end   = new K_POINT2D(0, 1);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }
      else  //  if (i == 3)  trim_curves[3]: s = 0 from (0, 1) to (0, 0)
      {
        d[0] = 1;
        d[1] = 0;
        P    = K_RATPOLY(2, d);

        p[0] = 1;
        p[1] = 0;
        P.get_coeff(p) = 1;

        P.reduce_deg();

        start = new K_POINT2D(0, 1);
        end   = new K_POINT2D(0, 0);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }

      trim_curves[i] = new K_CURVE(P, segment, 1);
      trim_curves[i]->ref_count++;

      if (!--segment[0]->ref_count)
        delete segment[0];

      delete [] segment;
    }
  }
  else  //  if (!is_plane)
  {
    //  A bilinear patch. Set trim_curves to be t = 0, s = 1, t = 1 and s = 0.

    for (i = 0; i < num_trim_curves; i++)
    {
      if (i == 0)  //  trim_curves[0]: t = 0 from (0, 0) to (1, 0)
      {
        d[0] = 0;
        d[1] = 1;
        P    = K_RATPOLY(2, d);

        p[0] = 0;
        p[1] = 1;
        P.get_coeff(p) = 1;

        P.reduce_deg();

        start = new K_POINT2D(0, 0);
        end   = new K_POINT2D(1, 0);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }
      else if (i == 1)  //  trim_curves[1]: - 1 + s = 0 from (1, 0) to (1, 1)
      {
        d[0] = 1;
        d[1] = 0;
        P    = K_RATPOLY(2, d);

        p[0] = p[1] = 0;
        P.get_coeff(p) = - 1;
        p[0] = 1;
        p[1] = 0;
        P.get_coeff(p) = 1;

        P.reduce_deg();

        start = new K_POINT2D(1, 0);
        end   = new K_POINT2D(1, 1);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }
      else if (i == 2)  //  trim_curves[2]: - 1 + t from (1, 1) to (0, 1)
      {
        d[0] = 0;
        d[1] = 1;
        P    = K_RATPOLY(2, d);

        p[0] = p[1] = 0;
        P.get_coeff(p) = - 1;
        p[0] = 0;
        p[1] = 1;
        P.get_coeff(p) = 1;

        P.reduce_deg();

        start = new K_POINT2D(1, 1);
        end   = new K_POINT2D(0, 1);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }
      else  //  if (i == 3)  trim_curves[3]: s = 0 from (0, 1) to (0, 0)
      {
        d[0] = 1;
        d[1] = 0;
        P    = K_RATPOLY(2, d);

        p[0] = 1;
        p[1] = 0;
        P.get_coeff(p) = 1;

        P.reduce_deg();

        start = new K_POINT2D(0, 1);
        end   = new K_POINT2D(0, 0);

        segment    = new K_SEGMENT* [1];
        segment[0] = new K_SEGMENT(start, end);
        segment[0]->ref_count++;
      }

      trim_curves[i] = new K_CURVE(P, segment, 1);
      trim_curves[i]->ref_count++;

      if (!--segment[0]->ref_count)
        delete segment[0];

      delete [] segment;
    }
  }

  patch = new K_PATCH(surf, trim_curves, num_trim_curves);

  for (i = 0; i < num_trim_curves; i++)
    if (!--trim_curves[i]->ref_count)
      delete trim_curves[i];

  delete [] trim_curves;

  return 0;
}

int planarize(const bigrational_vector inp_pts[], bigrational_vector pts[],
              const unsigned long      num_pts)
{
  assert(num_pts == 8);

  const unsigned long num_patches = 6;

  unsigned long      i, j;
  bigrational*       planes[6];
  bigrational_vector int_pts[8];
  bigrational        d, max_dev, t, tol;
  int                planarized;

  //
  //        7-------------6
  //       /|            /|
  //      / |           / |
  //     /  |          /  |
  //    4-------------5   |
  //    |   |         |   |
  //    |   |         |   |
  //    |   |         |   |
  //    |   0---------|---1
  //    |  /          |  /
  //    | /           | /
  //    |/            |/
  //    3-------------2
  //
  //  1.  Compute the implicit representations for
  //        plane[0], plane[1], ..., plane[5].

  assert(get_plane_coeffs(inp_pts[0], inp_pts[1], inp_pts[3], planes[0]));
  assert(get_plane_coeffs(inp_pts[0], inp_pts[1], inp_pts[7], planes[1]));
  assert(get_plane_coeffs(inp_pts[5], inp_pts[2], inp_pts[6], planes[2]));
  assert(get_plane_coeffs(inp_pts[5], inp_pts[4], inp_pts[2], planes[3]));
  assert(get_plane_coeffs(inp_pts[0], inp_pts[3], inp_pts[7], planes[4]));
  assert(get_plane_coeffs(inp_pts[5], inp_pts[4], inp_pts[6], planes[5]));

  //  2.  Compute 8 corners as intersections
  //        int_pts[0], int_pts[1], ..., int_pts[7]
  //      of 3 planes.

  intersect_planes(planes[0], planes[1], planes[4], int_pts[0]);
  intersect_planes(planes[0], planes[1], planes[2], int_pts[1]);
  intersect_planes(planes[0], planes[2], planes[3], int_pts[2]);
  intersect_planes(planes[0], planes[3], planes[4], int_pts[3]);
  intersect_planes(planes[3], planes[4], planes[5], int_pts[4]);
  intersect_planes(planes[2], planes[3], planes[5], int_pts[5]);
  intersect_planes(planes[1], planes[2], planes[5], int_pts[6]);
  intersect_planes(planes[1], planes[4], planes[5], int_pts[7]);

//  for (i = 0; i < num_pts; i++)
//  {
//    cerr << " genbox: planarize: inp_pts[" << i << "] =" << inp_pts[i] << endl << flush;
//    cerr << " genbox: planarize: int_pts[" << i << "] =" << int_pts[i] << endl << flush;
//  }
//  cerr << endl << flush;

  //
  //        7-------------6
  //       /|            /|
  //      / |           / |
  //     /  |          /  |
  //    4-------------5   |
  //    |   |         |   |
  //    |   |         |   |
  //    |   |         |   |
  //    |   0---------|---1
  //    |  /          |  /
  //    | /           | /
  //    |/            |/
  //    3-------------2
  //
  //  3.  See if its_pts[i] deviates from inp_pts[i] within a tolerance.

  max_dev = 0;

  for (i = 0; i < num_pts; i++)
    for (j = 0; j < 3; j++)
      if ((d = abs(int_pts[i][j] - inp_pts[i][j])) > max_dev)
        max_dev = d;

  tol = in_prod(inp_pts[7] - inp_pts[0], inp_pts[7] - inp_pts[0]) / 10000;

//  tol = in_prod(inp_pts[7] - inp_pts[0], inp_pts[7] - inp_pts[0]);
//
//  if (tol > (t = in_prod(inp_pts[3] - inp_pts[0], inp_pts[3] - inp_pts[0])))
//    tol = t;
//
//  if (tol > (t = in_prod(inp_pts[1] - inp_pts[0], inp_pts[1] - inp_pts[0])))
//    tol = t;
//
//  tol = tol / 16384;

  if (max_dev < tol)
  {
    for (i = 0; i < num_pts; i++)
      pts[i] = int_pts[i];

    planarized = 1;
  }
  else  //  if (max_dev >= tol)
  {
    for (i = 0; i < num_pts; i++)
      pts[i] = inp_pts[i];

    planarized = 0;
  }

//  cerr << " genbox: planarize: max_dev = " << max_dev << ", tol = " << tol << endl << flush;
//  cerr << endl << flush;

//  for (i = 0; i < num_pts; i++)
//    cerr << " genbox: planarize: pts[" << i << "] =" << pts[i] << endl << flush;
//  cerr << endl << flush;

  for (i = 0; i < num_patches; i++)
    delete [] planes[i];

  return planarized;
}

//  int perturb_box(const bigrational_vector pts[],
//                  bigrational_vector       pertrubed_pts[],
//                  const unsigned long      num_pts,
//                  const bigrational&       factor)
//    perturbs a box whose corners are pts to
//             the box whose corners are perturbed_pts
//      outward by "factor" if factor > 0 and
//      inward  by |factor| if factor < 0.

int perturb_box(const bigrational_vector pts[],
                bigrational_vector       pertrubed_pts[],
                const unsigned long      num_pts,
                const bigrational&       factor)
{
  assert(num_pts == 8);
  assert(abs(factor) < 1);

  unsigned long       i, j;
  bigrational_vector  center(3);
  bigrational_vector* from_center;
  bigrational         squared_dist_proto, squared_dist;

  for (i = 0; i < 3; i++)
    center[i] = 0;

  for (i = 0; i < num_pts; i++)
    center += pts[i];

  center = scalar_mul(bigrational(1, num_pts), center);

  from_center = new bigrational_vector [num_pts];  //  num_pts == 8

  from_center[0] = pts[0] - center;
//  squared_dist   = in_prod(from_center[0], from_center[0]);
//  j              = 0;

  for (i = 1; i < num_pts; i++)
  {
    from_center[i]     = pts[i] - center;
//    squared_dist_proto = in_prod(from_center[i], from_center[i]);
//
//    if (squared_dist < squared_dist_proto)
//    {
//      squared_dist = squared_dist_proto;
//      j            = i;
//    }
  }

//  center = center + scalar_mul(factor * factor * squared_dist, from_center[j]);

  for (i = 0; i < num_pts; i++)
    pertrubed_pts[i] = center + scalar_mul(1 + factor, from_center[i]);

  delete [] from_center;  //  num_pts == 8 => from_center != 0

  return 0;
}

//  K_SOLID gen_box(const bigrational_vector pts[],
//                  const unsigned long num_pts)
//    returns a box whose corners are pts.

K_SOLID gen_box(const bigrational_vector pts[], const unsigned long num_pts)
{
  assert(num_pts == 8);

  const unsigned long num_patches = 6;

  unsigned long      i;
  bigrational_vector four_pts[4];
  K_PATCH*           patches[6];
  K_SOLID            s;

  //
  //        7-------------6
  //       /|            /|
  //      / |           / |
  //     /  |          /  |
  //    4-------------5   |
  //    |   |         |   |
  //    |   |         |   |
  //    |   |         |   |
  //    |   0---------|---1
  //    |  /          |  /
  //    | /           | /
  //    |/            |/
  //    3-------------2
  //
  //  1.  Set pathces.

  //  patches[0]: (0 1 2 3)

  four_pts[0] = pts[0];
  four_pts[1] = pts[1];
  four_pts[2] = pts[2];
  four_pts[3] = pts[3];
  get_patch4(four_pts, patches[0]);

  //  patches[1]: (1 0 7 6)

  four_pts[0] = pts[1];
  four_pts[1] = pts[0];
  four_pts[2] = pts[7];
  four_pts[3] = pts[6];
  get_patch4(four_pts, patches[1]);

  //  patches[2]: (2 1 6 5)

  four_pts[0] = pts[2];
  four_pts[1] = pts[1];
  four_pts[2] = pts[6];
  four_pts[3] = pts[5];
  get_patch4(four_pts, patches[2]);

  //  patches[3]: (3 2 5 4)

  four_pts[0] = pts[3];
  four_pts[1] = pts[2];
  four_pts[2] = pts[5];
  four_pts[3] = pts[4];
  get_patch4(four_pts, patches[3]);

  //  patches[4]: (0 3 4 7)

  four_pts[0] = pts[0];
  four_pts[1] = pts[3];
  four_pts[2] = pts[4];
  four_pts[3] = pts[7];
  get_patch4(four_pts, patches[4]);

  //  patches[5]: (4 5 6 7)

  four_pts[0] = pts[4];
  four_pts[1] = pts[5];
  four_pts[2] = pts[6];
  four_pts[3] = pts[7];
  get_patch4(four_pts, patches[5]);

  //  2.  For each patch, set its adj_surfs and adj_patches.

  //  patches[0]:

  patches[0]->adj_surfs[0] = patches[1]->surf;
  patches[0]->adj_surfs[0]->ref_count++;
  patches[0]->adj_surfs[1] = patches[2]->surf;
  patches[0]->adj_surfs[1]->ref_count++;
  patches[0]->adj_surfs[2] = patches[3]->surf;
  patches[0]->adj_surfs[2]->ref_count++;
  patches[0]->adj_surfs[3] = patches[4]->surf;
  patches[0]->adj_surfs[3]->ref_count++;

  patches[0]->adj_patches[0] = patches[1];
  patches[0]->adj_patches[0]->ref_count++;
  patches[0]->adj_patches[1] = patches[2];
  patches[0]->adj_patches[1]->ref_count++;
  patches[0]->adj_patches[2] = patches[3];
  patches[0]->adj_patches[2]->ref_count++;
  patches[0]->adj_patches[3] = patches[4];
  patches[0]->adj_patches[3]->ref_count++;

  //  patches[1]:

  patches[1]->adj_surfs[0] = patches[0]->surf;
  patches[1]->adj_surfs[0]->ref_count++;
  patches[1]->adj_surfs[1] = patches[4]->surf;
  patches[1]->adj_surfs[1]->ref_count++;
  patches[1]->adj_surfs[2] = patches[5]->surf;
  patches[1]->adj_surfs[2]->ref_count++;
  patches[1]->adj_surfs[3] = patches[2]->surf;
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
  patches[2]->adj_surfs[1] = patches[1]->surf;
  patches[2]->adj_surfs[1]->ref_count++;
  patches[2]->adj_surfs[2] = patches[5]->surf;
  patches[2]->adj_surfs[2]->ref_count++;
  patches[2]->adj_surfs[3] = patches[3]->surf;
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
  patches[3]->adj_surfs[1] = patches[2]->surf;
  patches[3]->adj_surfs[1]->ref_count++;
  patches[3]->adj_surfs[2] = patches[5]->surf;
  patches[3]->adj_surfs[2]->ref_count++;
  patches[3]->adj_surfs[3] = patches[4]->surf;
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
  patches[4]->adj_surfs[1] = patches[3]->surf;
  patches[4]->adj_surfs[1]->ref_count++;
  patches[4]->adj_surfs[2] = patches[5]->surf;
  patches[4]->adj_surfs[2]->ref_count++;
  patches[4]->adj_surfs[3] = patches[1]->surf;
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

  patches[5]->adj_surfs[0] = patches[3]->surf;
  patches[5]->adj_surfs[0]->ref_count++;
  patches[5]->adj_surfs[1] = patches[2]->surf;
  patches[5]->adj_surfs[1]->ref_count++;
  patches[5]->adj_surfs[2] = patches[1]->surf;
  patches[5]->adj_surfs[2]->ref_count++;
  patches[5]->adj_surfs[3] = patches[4]->surf;
  patches[5]->adj_surfs[3]->ref_count++;

  patches[5]->adj_patches[0] = patches[3];
  patches[5]->adj_patches[0]->ref_count++;
  patches[5]->adj_patches[1] = patches[2];
  patches[5]->adj_patches[1]->ref_count++;
  patches[5]->adj_patches[2] = patches[1];
  patches[5]->adj_patches[2]->ref_count++;
  patches[5]->adj_patches[3] = patches[4];
  patches[5]->adj_patches[3]->ref_count++;

  //  3.  Associate edges.

  patches[0]->trim_curves[0]->assoc(patches[1]->trim_curves[0], - 1);
  patches[0]->trim_curves[1]->assoc(patches[2]->trim_curves[0], - 1);
  patches[0]->trim_curves[2]->assoc(patches[3]->trim_curves[0], - 1);
  patches[0]->trim_curves[3]->assoc(patches[4]->trim_curves[0], - 1);

  patches[5]->trim_curves[0]->assoc(patches[3]->trim_curves[2], - 1);
  patches[5]->trim_curves[1]->assoc(patches[2]->trim_curves[2], - 1);
  patches[5]->trim_curves[2]->assoc(patches[1]->trim_curves[2], - 1);
  patches[5]->trim_curves[3]->assoc(patches[4]->trim_curves[2], - 1);

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

K_SOLID read_box(istream& in_fs, const bigrational& perturb_factor)
{
  unsigned long       i, j;
//  unsigned long       type;
  unsigned long       num_pts;
  bigrational_vector* inp_pts;
  bigrational_vector* pts;
  bigrational_vector* perturbed_pts;
  K_SOLID             s;

//  in_fs >> type;
//  cerr << " genbox: read_box: type = " << type << endl << flush;
//  assert(type == 1);
  in_fs >> num_pts;
//  cerr << " genbox: read_box: num_pts = " << num_pts << endl << flush;
  cerr << endl << flush;

  //  1.  Read inp_pts from is_fs.

  if (num_pts > 0)
  {
    inp_pts = new bigrational_vector [num_pts];  //  num_pts > 0

    for (i = 0; i < num_pts; i++)
    {
      inp_pts[i] = bigrational_vector(3);

      cerr << " genbox: read_box: inp_pts[" << i << "] = ( ";
      for (j = 0; j < 3; j++)
      {
        in_fs >> inp_pts[i][j];

        cerr << inp_pts[i][j];
        if (j < 2)
          cerr << ", ";
      }
      cerr << " )" << endl << flush;
    }
    cerr << endl << flush;
  }
  else  //  if (num_pts == 0)
    inp_pts = 0;

  //  2.  Planarize inp_pts to pts.

  if (num_pts > 0)
  {
    pts = new bigrational_vector [num_pts];  //  num_pts > 0

    planarize(inp_pts, pts, num_pts);

    for (i = 0; i < num_pts; i++)
      cerr << " genbox: read_box: pts[" << i << "] =" << pts[i] << endl << flush;
    cerr << endl << flush;

    delete [] inp_pts;  //  num_pts > 0 => inp_pts != 0
  }
  else  //  if (num_pts == 0)
    pts = 0;

  //  3.  Perturb pts if necessarily and generate a box.

  if (!sgn(perturb_factor))
  {
    if (num_pts > 0)
      s = gen_box(pts, num_pts);

    if (num_pts > 0)
      delete [] pts;  //  num_pts > 0 => pts != 0
  }
  else  //  if (sgn(perturb_factor))
  {
    if (num_pts > 0)
    {
      perturbed_pts = new bigrational_vector [num_pts];  //  num_pts > 0

      perturb_box(pts, perturbed_pts, num_pts, perturb_factor);

      for (i = 0; i < num_pts; i++)
        cerr << " genbox: read_box: perturbed_pts[" << i << "] =" << perturbed_pts[i] << endl << flush;
      cerr << endl << flush;

      delete [] pts;  //  num_pts > 0 => pts != 0
    }
    else  //  if (num_pts == 0)
      perturbed_pts = 0;

    if (num_pts > 0)
      s = gen_box(perturbed_pts, num_pts);

    if (num_pts > 0)
      delete [] perturbed_pts;  //  num_pts > 0 => perturbed_pts != 0
  }

  return s;
}

K_SOLID read_BRLCAD_box(istream& in_fs, const bigrational& perturb_factor)
{
  unsigned long       i, j;
//  unsigned long       type;
  unsigned long       num_pts;
  bigrational_vector* inp_pts;
  bigrational_vector* pts;
  float               f;
  bigrational_vector* n;
  bigrational         o;
  bigrational_vector  p;
  bigrational_vector* perturbed_pts;
  K_SOLID             s;

//  in_fs >> type;
//  assert(type == 11);
//  cerr << " genbox: read_BRLCAD_box: type = " << type << endl << flush;
  num_pts = 8;
//  cerr << " genbox: read_BRLCAD_box: num_pts = " << num_pts << endl << flush;
  assert(num_pts == 8);
  cerr << endl << flush;

  //  1.  Read inp_pts from is_fs.

  if (num_pts > 0)
  {
    inp_pts = new bigrational_vector [num_pts];  //  num_pts == 8

    for (i = 0; i < num_pts; i++)
      inp_pts[i] = bigrational_vector(3);

    cerr << " genbox: read_BRLCAD_box: inp_pts[0] = ( ";
    for (j = 0; j < 3; j++)
    {
      in_fs >> f;
      inp_pts[0][j] = as_bigrational(f);

      cerr << inp_pts[0][j];
      if (j < 2)
        cerr << ", ";
    }
    cerr << " )" << endl << flush;

    for (i = 1; i < num_pts; i++)
    {
      cerr << " genbox: read_BRLCAD_box: inp_pts[" << i << "] = ( ";
      for (j = 0; j < 3; j++)
      {
        in_fs >> f;
        inp_pts[i][j] = as_bigrational(f) + inp_pts[0][j];

        cerr << inp_pts[i][j];
        if (j < 2)
          cerr << ", ";
      }
      cerr << " )" << endl << flush;
    }
    cerr << endl << flush;
  }
  else  //  if (num_pts == 0)
    inp_pts = 0;

  //  1.5.  Sort inp_pts clockwise looking from outside.

  if (num_pts > 0)
  {
    n = new bigrational_vector [3];

    n[0] = inp_pts[1] - inp_pts[0];
    n[1] = inp_pts[3] - inp_pts[0];
    n[2] = inp_pts[4] - inp_pts[0];

    o = (n[0][1] * n[1][2] - n[0][2] * n[1][1]) * n[2][0]
      + (n[0][2] * n[1][0] - n[0][0] * n[1][2]) * n[2][1]
      + (n[0][0] * n[1][1] - n[0][1] * n[1][0]) * n[2][2];

    if (sgn(o) < 0)
    {
      p          = inp_pts[1];
      inp_pts[1] = inp_pts[3];
      inp_pts[3] = p;

      p          = inp_pts[5];
      inp_pts[5] = inp_pts[7];
      inp_pts[7] = p;
    }

    p          = inp_pts[0];
    inp_pts[0] = inp_pts[7];
    inp_pts[7] = inp_pts[4];
    inp_pts[4] = p;

    p          = inp_pts[1];
    inp_pts[1] = inp_pts[6];
    inp_pts[6] = inp_pts[5];
    inp_pts[5] = p;

    delete [] n;  //  n != 0

    for (i = 0; i < num_pts; i++)
      cerr << " genbox: read_BRLCAD_box: inp_pts[" << i << "] = " << inp_pts[i] << endl << flush;
    cerr << endl << flush;
  }

  //  2.  Planarize inp_pts to pts.

  if (num_pts > 0)
  {
    pts = new bigrational_vector [num_pts];  //  num_pts == 8

    planarize(inp_pts, pts, num_pts);

    for (i = 0; i < num_pts; i++)
      cerr << " genbox: read_BRLCAD_box: pts[" << i << "] = " << pts[i] << endl << flush;
    cerr << endl << flush;

    delete [] inp_pts;  //  num_pts > 0 => inp_pts != 0
  }
  else  //  if (num_pts == 0)
    pts = 0;

  //  3.  Perturb pts if necessarily and generate a box.

  if (!sgn(perturb_factor))
  {
    if (num_pts > 0)
      s = gen_box(pts, num_pts);

    if (num_pts > 0)
      delete [] pts;  //  num_pts > 0 => pts != 0
  }
  else  //  if (sgn(perturb_factor))
  {
    if (num_pts > 0)
    {
      perturbed_pts = new bigrational_vector [num_pts];  //  num_pts > 0

      perturb_box(pts, perturbed_pts, num_pts, perturb_factor);

      for (i = 0; i < num_pts; i++)
        cerr << " genbox: read_BRLCAD_box: perturbed_pts[" << i << "] = " << perturbed_pts[i] << endl << flush;
      cerr << endl << flush;

      delete [] pts;  //  num_pts > 0 => pts != 0
    }
    else  //  if (num_pts == 0)
      perturbed_pts = 0;

    if (num_pts > 0)
      s = gen_box(perturbed_pts, num_pts);

    if (num_pts > 0)
      delete [] perturbed_pts;
  }

  return s;
}

