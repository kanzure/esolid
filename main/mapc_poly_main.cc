//  file:   mapc_poly_main.cc
//  update: 04/15/04

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <kratpoly.h>
#include <kfloatpoly.h>

using namespace std;

//int main(const int argc, const char* argv[])
int main()
{
//  ofstream of(argv[1]);
  
  unsigned long i, j;
  i = 2;
  j = 1;
  
  if (i > j)
    cerr << " i > j " << endl << flush;
  
//  long        xp[] = { 1, 1 };
//  bigrational xc[] = { 1, 1, 1, 1 };
  long        xp[] = { 1, 2 };
  bigrational xc[] = { 1, 1, 1, 0, 1, 1 };
  K_RATPOLY   X(2, xp, 6, xc);
  cerr << " X = " << endl << flush;
  cerr << X << endl << flush;
  
//  long        yp[] = { 3 };
//  bigrational yc[] = { bigrational(1, 2), 0, 0, bigrational(1, 2) };
//  bigrational yc[] = { bigrational(1, 2), 0, 0, 0 };  
  long        yp[] = { 4 };
  bigrational yc[] = { bigrational(1, 2), - bigrational(1, 2), - bigrational(3, 2), bigrational(5, 2), - 1 };
  K_RATPOLY   Y(1, yp, 5, yc);
  cerr << " Y = " << endl << flush;
  cerr << Y << endl << flush;
  
  long        zp[] = { 1 };
  bigrational zc[] = { bigrational(1, 2), bigrational(1, 3) };
//  long        zp[] = { 0 };
//  bigrational zc[] = { bigrational(1, 2) };  
  K_RATPOLY   Z(1, zp, 2, zc);
  cerr << " Z = " << endl << flush;
  cerr << Z << endl << flush;
  
//  K_RATPOLY W = X.subst_expr(0, Y, Z);
  K_RATPOLY W = X.subst_expr(0, Z);  
  cerr << " W = " << endl << flush;
  cerr << W << endl << flush;
  
//  K_RATPOLY Q, R;
//  Q = div(Y, Z, R);
  
//  cerr << " Q = " << endl << flush;
//  cerr << Q << endl << flush;
//  cerr << " R = " << endl << flush;
//  cerr << R << endl << flush;

//  K_RATPOLY G = gcd(Z, Y);
//  cerr << " G = " << endl << flush;
//  cerr << G << endl << flush;
  
  bigrational v   = - 1;
  bigrational V[] = { - 1 };
  cerr << " Y(" << v << ") = " << Y.evaluate(v) << endl << flush;
  cerr << " Y(" << V[0] << ") = " << Y.evaluate(V) << endl << flush;
  cerr << " sgn(Y(" << V[0] << ")) = " << Y.sgn_at(V) << endl << flush;

  Y.set_Sturm_seq();
  cerr << Y.num_Sturm_seq_perm(- 1) << endl << flush;
  cerr << Y.num_Sturm_seq_perm(0) << endl << flush;
  cerr << Y.num_Sturm_seq_perm(1) << endl << flush;
  
//  long        fp[]   = { 3 };
//  double      fc[] = { 1.0, 0.0, - 1.0, 0.0 };
//  K_FLOATPOLY F(1, fp, fc);
//  cerr << " F = " << endl << flush;
//  cerr << F << endl << flush;
//  
//  unsigned long l_FR;
//  double*       FR;
//  l_FR = F.gen_fp_roots(-2.0, 2.0, &FR);
//  cerr << " F has " << l_FR << " real roots " << endl << flush;
//  
//  for (i = 0; i < l_FR; i++)
//    cerr << FR[i] << endl << flush;
//  
//  if (l_FR)
//    delete [] FR;
  
  return 0;
}
