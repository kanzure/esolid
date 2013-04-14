#define DEG 64

#include <cassert>
#include <cstdlib>
//#include <fstream>
#include <iostream>

#include <kratpoly.h>
#include <root1.h>

using namespace std;

int main(int argc, char* argv[])
{
//  ifstream in_file;
//  ofstream out_file;

  unsigned long i;
  long          t;

  K_RATPOLY P, Q;
  P = K_RATPOLY(1, 0, 1);
  for (i = 2; i <= DEG; i++)
  {
    Q = K_RATPOLY(1, 0, i);
    P = P * Q;
  }
  cerr << " P = " << endl << P << endl << flush;

  bigrational tol = bigrational(1, 65536);

  ROOT1*   R_0;
  unsigned long num_R_0;

  ROOT1 R_0_proto(P, - 1, DEG + 1);
//  num_R_0 = R_0_proto.isolate_roots(R_0, tol, 1);
  num_R_0 = R_0_proto.isolate_roots(R_0, tol);

  cerr << "num_roots: " << num_R_0 << endl << flush;
  cerr << "roots: " << endl << flush;
  for (i = 0; i < num_R_0; i++)
    cerr <<  R_0[i] << endl << flush;
  cerr << "----------------------------------------" << endl << flush;

  if (num_R_0 > 0)
    delete [] R_0;

  ROOT1*   R_1;
  unsigned long num_R_1;

  ROOT1 R_1_proto(P, - 1, DEG + 1);
//  num_R_1 = R_1_proto.isolate_roots(R_1, tol, 0);
  num_R_1 = R_1_proto.isolate_roots(R_1, tol);

  cerr << "num_roots: " << num_R_1 << endl << flush;
  cerr << "roots: " << endl << flush;
  for (i = 0; i < num_R_1; i++)
    cerr <<  R_1[i] << endl << flush;
  cerr << "----------------------------------------" << endl << flush;

  if (num_R_1 > 0)
    delete [] R_1;

  return 0;
}

