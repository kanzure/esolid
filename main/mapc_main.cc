//  file:  msri_main.cc
//  update:  04/13/04

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <timer.h>

#include <config.h>
#include <kpoint1d.h>
#include <kpoint2d.h>

using namespace std;

unsigned long num_kratpoly_sgn_at       = 0;
unsigned long num_kratpoly_exact_sgn_at = 0;
unsigned long num_root1_isolate         = 0;
unsigned long num_root1_exact_isolate   = 0;
unsigned long num_kpoint1d_get_pts      = 0;
unsigned long num_kpoint2d_get_pts      = 0;
unsigned long num_ksolid_boolean        = 0;

int main(const int argc, const char* argv[])
{
  ifstream in_file;
  
  long          i;
  K_RATPOLY     P, Q;
  bigrational   l_s, h_s, l_t, h_t;
  unsigned long num_X;
  K_POINT2D**   X;
  
  long lap;
  
  in_file.open(argv[1]);
  P = read_poly(in_file);
  in_file.close();
  
  in_file.open(argv[2]);
  Q = read_poly(in_file);
  in_file.close();
  
  cerr << " P = " << endl << P << endl << flush;
  cerr << " Q = " << endl << Q << endl << flush;
  
  CLOCK_START();  
  num_X = get_all_pts(P, Q, X, 0);
  CLOCK_STOP(&lap);
  
  cerr << " num_X = " << num_X << endl << flush;
  for (i = 0; i < num_X; i++)
    cerr << " X[" << i << "] = " << *X[i] << endl << flush;
  cerr << " mapc_main: main: time_taken: " << lap << " micro-seconds. " << endl << flush;
  
  for (i = 0; i < num_X; i++)
    delete X[i];
  
  if (num_X > 0)
    delete [] X;
  
  return 0;
}

