#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <timer.h>

#include <ksolid.h>
#include <genbox.h>
#include <gencyl.h>
#include <genell.h>
#include <gentor.h>

using namespace std;

#ifdef _EXPERIMENT
unsigned long num_kratpoly_sgn_at       = 0;
unsigned long num_kratpoly_exact_sgn_at = 0;
unsigned long num_root1_isolate         = 0;
unsigned long num_root1_exact_isolate   = 0;
unsigned long num_kpoint1d_get_pts      = 0;
unsigned long num_kpoint2d_get_pts      = 0;
unsigned long num_ksolid_boolean        = 0;
#endif

int main(int argc, char* argv[])
{
  assert(argc >= 4);

  ifstream    in_file;
  ofstream    out_file;
  char        op;
  long        lap;
  bigrational perturb_factor1, perturb_factor2;
  bigrational perturb_factor;
  K_SOLID     s1, s2, s3;

  static const bigrational perturb_factor_proto = 0;
//  static const bigrational perturb_factor_proto = bigrational(1, 512);
//  static const bigrational perturb_factor_proto = bigrational(1, 1024);

  op = argv[1][1];
  assert(op == 'D' || op == 'I' || op == 'U' || op == 'O' || op == 'M');

  CLOCK_START();

  if (op == 'D' || op == 'I' || op == 'U')
  {
    assert(argc == 5);

    if (op == 'D')
    {
      perturb_factor1 = perturb_factor_proto;
      perturb_factor2 = - perturb_factor_proto;
    }
    else if (op == 'I')
    {
      perturb_factor1 = - perturb_factor_proto;
      perturb_factor2 = - perturb_factor_proto;
    }
    else  //  if (op == 'U')
    {
      perturb_factor1 = perturb_factor_proto;
      perturb_factor2 = perturb_factor_proto;
    }

    in_file.open(argv[2]);
    s1 = read_solid(in_file, perturb_factor1);
    in_file.close();

    in_file.open(argv[3]);
    s2 = read_solid(in_file, perturb_factor2);
    in_file.close();

    s3 = s1.boolean(s2, op);

    out_file.open(argv[4]);
    s3.Bezier_output(out_file);
    out_file.close();
  }
  else if (op == 'O')
  {
    assert(argc == 4);

    perturb_factor = perturb_factor_proto;

    in_file.open(argv[2]);
    s3 = read_solid(in_file, perturb_factor);
    in_file.close();

    out_file.open(argv[3]);
    s3.Bezier_output(out_file);
    out_file.close();
  }
  else if (op == 'M')
  {
    assert(argc == 4);

    perturb_factor = perturb_factor_proto;

    s3 = read_CSG(argv[2], argv[3], perturb_factor);

    out_file.open(argv[3]);
    s3.Bezier_output(out_file);
    out_file.close();
  }

  CLOCK_STOP(&lap);
  cerr << " esolid_main: main: lap = " << lap << " micro-seconds. " << endl << flush;
#ifdef _EXPERIMENT
  cerr << " esolid_main: main: num_kratpoly_sgn_at       = " << num_kratpoly_sgn_at << endl << "                    num_kratpoly_exact_sgn_at = " << num_kratpoly_exact_sgn_at << endl << flush;
  cerr << " esolid_main: main: num_root1_isolate         = " << num_root1_isolate << endl << "                    num_root1_exact_isolate   = " << num_root1_exact_isolate << endl << flush;
  cerr << " esolid_main: main: num_kpoint1d_get_pts      = " << num_kpoint1d_get_pts << endl << flush;
  cerr << " esolid_main: main: num_kpoint2d_get_pts      = " << num_kpoint2d_get_pts << endl << flush;
  cerr << " esolid_main: main: num_ksolid_boolean        = " << num_ksolid_boolean << endl << flush;
#endif

  return 0;
}

