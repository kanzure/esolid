//  file:  pt1_main.cc
//  update:  09/25/02

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <fpconversion.h>
#include <kpoint1d.h>

using namespace std;

//int main(const int argc, const char* argv[])
int main()
{
//  ofstream of(argv[1]);
  
  long i, j;
  
  long        pxp[] = { 2 };
  bigrational pxc[] = { 1, 0, - bigrational(1, 4) };
//  long        pxp[] = { 1 };
//  bigrational pxc[] = { 1, - bigrational(1, 2) };
  K_RATPOLY   PX(1, pxp, 3, pxc);
  cerr << " PX = " << endl << flush;
  cerr << PX << endl << flush;
  
  long        pyp[] = { 2 };
  bigrational pyc[] = { 1, 0, - bigrational(1, 9) };
//  long        pyp[] = { 1 };
//  bigrational pyc[] = { 1, - bigrational(1, 3) };
  K_RATPOLY   PY(1, pyp, 3, pyc);
  cerr << " PY = " << endl << flush;
  cerr << PY << endl << flush;
  
  ROOT1 rx(PX, - bigrational(1, 4), 1);
  cerr << " rx = " << rx << endl << flush;
  
  ROOT1 ry(PY, - bigrational(1, 9), 1);
  cerr << " ry = " << ry << endl << flush;

  K_POINT1D x(rx);
  cerr << " x = " << x << endl << flush;
  
  K_POINT1D y(ry);
//  K_POINT1D y(bigrational(3));
  cerr << " y = " << y << endl << flush;
  
  K_POINT1D z = x / y;
  cerr << " z = " << z << endl << flush;

//  long        pwp[] = { 3 };
//  bigrational pwc[] = { 0, 1, 2, 0 };
//  K_RATPOLY   PW(1, pwp, 4, pwc);
//  cerr << " PW = " << endl << flush;
//  cerr << PW << endl << flush;
//  
//  PW.reduce_deg();
//  cerr << " PW = " << endl << flush;
//  cerr << PW << endl << flush;
//  
//  K_RATPOLY PT = PW.rm_x_factor();
//  cerr << " PT = " << endl << flush;
//  cerr << PT << endl << flush;
  
  return 0;
}

