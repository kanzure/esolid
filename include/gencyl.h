#ifndef _GENCYL_H
#define _GENCYL_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational_vector.h>

#include <ksolid.h>

using namespace std;

K_SOLID read_cyl(istream&, const bigrational& = 0);
K_SOLID read_BRLCAD_cyl(istream&, const bigrational& = 0);

#endif

