#ifndef _GENELL_H
#define _GENELL_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <ksolid.h>

using namespace std;

K_SOLID read_ell(istream&, const bigrational& = 0);
K_SOLID read_BRLCAD_ell(istream&, const bigrational& = 0);

#endif

