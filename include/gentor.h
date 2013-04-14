#ifndef _GENTOR_H
#define _GENTOR_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <ksolid.h>

using namespace std;

K_SOLID read_tor(istream&, const bigrational& = 0);
K_SOLID read_BRLCAD_tor(istream&, const bigrational& = 0);

#endif

