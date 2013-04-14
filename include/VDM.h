#ifndef _VDM_H
#define _VDM_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <bigrational_vector.h>
#include <bigrational_matrix.h>

using namespace std;

bigrational_vector Solve_VDM_GE(const bigrational_vector& Val,
                                const bigrational_vector& RHS);

#endif

