//  file:    config.cc
//  update:  10/14/02

#include <bigrational.h>

bigrational shrink_step = bigrational(1, 128);
bigrational init_tol    = bigrational(1, 1024);
bigrational epsilon     = bigrational(1, 8);
bigrational s_epsilon   = bigrational(1, 8);
bigrational t_epsilon   = bigrational(1, 16);

