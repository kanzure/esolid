//  file:    config.h
//  update:  10/14/02

#ifndef _CONFIG_H
#define _CONFIG_H

#include <bigrational.h>

//  Preprocessor Symbols depending "only" on the system.
////    In makefile, put
////                      CFLAGS += -D_Linux_i386_
////    or some such.
//#define _Linux_i386_
//#define _SunOS_
//#define _IRIX_

#define DIM 3
#define MAX_NUM_GOOD_FP_BITS 30
#define IN  0
#define OUT 1

extern bigrational shrink_step;
extern bigrational init_tol;
extern bigrational epsilon;
extern bigrational s_epsilon;
extern bigrational t_epsilon;

#endif

