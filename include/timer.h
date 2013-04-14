//  file :  timer.h
//  last update:  07/09/01

#ifndef __TIMER_H
#define __TIMER_H

#include <cstdlib>
#include <sys/time.h>

static struct timeval  TV;
static struct timezone TZ;
static long   s_utime, s_time, u_utime, u_time;

void CLOCK_START();
void CLOCK_STOP(long*);

#endif

