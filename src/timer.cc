//  file :  timer.cc
//  last update:  07/09/01

#include <cstdlib>

#include <timer.h>

void CLOCK_START()
{
  gettimeofday(&TV, &TZ);
  s_utime = TV.tv_usec;
  s_time  = TV.tv_sec;
}

void CLOCK_STOP(long* result)
{
  gettimeofday(&TV, &TZ);
  u_utime = TV.tv_usec;
  u_time  = TV.tv_sec;
  
  *result = (u_utime + u_time * 1000000) - (s_utime + s_time * 1000000);
//  *result = u_time - s_time;
}

