#ifndef _PASCAL_H
#define _PASCAL_H

#include <cassert>
#include <cstdlib>
#include <iostream>

using namespace std;

class Pascal
{
  unsigned long num_row;
  unsigned long num_col;
  long*         tab;
  
public:
  
  Pascal();
  Pascal(const unsigned long, const unsigned long);
  ~Pascal();
  
  int  resize(const unsigned long, const unsigned long);
  long get_Pascal(const unsigned long, const unsigned long);
};

#endif

