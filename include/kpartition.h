//  file:    kpartition.h
//  update:  01/15/03

#ifndef _KPARTITION_H
#define _KPARTITION_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <kpatch.h>

using namespace std;

class K_GRAPH;

class K_PARTITION
{
friend class K_PATCH;
friend class K_SOLID;
  
  unsigned long ID;
  
  K_PATCH* from;  //  patch from where *this comes
  
  bool is_head;  //  true if *this is the head of surf
  
  unsigned long num_trim_curves;
  K_CURVE**     trim_curves;
  long*         adj_curves;
  K_SURF**      adj_surfs;
  K_PARTITION** adj_partitions;
  int*          rev_curves;
  
  unsigned long ref_count;
  
  //  stream
  
  ostream& output(ostream&) const;
  
public:
  
  K_PARTITION();
  K_PARTITION(K_PATCH* const);
  
  K_PARTITION(const K_PARTITION&);
  K_PARTITION& operator =(const K_PARTITION&);
  
  ~K_PARTITION();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_PARTITION&);
  
  //  other functions
  
  K_POINT2D get_pt_in() const;
  
  friend unsigned long gen_partitions(K_PATCH* const, K_PARTITION**&);
  
  friend unsigned long gen_adjacency(K_PARTITION** const, const unsigned long,
                                     K_GRAPH*&,
                                     long*&, K_GRAPH*&);
  
  friend unsigned long select_relevant_partitions(K_PARTITION** const,
                                                  const unsigned long,
                                                  K_GRAPH* const,
                                                  const long* const,
                                                  const int* const,
                                                  const int,
                                                  K_PARTITION**&,
                                                  K_GRAPH*&);
  
//  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
//                               K_PARTITION**, const unsigned long);
  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
                               K_PARTITION**, const unsigned long,
                               const char);
};

#endif

