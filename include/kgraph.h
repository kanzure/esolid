#ifndef _KGRAPH_H
#define _KGRAPH_H

#include <cassert>
#include <cstdlib>
#include <iostream>

using namespace std;

class K_PARTITION;

class K_GRAPH
{
  unsigned long  num_vertices;  //  number of vertices
  unsigned long* vertex_ID;     //  array of ID's of vertices
  int**          edges;         //  adjacency matrix
  
  unsigned long ref_count;
  
  //  stream
  
  ostream& output(ostream&) const;

  //  primitives
  
  unsigned long ID_to_index(const unsigned long) const;
  int           add_edge(const unsigned long, const unsigned long, const int);
  int           del_edge(const unsigned long, const unsigned long);
  
  int DFS(const unsigned long, int* const, long* const, const unsigned long);
  int DFS_for_PC(const unsigned long, int* const, int* const, const int);
  
public:
  
  K_GRAPH();
  K_GRAPH(const unsigned long n);
  K_GRAPH(const unsigned long n, const unsigned long* const V);
  
  K_GRAPH(const K_GRAPH&);
  K_GRAPH& operator =(const K_GRAPH&);
  
  ~K_GRAPH();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_GRAPH&);
  
  //  other functions
  
  K_GRAPH gen_subgraph(const unsigned long, const unsigned long* const);
  
  unsigned long get_connected_components(long*&, K_GRAPH*&);
  int           propagate_color(const unsigned long, int*&, const int);
  
  friend unsigned long gen_adjacency(K_PARTITION** const, const unsigned long,
                                     K_GRAPH*&,
                                     long*&, K_GRAPH*&);
};

#endif

