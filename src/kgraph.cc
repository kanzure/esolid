#include <config.h>

#include <kgraph.h>

K_GRAPH :: K_GRAPH()
{
  num_vertices = 0;
  vertex_ID    = 0;
  edges        = 0;
}

K_GRAPH :: K_GRAPH(const unsigned long n)
{
  unsigned long v, w;
  
  if ((num_vertices = n) > 0)
  {
    vertex_ID = new unsigned long [num_vertices];
    edges     = new int* [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      edges[v] = new int [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      vertex_ID[v] = v;
    
    for (v = 0; v < num_vertices; v++)
      for (w = 0; w < num_vertices; w++)
        edges[v][w] = 0;
  }
  else  //  if (!num_vertices)
  {
    vertex_ID = 0;
    edges     = 0;
  }
}

K_GRAPH :: K_GRAPH(const unsigned long n, const unsigned long* const V)
{
  unsigned long v, w;
  
  if ((num_vertices = n) > 0)
  {
    vertex_ID = new unsigned long [num_vertices];
    edges     = new int* [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      edges[v] = new int [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      vertex_ID[v] = V[v];
    
    for (v = 0; v < num_vertices; v++)
      for (w = 0; w < num_vertices; w++)
        edges[v][w] = 0;
  }
  else  //  if (!num_vertices)
  {
    vertex_ID = 0;
    edges     = 0;
  }
}

K_GRAPH :: K_GRAPH(const K_GRAPH& G)
{
  unsigned long v, w;
  
  if ((num_vertices = G.num_vertices) > 0)
  {
    vertex_ID = new unsigned long [num_vertices];
    edges     = new int* [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      edges[v] = new int [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      vertex_ID[v] = G.vertex_ID[v];
    
    for (v = 0; v < num_vertices; v++)
      for (w = 0; w < num_vertices; w++)
        edges[v][w] = G.edges[v][w];
  }
  else  //  if (!num_vertices)
  {
    vertex_ID = 0;
    edges     = 0;
  }
}

K_GRAPH& K_GRAPH :: operator =(const K_GRAPH& G)
{
  unsigned long v, w;
  
  if (num_vertices > 0)
  {
    delete [] vertex_ID;
    
    for (v = 0; v < num_vertices; v++)
      delete [] edges[v];
    
    delete [] edges;
  }
  
  if ((num_vertices = G.num_vertices) > 0)
  {
    vertex_ID = new unsigned long [num_vertices];
    edges     = new int* [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      edges[v] = new int [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
      vertex_ID[v] = G.vertex_ID[v];
    
    for (v = 0; v < num_vertices; v++)
      for (w = 0; w < num_vertices; w++)
        edges[v][w] = G.edges[v][w];
  }
  else  //  if (!num_vertices)
  {
    vertex_ID = 0;
    edges     = 0;
  }
}

K_GRAPH :: ~K_GRAPH()
{
  unsigned long v;
  
  if (num_vertices > 0)
  {
    delete [] vertex_ID;
    
    for (v = 0; v < num_vertices; v++)
      delete [] edges[v];
    
    delete [] edges;
  }
}

ostream& K_GRAPH :: output(ostream& o) const
{
  unsigned long v, w;
  
  if (num_vertices > 0)
  {
    o << endl << " ------";
    
    for (v = 1; v < num_vertices; v++)
      o << "-------";
    
    o << "-- " << endl;
    
    for (v = 0; v < num_vertices; v++)
    {
      for (w = 0; w < num_vertices - 1; w++)
      {
        o.width(2);
        o << " | " << edges[v][w];
      }
      
      o.width(2);
      o << " | " << edges[v][num_vertices - 1] << " | " << endl;
    }
    
    o << endl << " ------";
    
    for (v = 1; v < num_vertices; v++)
      o << "-------";
    
    o << "-- " << endl << flush;
  }
  else
    o << " NULL " << flush;
  
  return o;
}

ostream& operator <<(ostream& o, const K_GRAPH& G)
{
  return G.output(o);
}

unsigned long K_GRAPH :: ID_to_index(const unsigned long x) const
{
  assert(num_vertices > 0);
  
  unsigned long v;
  
  for (v = 0; v < num_vertices && vertex_ID[v] != x; v++)
    ;
  
  return v;
}

int K_GRAPH :: add_edge(const unsigned long x, const unsigned long y,
                        const int ps)
{
  unsigned long v, w;
  
  v = ID_to_index(x);
  w = ID_to_index(y);
  
  assert(v < num_vertices && w < num_vertices);
  
  if (!ps)
    edges[v][w] = edges[w][v] = 1;
  else  //  if (ps)
    edges[v][w] = edges[w][v] = - 1;
  
  return 0;
}

int K_GRAPH :: del_edge(const unsigned long x, const unsigned long y)
{
  unsigned long v, w;
  
  v = ID_to_index(x);
  w = ID_to_index(y);
  
  assert(v < num_vertices && w < num_vertices);
  edges[v][w] = edges[w][v] = 0;
  
  return 0;
}

K_GRAPH K_GRAPH :: gen_subgraph(const unsigned long n,
                                const unsigned long* const V)
{
  unsigned long v, w, v_o, w_o;
  K_GRAPH       S = K_GRAPH(n, V);
  
  for (v = 0; v < n; v++)
  {
    assert((v_o = ID_to_index(V[v])) < num_vertices);
    
    for (w = 0; w < n; w++)
    {
      assert((w_o = ID_to_index(V[w])) < num_vertices);
      
      if (edges[v_o][w_o])
        S.edges[v][w] = S.edges[w][v] = 1;
    }
  }
  
  return S;
}

int K_GRAPH :: DFS(const unsigned long v, int* const visited,
                   long* const C, const unsigned long c)
{
  unsigned long w;
  
  visited[v] = 1;
  C[v]       = c;
  
  for (w = 0; w < num_vertices; w++)
    if (edges[v][w] > 0 && !visited[w])
      DFS(w, visited, C, c);
  
  return 0;
}

unsigned long K_GRAPH :: get_connected_components(long*& C, K_GRAPH*& G)
{
  unsigned long v, w;
  int*          visited;
  unsigned long num_components;
  
  if (num_vertices > 0)
  {
    visited = new int [num_vertices];
    C       = new long [num_vertices];
    
    for (v = 0; v < num_vertices; v++)
    {
      visited[v] = 0;
      C[v]       = - 1;
    }
    
    v              = 0;
    num_components = 0;
    
    while (v < num_vertices)
    {
      DFS(v, visited, C, num_components);
      num_components++;
      
      for (++v; v < num_vertices && visited[v]; v++)
        ;
    }
    
    delete [] visited;
    
//    for (v = 0; v < num_vertices; v++)
//      cerr << " kgraph: get_connected_components: C[" << v << "] = " << C[v] << endl << flush;
    
    G = new K_GRAPH(num_components);
    
    for (v = 0; v < num_vertices; v++)
      for (w = v + 1; w < num_vertices; w++)
        if (edges[v][w] < 0)
          G->add_edge(C[v], C[w], 0);
  }
  else  //  if (!num_vertices)
  {
    C              = 0;
    G              = 0;
    num_components = 0;
  }
  
  return num_components;
}

int K_GRAPH :: DFS_for_PC(const unsigned long v, int* const visited,
                          int* const C, const int cc)
{
  assert(cc == IN || cc == OUT);
  
  unsigned long w;
  int           nc;
  
  visited[v] = 1;
  C[v]       = cc;
  
  if (cc == IN)
    nc = OUT;
  else  //  if (cc == OUT)
    nc = IN;
  
  for (w = 0; w < num_vertices; w++)
//    if (edges[v][w] > 0 && !visited[w])
//      DFS_for_PC(w, visited, C, nc);
    if (edges[v][w] > 0)
      if (!visited[w])
        DFS_for_PC(w, visited, C, nc);
      else  //  if w has been visited
        assert(C[w] == nc);
  
  return 0;
}

int K_GRAPH :: propagate_color(const unsigned long x, int*& C, const int c)
{
  assert(num_vertices > 0);
  
  unsigned long v, w;
  int*          visited;
  unsigned long num_components;
  
  v = ID_to_index(x);
  assert(v < num_vertices);
  
  visited = new int [num_vertices];
  C       = new int [num_vertices];
  
  for (w = 0; w < num_vertices; w++)
    visited[w] = 0;
  
  DFS_for_PC(v, visited, C, c);
  
  delete [] visited;  //  num_vertcies > 0
  
//  for (w = 0; w < num_vertices; w++)
//    cerr << " kgraph: propagate_color: C[" << w << "] = " << C[w] << endl << flush;
  
  return 0;
}

