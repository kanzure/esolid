//  file:    ksurf.h
//  update:  11/18/02

#ifndef _KSURF_H
#define _KSURF_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <bigrational_vector.h>
#include <kratpoly.h>
#include <kpoint2d.h>
#include <kbox3d.h>

using namespace std;

class K_SOLID;

class K_SURF
{
friend class PT_SURF_ASSOC;  
friend class K_PATCH;
friend class K_PARTITION;
friend class K_SOLID;
  
  K_RATPOLY* Impl;
  int        Impl_ok;
  
  K_RATPOLY* X;
  K_RATPOLY* Y;
  K_RATPOLY* Z;
  K_RATPOLY* W;
  int        mon_ok;
  
  K_RATPOLY* X_Bern;
  K_RATPOLY* Y_Bern;
  K_RATPOLY* Z_Bern;
  K_RATPOLY* W_Bern;
  int        Bern_ok;
  
  unsigned long ref_count;
  
  //  stream
  
  ostream& output(ostream&) const;
  
public:
  
  //  constructors, assignment and destructor
  
  K_SURF();
  K_SURF(const K_RATPOLY&);
  K_SURF(K_RATPOLY* const,
         K_RATPOLY* const, K_RATPOLY* const, K_RATPOLY* const,
         K_RATPOLY* const);
  
  K_SURF(const K_SURF&);
  K_SURF& operator =(const K_SURF&);
  
  ~K_SURF();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_SURF&);
  
  int Bezier_output(ostream&,
                    const bigrational&, const bigrational&,
                    const bigrational&, const bigrational&) const;
  
  //  other functions
  
  //  int param_to_coord(const bigrational& s, const bigrational& t,
  //                     bigrational& x, bigrational& y, bigrational& z) const
  //    compute (x, y, z) =
  //              (X(s, t) / W(s, t), Y(s, t) / W(s, t), Z(s, t) / W(s, t)).
  
  int param_to_coord(const bigrational&, const bigrational&,
                     bigrational&, bigrational&, bigrational&) const;
  
  //  int param_to_coord(const bigrational_vector& p,
  //                     bigrational_vector& c) const
  //    compute
  //      c = (X(p) / W(p), Y(p) / W(p), Z(p) / W(p)).
  
  int param_to_coord(const bigrational_vector&, bigrational_vector&) const;
  
  K_BOX3D get_range(const bigrational&, const bigrational&,
                    const bigrational&, const bigrational&) const;
  K_BOX3D get_range(const bigrational_vector&,
                    const bigrational_vector&) const;
  
  friend unsigned long get_all_int_pts(const K_RATPOLY&,
                                       const K_SURF&, const K_SURF&,
                                       K_POINT2D**&);
  friend int match_pts(K_POINT2D** const, const unsigned long, K_SURF* const,
                       K_POINT2D** const, const unsigned long, K_SURF* const);
  
  K_SURF split_surf(const bigrational&, const unsigned long) const;
  
  friend K_SOLID gen_new_solid(K_PARTITION**, const unsigned long,
                               K_PARTITION**, const unsigned long);
  
  friend int get_patch4(const bigrational_vector [4], K_PATCH*&);
  
  friend K_SOLID gen_box(const bigrational_vector [], const unsigned long);
  
  friend K_SOLID gen_cyl(const bigrational_vector&, const bigrational_vector&,
                         const bigrational_vector&, const bigrational_vector&,
                         const bigrational_vector&, const bigrational_vector&);
  
  friend K_SOLID gen_ell(const bigrational_vector [4]);
  
  friend K_SOLID gen_tor(const bigrational_vector&,
                         const bigrational_vector&,
                         const bigrational_vector&,
                         const bigrational_vector&,
                         const bigrational&, const bigrational&);
};

unsigned long get_plane_coeffs(const bigrational_vector&,
                               const bigrational_vector&,
                               const bigrational_vector&,
                               bigrational*&);

#endif

