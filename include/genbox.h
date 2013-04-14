#ifndef _GENBOX_H
#define _GENBOX_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational_vector.h>

#include <ksolid.h>

using namespace std;

//  int get_impl_plane_bilin(const bigrational_vector pts[4], K_RATPOLY*& impl)
//    computes the implicit formula for the plane of the bilinear surface
//      that passes 4 points pts.
//    returns 1 if 4 points pts are coplanar and
//            0 otherwise, i.e., 4 points pts are on some bilinear surface.

int get_impl_plane_bilin(const bigrational_vector [4], K_RATPOLY*&);

//  int get_param_plane(const bigrational_vector& x,
//                      const bigrational_vector& y,
//                      const bigrational_vector& z,
//                      K_RATPOLY*& X,
//                      K_RATPOLY*& Y,
//                      K_RATPOLY*& Z,
//                      K_RATPOLY*& W)
//    computes the parametric formula (X/W, Y/W, Z/W) for the plane
//      that passes 3 points x, y and z.

int get_param_plane(const bigrational_vector&,
                    const bigrational_vector&,
                    const bigrational_vector&,
                    K_RATPOLY*&,
                    K_RATPOLY*&,
                    K_RATPOLY*&,
                    K_RATPOLY*&);

K_SOLID read_box(istream&, const bigrational& = 0);
K_SOLID read_BRLCAD_box(istream&, const bigrational& = 0);

#endif

