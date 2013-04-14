//  file:    kratpoly.h
//  update:  10/04/02

#ifndef _KRATPOLY_H
#define _KRATPOLY_H

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <bigrational.h>
#include <bigrational_vector.h>
#include <bigrational_matrix.h>

#include <kfloatpoly.h>

using namespace std;

class K_RATPOLY
{
//  A class for representing polynomials with bigrational coefficients
  
//friend class PseudorootList;
friend class ROOT1;
friend class K_POINT1D;
friend class K_POINT2D;
friend class K_CURVE;
friend class K_SURF;
friend class K_PATCH;
friend class K_SOLID;
  
  unsigned long num_vars;    //  number of variables
  long*         deg;         //  max. degrees in variables
  
  unsigned long num_coeffs;  //  number of coefficients
  bigrational*  coeffs;
  
  K_RATPOLY* Sturm_seq;      //  pointer to the Sturm sequence
  
  unsigned long ref_count;   //  reference counter
  
  //  stream
  
  ostream& output(ostream&) const;
  
  //  primitive
  
  //  long* index_to_powers(const unsigned long i) const
  //    returns p = (p[0], p1[1], ..., p[num_vars - 1]) s.t.
  //      coeffs[i] is the coefficient of the monomial X^p of *this.
  
  long* index_to_powers(const unsigned long) const;
  
  //  unsigned long index_to_total_deg(const unsigned long i) const
  //    returns t s.t.
  //      coeffs[i] is the coefficient of the monomial of *this
  //                                        of total degree t.
  
  unsigned long index_to_total_deg(const unsigned long) const;
  
  //  unsigned long powers_to_index(const long* const p) const
  //    returns i s.t.
  //      coeffs[i] is the coefficient of the monomial X^p of *this.
  
  unsigned long powers_to_index(const long* const) const;
  
  //  K_RATPOLY add_var(const unsigned long i) const
  //    returns a polynomial of num_vars + 1 variables of degree
  //      deg[0], deg[1], ..., deg[i - 1], 0, deg[i], ..., deg[num_vars - 1].
  //    e.g. add_variable(1) applied to P(X, Y) returns P(X, Z, Y).
  
  K_RATPOLY add_var(const unsigned long) const;
  
  //  K_RATPOLY remove_var(const unsigned long i) const
  //    PROVIDED deg[i] == 0,
  //    returns a polynomial of num_vars - 1 variables of degree
  //      deg[0], deg[1], ..., deg[i - 1], deg[i + 1], ..., deg[num_vars - 1].
  //  e.g. remove_variable(2) applied to P(X, Y, Z) returns P(X, Y)
  //       provided deg_Z P(X, Y, Z) == 0.
  
  K_RATPOLY remove_var(const unsigned long) const;
  
  //  arithmetic
  
  K_RATPOLY add(const K_RATPOLY&) const;
  K_RATPOLY sub(const K_RATPOLY&) const;
  K_RATPOLY mul(const K_RATPOLY&) const;
  K_RATPOLY mul(const bigrational&) const;
  
  K_RATPOLY neg() const;
  
  K_RATPOLY exact_div(const K_RATPOLY&) const;
  
  //  comparison
  
  //  int cmp(const K_RATPOLY& P) const
  //    returns 0     if *this and P are the same. i.e.,
  //                       deg and coeffs of *this and P are the same, and
  //    returns non-0 otherwise.
  
  int cmp(const K_RATPOLY&) const;
  
  //  evaluation
  
  //  K_RATPOLY subst_val_first_var_proto(const bigrational& b,
  //                                      const int reduced) const
  //    returns a polynomial of num_vars - 1 variables obtained by
  //    computing *this(b, X_1, ..., X_{num_vars - 1}), and
  //    renaming X_1, ..., X_{num_vars - 1} to X_0, X_1, ..., X_{num_vars - 2}.
  //    reduce_deg() is applied if reduced == 1.
  
  K_RATPOLY subst_val_first_var_proto(const bigrational&, const int) const;
  
  //  int fp_sgn_at(const bigrational& b) const
  //    PROVIDED *this is a univariate polynomial,
  //    returns   1 if *this(b) is certainly positive,
  //            - 1 if *this(b) is certainly negative, and
  //              0 otherwise.
  
  int fp_sgn_at(const bigrational&) const;
  
  //  K_RATPOLY subst_val_proto(const unsigned long i, const bigrational& b,
  //                            const int reduced) const
  //    returns a polynomial of num_vars - 1 variables obtained by
  //    computing
  //      *this(X_0, X_1, ..., X_{i - 1}, b, X_{i + 1}, ..., X_{num_vars - 1})
  //    and, renaming X_{i + 1}, ..., X_{num_vars - 1} to
  //                    X_i, ..., X_{num_vars - 2}.
  //    reduce_deg() is applied if reduced == 1.
  
  K_RATPOLY subst_val_proto(const unsigned long, const bigrational&,
                            const int) const;
  
  //  univariate algebra
  
  //  K_RATPOLY div(const K_RATPLOY& P, K_RATPOLY& R) const
  //    PROVIDED *this and P are univariate polynomials,
  //    computes univariate polynomials Q and R s.t.
  //      *this = P * Q + R with deg R < deg P
  //    returns Q.
  
  K_RATPOLY div(const K_RATPOLY&, K_RATPOLY&) const;
  
  //  K_RATPOLY* set_Sturm_seq()
  //    PROVIDED *this is a univariate polynomial,
  //    generates Sturm sequence of *this.
  
  K_RATPOLY*  set_Sturm_seq();
  
  //  K_RATPOLY gcd1(const K_RATPOLY& P) const
  //    PROVIDED *this and P are univariate polynomials,
  //    returns their gcd computed by using Euclidean algorithm.
  
  K_RATPOLY gcd1(const K_RATPOLY&) const;
  
  //  bigrational Sylvester1(const K_RATPOLY& P) const
  //    PROVIDED *this and P are univariate polynomials,
  //    returns Sylvester resultant for *this and P.
  
  bigrational Sylvester1(const K_RATPOLY&) const;
  
  K_RATPOLY monic() const;
  K_RATPOLY monic_gcd1(const K_RATPOLY&) const;
  
  //  bivariate algebra
  
  //  K_RATPOLY Sylvester2(const K_RATPOLY& P, const unsigned long i) const
  //    PROVIDED *this and P are bivariate polynomials,
  //    returns Sylvester resultant for *this and P w.r.t. X_i.
  
  K_RATPOLY Sylvester2(const K_RATPOLY&, const unsigned long) const;
  
  K_RATPOLY gcd2_pp(const K_RATPOLY&) const;
  K_RATPOLY gcd2(const K_RATPOLY&) const;
  
  //  K_RATPOLY GoodSylvester(const K_RATPOLY& P, const unsigned long i) const
  //    PROVIDED *this and P are bivariate polynomials,
  //    returns Sylvetser resultant for *this and P w.r.t. X_i.
  
  K_RATPOLY GoodSylvester(const K_RATPOLY&, const unsigned long) const;
  
public:
  
  //  constructors, assignment and destructor
  
  //  K_RATPOLY()
  //    constructs a zero polynomial.
  
  K_RATPOLY();
  
  //  K_RATPOLY(const unsigned long nv, const long* d)
  //    constructs a polynomial
  //      of nv variables of degree d[0], d[1], ..., d[nv - 1]
  //      with anonymous coefficients.
  
  K_RATPOLY(const unsigned long, const long* const);
  
  //  K_RATPOLY(const unsigned long nv, const long* const d,
  //            const unsigned long nc, const bigrational* const c)
  //    constructs a polynomial
  //      of nv variables of degree d[0], d[1], ..., d[nv - 1]
  //      with coefficients c[0], c[1], ..., c[nc - 1].
  
  K_RATPOLY(const unsigned long, const long* const,
            const unsigned long, const bigrational* const);
  
  //  K_RATPOLY(const unsigned long nv,
  //            const unsigned long i, const bigrational b)
  //    constructs a polynomial X_i - b of nv variables.
  
  K_RATPOLY(const unsigned long n, const unsigned long i,
            const bigrational& b);
  
  K_RATPOLY(const K_RATPOLY&);
  K_RATPOLY& operator =(const K_RATPOLY&);
  
  ~K_RATPOLY();
  
  //  stream
  
  friend ostream& operator <<(ostream&, const K_RATPOLY&);
  
  //  primitive
  
  //  unsigned long get_num_vars() const
  //    returns the number of variables of *this.
  
  unsigned long get_num_vars() const;
  
  //  unsigned long get_total_deg() const
  //    returns the total degree of *this.
  
  unsigned long get_total_deg() const;
  
  //  bigrational& get_coeff(const long* const p) const
  //    returns the coefficient of the monomial X^p of *this.
  
  bigrational& get_coeff(const long* const) const;
  
  //  int reduce_deg()
  //    reduces deg[0], deg[1], ..., deg[num_vars - 1]
  //      to max. degrees in X_0, X_1, ..., X_{num_vars - 1}, resp.
  //    returns 1 if some deg[i] is reduced, and
  //            0 otherwise.
  
  int reduce_deg();
  
  int reduce_num_coeffs();
  //  Make all the coefficients to integers.
  
  int reduce_coeffs();
  //  Make all the coefficients to integers.
  
  //  arithmetic
  
  friend K_RATPOLY operator +(const K_RATPOLY&, const K_RATPOLY&);
  friend K_RATPOLY operator -(const K_RATPOLY&, const K_RATPOLY&);
  friend K_RATPOLY operator *(const K_RATPOLY&, const K_RATPOLY&);
  friend K_RATPOLY operator *(const K_RATPOLY&, const bigrational&);
  friend K_RATPOLY operator *(const bigrational&, const K_RATPOLY&);
  
  friend K_RATPOLY operator -(const K_RATPOLY&);
  
  //  K_RATPOLY exact_div(const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    returns P1 / P2 PROVIDED P1 is divisible by P2.
  
  friend K_RATPOLY exact_div(const K_RATPOLY&, const K_RATPOLY&);
  
  //  K_RATPOLY derivative(const unsigned long i) const
  //     returns the partial derivative of *this w.r.t. X_i.
  
  K_RATPOLY derivative(const unsigned long) const;
  
  //  comparison
  
  friend int operator ==(const K_RATPOLY&, const K_RATPOLY&);
  
  //  int is_zero() const
  //    returns 1 if *this is a zero polynomial, and
  //            0 otherwise.
  
  int is_zero() const;
  
  //  int eq_upto_const(const K_RATPOLY& P) const
  //    returns 1 if *this == P * b for some bigrational b, and
  //            0 otherwise.
  
  int eq_upto_const(const K_RATPOLY&) const;
  
  //  evaluation
  
  //  K_RATPOLY subst_val_first_var(const bigrational& b) const
  //    returns a polynomial of num_vars - 1 variables obtained by
  //    computing *this(b, X_1, ..., X_{num_vars - 1}), and
  //    renaming X_1, ..., X_{num_vars - 1} to X_0, X_1, ..., X_{num_vars - 2}.
  
  K_RATPOLY subst_val_first_var(const bigrational&) const;
  
  //  bigrational evaluate(const bigrational& b) const
  //    returns *this(b) PROVIDED *this is a univariate polynomial.
  
  bigrational evaluate(const bigrational&) const;
  
  //  bigrational evaluate(const bigrational* const B) const
  //    returns *this(B[0], B[1], ..., B[num_vars - 1]).
  
  bigrational evaluate(const bigrational* const) const;
  
  //  bigrational evaluate(const bigrational_vector& B) const
  //    returns *this(B[0], B[1], ..., B[num_vars - 1]).
  
  bigrational evaluate(const bigrational_vector&) const;
  
  //  int sgn_at(const bigrational& b) const
  //    PROVIDED *this is a univariate polynomial,
  //    returns the sign of *this(b).
  
  int sgn_at(const bigrational&) const;
  
  //  int sgn_at(const bigrational* const B) const
  //    PROVIDED *this is a univariate polynomial,
  //    returns the sign of *this(B[0], B[1], ..., B[num_vars - 1]).
  
  int sgn_at(const bigrational* const) const;
  
  //  int sgn_at(const bigrational_vector& B) const
  //    PROVIDED *this is a univariate polynomial,
  //    returns the sign of *this(B[0], B[1], ..., B[num_vars - 1]).
  
  int sgn_at(const bigrational_vector&) const;
  
  //  K_RATPOLY subst_val(const unsigned long i, const bigrational& b) const
  //    returns a polynomial of num_vars - 1 variables obtained by
  //    computing
  //      *this(X_0, X_1, ..., X_{i - 1}, b, X_{i + 1}, ..., X_{num_vars - 1})
  //    and, renaming
  //           X_{i + 1}, ..., X_{num_vars - 1} to X_i, ..., X_{num_vars - 2}.
  
  K_RATPOLY subst_val(const unsigned long, const bigrational&) const;
  
  //  K_RATPOLY subst_expr(const unsigned long i, const K_RATPOLY& P) const
  //    returns a polynomial of num_vars - 1 variables obtained by
  //    computing
  //      *this(X_0, X_1, ..., X_{i - 1},
  //            P (X_0, X_1, ..., X_{i - 1}, X_{i + 1}, X_{num_vars - 1}),
  //            X_{i + 1}, ..., X_{num_vars - 1})
  //    and, renaming
  //           X_{i + 1}, ..., X_{num_vars - 1} to X_i, ..., X_{num_vars - 2}.
  
  K_RATPOLY subst_expr(const unsigned long, const K_RATPOLY&) const;
  
  //  K_RATPOLY subst_expr(const unsigned long i,
  //                       const K_RATPOLY& N, const K_RATPOLY& D) const
  //    returns a polynomial of num_vars - 1 variables obtained by
  //    computing
  //      *this(X_0, X_1, ..., X_{i - 1},
  //            N/D (X_0, X_1, ..., X_{i - 1}, X_{i + 1}, X_{num_vars - 1}),
  //            X_{i + 1}, ..., X_{num_vars - 1})
  //    and, renaming
  //           X_{i + 1}, ..., X_{num_vars - 1} to X_i, ..., X_{num_vars - 2}.
  
  K_RATPOLY subst_expr(const unsigned long,
                       const K_RATPOLY&, const K_RATPOLY&) const;
  
  //  K_RATPOLY subst_param_expr(const K_RATPOLY& X,
  //                             const K_RATPOLY& Y,
  //                             const K_RATPOLY& Z,
  //                             const K_RATPOLY& W) const
  //    PROVIDED *this is a 3-variate polynomial and
  //             X, Y, Z, W are bivariate polynomials,
  //    returns *this(X/W, Y/W, Z/W).
  
  K_RATPOLY subst_param_expr(const K_RATPOLY& X,
                             const K_RATPOLY& Y,
                             const K_RATPOLY& Z,
                             const K_RATPOLY& W) const;
  
  //  other
  
  //  int eval_range(const bigrational* const low_in,
  //                 const bigrational* const high_in,
  //                 bigrational& low_out,
  //                 bigrational& high_out)
  //    computes the range [low_out, high_out]
  //    which *this is evaluated at some value in [low_in, high_in]
  //    using affine arithmetic.
  
  int eval_range(const bigrational* const low_in,
                 const bigrational* const high_in,
                 bigrational& low_out,
                 bigrational& high_out);
  
  //  int eval_range(const bigrational_vector& low_in,
  //                 const bigrational_vector& high_in,
  //                 bigrational& low_out,
  //                 bigrational& high_out)
  //    computes the range [low_out, high_out]
  //    which *this is evaluated at some value in [low_in, high_in]
  //    using affine arithmetic.
  
  int eval_range(const bigrational_vector& low_in,
                 const bigrational_vector& high_in,
                 bigrational& low_out,
                 bigrational& high_out);
  
  //  K_FLOATPOLY as_FLOATPOLY() const
  //    returns a polynomial whose coefficients are
  //                           floating point approximations of those of *this.
  
  K_FLOATPOLY as_FLOATPOLY() const;
  
  K_RATPOLY conv_to_Bernstein(const long, const long) const;
  
  //  univariate algebra
  
  //  K_RATPOLY div(const K_RATPOLY& P1, const K_RATPLOY& P2,
  //                K_RATPOLY& R) const
  //    PROVIDED P1 and P2 are univariate polynomials,
  //    computes univariate polynomials Q and R s.t.
  //      P1 = P2 * Q + R with deg R < deg P2
  //    returns Q.
  
  friend K_RATPOLY div(const K_RATPOLY&, const K_RATPOLY&, K_RATPOLY&);
  
  friend K_RATPOLY operator /(const K_RATPOLY&, const K_RATPOLY&);
  friend K_RATPOLY rem(const K_RATPOLY&, const K_RATPOLY&);
  
  //  long num_Sturm_seq_perm(const bigrational& b) const
  //    PROVIDED *this is a univariate polynomial,
  //    returns
  //      the number of sign permanencies in the Sturm sequence of *this at b.
  
  long num_Sturm_seq_perm(const bigrational&) const;
  
  //  bigrational get_Mignotte_bd() const
  //    PROVIDED *this is a univariate polynomial,
  //    returns an upper bound on the size of the largest root of *this.
  
  bigrational get_Mignotte_bd() const;
  
  //  univariate and bivariate algebra
  
  //  K_RATPOLY gcd(const K_RATPOLY& P1, const K_RATPOLY& P2)
  //    PROVIDED P1 and P2 are uni or bivariate polynomials,
  //    returns their gcd.
  
  friend K_RATPOLY gcd(const K_RATPOLY&, const K_RATPOLY&);
  
  //  bigrational Sylvester1(const K_RATPOLY& P1, const K_RATPOLY& P2) const
  //    PROVIDED P1 and P2 are univariate polynomials,
  //    returns Sylvester resultant for P1 and P2.
  
  friend bigrational Sylvester1(const K_RATPOLY&, const K_RATPOLY&);
  
  //  bivariate algebra
  
  //  K_RATPOLY GoodSylvester(const K_RATPOLY& P1, const K_RATPOLY& P2,
  //                          const unsigned long i) const
  //    PROVIDED P1 and P2 are bivariate polynomials,
  //    returns Sylvetser resultant for P1 and P2 w.r.t. X_i.
  
  friend K_RATPOLY GoodSylvester(const K_RATPOLY&, const K_RATPOLY&,
                                 const unsigned long i);
  
  //  get_pts(): get algebraic points
  
  friend unsigned long get_pts(const bigrational&, const bigrational&,
                               const K_RATPOLY&,
                               K_POINT1D**&,
                               const bigrational&, const int);
  
  friend unsigned long get_pts(const bigrational&, const bigrational&,
                               const bigrational&, const bigrational&,
                               const K_RATPOLY&, const K_RATPOLY&,
                               K_POINT2D**&,
                               const bigrational&, const int);
  
  friend unsigned long get_pts_interior(const bigrational&, const bigrational&,
                                        const bigrational&, const bigrational&,
                                        const K_RATPOLY&, const K_RATPOLY&,
                                        K_POINT2D**&,
                                        const bigrational&);
  
  friend int refine_interior(K_POINT1D* const, K_POINT1D* const,
                             const K_RATPOLY&, const K_RATPOLY&,
                             K_POINT2D*&,
                             const bigrational&);
  
  friend unsigned long get_pts_proto(const bigrational&, const bigrational&,
                                     const K_RATPOLY&,
                                     const bigrational&,
                                     const K_RATPOLY&, const K_RATPOLY&,
                                     K_POINT2D**&,
                                     const bigrational&, const int);
  
  friend unsigned long get_pts_proto(const bigrational&,
                                     const bigrational&, const bigrational&,
                                     const K_RATPOLY&,
                                     const K_RATPOLY&, const K_RATPOLY&,
                                     K_POINT2D**&,
                                     const bigrational&, const int);
  
  //  gen_curve_topo()
  
  friend unsigned long gen_curve_topo(const K_RATPOLY&,
                                      const bigrational&, const bigrational&,
                                      const bigrational&, const bigrational&,
                                      K_CURVE**&);
  
  friend unsigned long gen_curve_topo_proto(const K_RATPOLY&,
                                            const bigrational&,
                                            const bigrational&,
                                            const bigrational&,
                                            const bigrational&,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_POINT2D** const,
                                            const unsigned long,
                                            K_CURVE**&);
  
  friend int implicitize(const K_RATPOLY&,
                         const K_RATPOLY&,
                         const K_RATPOLY&,
                         const K_RATPOLY&,
                         const long,
                         K_RATPOLY*&);
  
  //  gen_box(), gen_cyl(), gen_ell(), gen_tor()
  
  friend int get_param_tor(const bigrational_vector&,
                           const bigrational_vector&,
                           const bigrational_vector&,
                           const bigrational_vector&,
                           const bigrational&, const bigrational&,
                           K_RATPOLY*&, K_RATPOLY*&, K_RATPOLY*&, K_RATPOLY*&);
  
  K_RATPOLY transform_Impl(const bigrational_matrix&) const;
  
  friend K_RATPOLY read_poly(istream&);
};

#endif

