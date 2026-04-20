#ifndef NEWTON_H_
#define NEWTON_H_
#include "types.h"



/*
   BIHARMONIC GRID GENERATION
*/

/* Compute Dx/Dy and laplacians of xi and eta */
point_2D compute_Dx_Dy_ij (coeffs_1, grid_dder_2D);
point_2D compute_laplacians_ij (grid_der_2D, coeffs_1, grid_dder_2D);

/* Compute f(u(i, j)) components in the
   interior of the domain as well as along
   the boundaries */
void compute_interior_fx_fy_2D (grid_der_2D, coeffs_1, grid_dder_2D, point_2D,
                                long double *, long double *);
void compute_interior_fp_fq_2D (coeffs_1, grid_dder_2D, long double *, long double *);
void compute_xi_01_fp_fq_2D (long double, coeffs_1, grid_der_2D, grid_dder_2D, point_2D,
                             point_2D, long double *, long double *);
void compute_eta_01_fp_fq_2D (long double, coeffs_1, grid_der_2D, grid_dder_2D, point_2D,
                              point_2D, long double *, long double *);

/* Construct LHS matrix and RHS vector
   for linearized Newton's method */
void bh_construct_newton_matrix (int, int, long double, long double, point_2D **,
                                 grid_der_2D **, coeffs_1 **, grid_dder_2D **, grid_dder_2D **,
                                 long double, long double, bh_jacobian_2D ***);
void bh_construct_newton_rhs_vector (int, int, point_2D **, point_2D **, grid_der_2D **,
                                     coeffs_1 **, grid_dder_2D **, grid_dder_2D **,
                                     long double, long double, bh_sol_2D **);

/* Convert bh_sol_2D to point_2D arrays
   and update solution */
void bh_update_solution (int, int, long double, bh_sol_2D *, point_2D ***, point_2D ***);



/*
   ELLIPTIC GRID GENERATION
*/

/* Compute f(u(i, j)) components in the
   interior of the domain as well as along
   the boundaries */
point_2D ell_compute_interior_fx_fy_2D (coeffs_1, grid_dder_2D);

/* Construct LHS matrix and RHS vector
   for linearlized Newton's method */
void ell_construct_newton_matrix (int, int, long double, long double, grid_der_2D **,
                                  coeffs_1 **, grid_dder_2D **, ell_jacobian_2D ***);
void ell_construct_newton_rhs_vector (int, int, point_2D **, coeffs_1 **, grid_dder_2D **,
                                      point_2D **, int);

/* Convert 1D point_2D to 2D point_2D arrays
   and update elliptic grid generation solution */
void ell_update_solution (int, int, point_2D *, point_2D ***);



/*
   ELLIPTIC GRID GENERATION WITH CONTROL FUNCTIONS
*/

/* Compute F(u(i, j)) components in the
   interior of the domain as well as along
   the boundaries */
point_2D psn_compute_interior_fx_fy_2D (grid_der_2D, coeffs_1, grid_dder_2D, point_2D);

/* Construct LHS matrix and RHS vector
   for linearized Newton's method */
void psn_construct_newton_matrix (int, int, long double, long double, grid_der_2D **,
                                  coeffs_1 **, grid_dder_2D **, point_2D **, ell_jacobian_2D ***);
void psn_construct_newton_rhs_vector (int, int, point_2D **, grid_der_2D **, coeffs_1 **,
                                      grid_dder_2D **, point_2D **, point_2D **);



#endif // NEWTON_H_
