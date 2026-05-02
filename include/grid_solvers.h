#ifndef GRID_SOLVERS_H_
#define GRID_SOLVERS_H_
#include "types.h"



/* Compute a computational grid in the (xi, eta) space */
point_2D **computational_grid (int, int, int);

/* Compute an interpolated 2D algebraic grid given boundary points */
point_2D **algebraic_grid_2D (int, int, point_2D **, point_2D **);

/* Compute first derivatives of all (x, y) grid points wrt (xi, eta) */
void grid_first_ders_2D (int, int, point_2D **, grid_der_2D ***);

/* Compute coefficients and Jacobian from first derivatives */
void first_der_coefficients (int, int, grid_der_2D **, coeffs_1 ***);

/* Compute finite difference stencil indices */
int **finite_difference_2 (int, int, int);

/* Compute finite difference approximation at a point */
point_2D compute_finite_diff_2 (int, int, int, point_2D **, int);

/* Compute second derivatives of all (x, y) grid points wrt (xi, eta) */
void grid_second_ders_2D (int, int, point_2D **, grid_dder_2D ***);

/* Compute second derivatives of pq wrt (xi, eta) */
void pq_second_ders_2D (int, int, point_2D **, grid_dder_2D ***);

/* Compute Neumann boundary conditions on (xi, eta) boundaries */
point_2D **compute_xbounds_dn (int, int, int, grid_der_2D **, coeffs_1 **, point_2D **);
point_2D **compute_ybounds_dn (int, int, int, grid_der_2D **, coeffs_1 **, point_2D **);

/* Compute pq boundary conditions */
void compute_pq_boundaries (int, int, int, grid_der_2D **, coeffs_1 **,
                            grid_dder_2D **, long double, long double, point_2D ***);
void compute_pq_boundaries_TFI (int, int, grid_der_2D **, coeffs_1 **,
                                grid_dder_2D **, long double, long double, point_2D ***);

/* Initialize (p, q) */
void initial_pq_2D (int, int, int, long double, long double, point_2D **, point_2D ***);
void initial_pq_2D_TFI (int, int, long double, long double, point_2D **, point_2D ***);
void initial_pq_2D_laplace (int, int, long double, long double, point_2D **, point_2D ***);

/* Compute maximum iteration residual */
long double compute_max_residual (int, int, point_2D **, point_2D **);

/* Print unit normals at grid boundaries */
void print_normals (int, int, grid_der_2D **);

/* Compute Thomas-Middlecoff control functions */
void compute_thomas_middlecoff (int, int, grid_der_2D **, grid_dder_2D **, point_2D ***);

/* Compute Steger-Sorenson control functions */
void compute_steger_sorenson (int, int, int, coeffs_1 **, grid_der_2D **, grid_dder_2D **, point_2D ***);

/* Compute residual norm */
long double compute_residual_norm (int, int, point_2D *);

/* Compute 2D elliptic grid */
point_2D **elliptic_grid_2D (int, int, point_2D **, int);
point_2D **elliptic_grid_2D_point (int, int, point_2D **, int);

/* Compute 2D Poisson grid */
point_2D **poisson_grid_2D (int, int, point_2D **, int);
point_2D **poisson_grid_2D_point (int, int, point_2D **, int);
point_2D **poisson_grid_2D_point_lim (int, int, point_2D **, int);

/* Compute 2D biharmonic grid */
point_2D **biharmonic_grid_2D (int, int, point_2D **, int);



#endif // GRID_SOLVERS_H_
