#ifndef JACOBIAN_DERS_H_
#define JACOBIAN_DERS_H_
#include "types.h"



/*-----------------------------------------------------------*/
/*                                                           */
/*     Jacobian derivatives for biharmonic grid generation   */
/*                                                           */
/*-----------------------------------------------------------*/

/* Interior point Jacobians */
bh_jacobian_2D interior_point_ders (int, int, int, int, long double, long double,
                                    point_2D, grid_der_2D, coeffs_1, grid_dder_2D,
                                    grid_dder_2D);

/* xi = 0 Jacobians */
bh_jacobian_2D xi0_point_ders (int, int, int, int, long double, long double,
                               long double, grid_der_2D, coeffs_1, grid_dder_2D,
                               point_2D, point_2D, point_2D);

/* xi = 1 Jacobians */
bh_jacobian_2D xi1_point_ders (int, int, int, int, long double, long double,
                               long double, grid_der_2D, coeffs_1, grid_dder_2D,
                               point_2D, point_2D, point_2D);

/* eta = 0 Jacobians */
bh_jacobian_2D eta0_point_ders (int, int, int, int, long double, long double,
                                long double, grid_der_2D, coeffs_1, grid_dder_2D,
                                point_2D, point_2D, point_2D);

/* eta = 1 Jacobians */
bh_jacobian_2D eta1_point_ders (int, int, int, int, long double, long double,
                                long double, grid_der_2D, coeffs_1, grid_dder_2D,
                                point_2D, point_2D, point_2D);

/* xi = 0, eta = 0 Jacobians */
bh_jacobian_2D xi0_eta0_point_ders (int, int, int, int, long double, long double,
                                    long double, grid_der_2D, coeffs_1, grid_dder_2D,
                                    point_2D, point_2D, point_2D);

/* xi = 1, eta = 0 Jacobians */
bh_jacobian_2D xi1_eta0_point_ders (int, int, int, int, long double, long double,
                                    long double, grid_der_2D, coeffs_1, grid_dder_2D,
                                    point_2D, point_2D, point_2D);

/* xi = 0, eta = 1 Jacobians */
bh_jacobian_2D xi0_eta1_point_ders (int, int, int, int, long double, long double,
                                    long double, grid_der_2D, coeffs_1, grid_dder_2D,
                                    point_2D, point_2D, point_2D);

/* xi = 1, eta = 1 Jacobians */
bh_jacobian_2D xi1_eta1_point_ders (int, int, int, int, long double, long double,
                                    long double, grid_der_2D, coeffs_1, grid_dder_2D,
                                    point_2D, point_2D, point_2D);



/*-----------------------------------------------------------*/
/*                                                           */
/*     Jacobian derivatives for elliptic grid generation     */
/*     Equivalent to setting (p, q) as 0                     */
/*                                                           */
/*-----------------------------------------------------------*/

/* Interior point derivatives */
ell_jacobian_2D ell_interior_point_ders (int, int, int, int, long double, long double,
                                         grid_der_2D, coeffs_1, grid_dder_2D);



/*-----------------------------------------------------------*/
/*                                                           */
/*     Jacobian derivatives for elliptic grid generation     */
/*     with non-zero control functions (p, q)                */
/*                                                           */
/*-----------------------------------------------------------*/

/* Interior point derivatives */
ell_jacobian_2D psn_interior_point_ders (int, int, int, int, long double, long double,
                                         grid_der_2D, coeffs_1, grid_dder_2D, point_2D);



#endif  // JACOBIAN_DERS_H_
