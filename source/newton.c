#include <stdio.h>
#include "globvar.h"
#include "geometry.h"
#include "newton.h"
#include "jacobian_ders.h"
#include "file_operations.h"
#include "aux_functions.h"
#include "types.h"



/* ----------------------------------------------------------- */
/*                                                             */
/* Functions needed for defining Newton's method for           */
/* biharmonic grid generation                                  */
/*                                                             */
/* ----------------------------------------------------------- */

/*
   Compute Dx and Dy used to compute (p, q)
   at a boundary point (i, j)

   Input parameters: coeff_ij   - first derivative coefficients at
                                  point (i, j)
                     d2grid_ij  - second order derivatives at point (i, j)
*/
point_2D compute_Dx_Dy_ij (coeffs_1 coeff_ij, grid_dder_2D d2grid_ij)
{
    /* Return Dx_Dy */
    point_2D                Dx_Dy;


    /* Dx */
    Dx_Dy.x                             = -((coeff_ij.alpha * d2grid_ij.x_der.x) -
                                            (TWO * coeff_ij.beta * d2grid_ij.xy_der.x) +
                                            (coeff_ij.gamma * d2grid_ij.y_der.x))/
                                            (coeff_ij.J * coeff_ij.J);

    /* Dy */
    Dx_Dy.y                             = -((coeff_ij.alpha * d2grid_ij.x_der.y) -
                                            (TWO * coeff_ij.beta * d2grid_ij.xy_der.y) +
                                            (coeff_ij.gamma * d2grid_ij.y_der.y))/
                                            (coeff_ij.J * coeff_ij.J);


    return Dx_Dy;
}






/*
   Compute Laplacians of xi and eta
   at a boundary point (i, j)

   Input parameters: dgrid_ij   - first derivatives at point (i, j)
                     coeff_ij   - first derivative coefficients at
                                  point (i, j)
                     d2grid_ij  - second order derivatives at point (i, j)
*/
point_2D compute_laplacians_ij (grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                grid_dder_2D d2grid_ij)
{
    /* Return nablas */
    point_2D                nablas;

    /* Local variables */
    point_2D                Dx_Dy;


    /* Compute Dx and Dy */
    Dx_Dy                               = compute_Dx_Dy_ij (coeff_ij, d2grid_ij);


    /* Compute Laplacians */
    nablas.x                            = ((dgrid_ij.y_der.y * Dx_Dy.x) -
                                           (dgrid_ij.y_der.x * Dx_Dy.y))/coeff_ij.J;
    nablas.y                            = ((dgrid_ij.x_der.x * Dx_Dy.y) -
                                           (dgrid_ij.x_der.y * Dx_Dy.x))/coeff_ij.J;


    return nablas;
}




/*
   Compute F(x(i, j)) and F(y(i, j)) of the
   right hand side vector for Newton's method
   in the interior of the grid generation
   domain (2D)

   Input parameters: dgrid_ij   - first derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
                     coeffs_ij  - coefficients derivaed from first
                                  derivatives at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
                     pq_ij      - control function (p, q) values at
                                  point (i, j)
*/
void compute_interior_fx_fy_2D (grid_der_2D dgrid_ij, coeffs_1 coeffs_ij, grid_dder_2D d2grid_ij,
                                point_2D pq_ij, long double *fx, long double *fy)
{
    (*fx)                               = (coeffs_ij.alpha * d2grid_ij.x_der.x) -
                                          (TWO * coeffs_ij.beta * d2grid_ij.xy_der.x) +
                                          (coeffs_ij.gamma * d2grid_ij.y_der.x) +
                                          (coeffs_ij.J * coeffs_ij.J * ((pq_ij.x * dgrid_ij.x_der.x) +
                                          (pq_ij.y * dgrid_ij.y_der.x)));
    (*fy)                               = (coeffs_ij.alpha * d2grid_ij.x_der.y) -
                                          (TWO * coeffs_ij.beta * d2grid_ij.xy_der.y) +
                                          (coeffs_ij.gamma * d2grid_ij.y_der.y) +
                                          (coeffs_ij.J * coeffs_ij.J * ((pq_ij.x * dgrid_ij.x_der.y) +
                                          (pq_ij.y * dgrid_ij.y_der.y)));
}






/*
   Compute F(p(i, j)) and F(q(i, j)) of the
   right hand side vector for Newton's method
   in the interior of the grid generation
   domain (2D)

   Input parameters: coeffs_ij  - coefficients derivaed from first
                                  derivatives at point (i, j)
                     d2pq_ij    - second order derivatives of (p, q) with
                                  respect to (xi, eta) at point (i, j)
*/
void compute_interior_fp_fq_2D (coeffs_1 coeffs_ij, grid_dder_2D d2pq_ij,
                                long double *fp, long double *fq)
{
    (*fp)                               = (coeffs_ij.alpha * d2pq_ij.x_der.x) -
                                          (TWO * coeffs_ij.beta * d2pq_ij.xy_der.x) +
                                          (coeffs_ij.gamma * d2pq_ij.y_der.x);
    (*fq)                               = (coeffs_ij.alpha * d2pq_ij.x_der.y) -
                                          (TWO * coeffs_ij.beta * d2pq_ij.xy_der.y) +
                                          (coeffs_ij.gamma * d2pq_ij.y_der.y);
}






/*
   Compute F(x(i, j)) and F(y(i, j)) of the
   right hand side vector for Newton's method
   on xi = 0 and xi = 1 boundaries (2D)

   Input parameters: coeffs_ij  - coefficients derivaed from first
                                  derivatives at point (i, j)
                     dgrid_ij   - first derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
                     pq_ij      - control function (p, q) values at
                                  point (i, j)
                     ny_ij      - outward normal at point (i, j)
*/
void compute_xi_01_fp_fq_2D (long double c, coeffs_1 coeffs_ij, grid_der_2D dgrid_ij,
                             grid_dder_2D d2grid_ij, point_2D pq_ij, point_2D ny_ij,
                             long double *fp, long double *fq)
{
    /* Local variables */
    long double             Dx, Dy, ptemp, qtemp;


    /* Compute Dx and Dy used to compute
       (p, q) at the boundaries */
    Dx                                  = -((coeffs_ij.alpha * d2grid_ij.x_der.x) -
                                            (TWO * coeffs_ij.beta * d2grid_ij.xy_der.x) +
                                            (coeffs_ij.gamma * d2grid_ij.y_der.x))/(coeffs_ij.J *
                                            coeffs_ij.J);
    Dy                                  = -((coeffs_ij.alpha * d2grid_ij.x_der.y) -
                                            (TWO * coeffs_ij.beta * d2grid_ij.xy_der.y) +
                                            (coeffs_ij.gamma * d2grid_ij.y_der.y))/(coeffs_ij.J *
                                            coeffs_ij.J);

    /* Compute Laplacians of xi and eta */
    ptemp                               = ((dgrid_ij.y_der.y * Dx) - (dgrid_ij.y_der.x * Dy))/
                                          coeffs_ij.J;
    qtemp                               = ((dgrid_ij.x_der.x * Dy) - (dgrid_ij.x_der.y * Dx))/
                                          coeffs_ij.J;


    /* Compute f(p) and f(q) at point (0, j) */
    (*fp)                               = pq_ij.x - ptemp;
    (*fq)                               = pq_ij.y - qtemp + (c * dgrid_ij.x_der.y * ny_ij.x) -
                                          (c * dgrid_ij.x_der.x * ny_ij.y);
}






/*
   Compute F(x(i, j)) and F(y(i, j)) of the
   right hand side vector for Newton's method
   on eta = 0 and eta = 1 boundaries (2D)

   Input parameters: coeffs_ij  - coefficients derivaed from first
                                  derivatives at point (i, j)
                     dgrid_ij   - first derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
                     pq_ij      - control function (p, q) values at
                                  point (i, j)
                     nx_ij      - outward normal at point (i, j)
*/
void compute_eta_01_fp_fq_2D (long double c, coeffs_1 coeffs_ij, grid_der_2D dgrid_ij,
                              grid_dder_2D d2grid_ij, point_2D pq_ij, point_2D nx_ij,
                              long double *fp, long double *fq)
{
    /* Local variables */
    long double             Dx, Dy, ptemp, qtemp;


    /* Compute Dx and Dy used to compute
       (p, q) at the boundaries */
    Dx                                  = -((coeffs_ij.alpha * d2grid_ij.x_der.x) -
                                            (TWO * coeffs_ij.beta * d2grid_ij.xy_der.x) +
                                            (coeffs_ij.gamma * d2grid_ij.y_der.x))/(coeffs_ij.J *
                                            coeffs_ij.J);
    Dy                                  = -((coeffs_ij.alpha * d2grid_ij.x_der.y) -
                                            (TWO * coeffs_ij.beta * d2grid_ij.xy_der.y) +
                                            (coeffs_ij.gamma * d2grid_ij.y_der.y))/(coeffs_ij.J *
                                            coeffs_ij.J);

    /* Compute Laplacians of xi and eta */
    ptemp                               = ((dgrid_ij.y_der.y * Dx) - (dgrid_ij.y_der.x * Dy))/
                                          coeffs_ij.J;
    qtemp                               = ((dgrid_ij.x_der.x * Dy) - (dgrid_ij.x_der.y * Dx))/
                                          coeffs_ij.J;

    /* Compute f(p) and f(q) at point (0, j) */
    (*fp)                               = pq_ij.x - ptemp - (c * dgrid_ij.y_der.y * nx_ij.x) +
                                          (c * dgrid_ij.y_der.x * nx_ij.y);
    (*fq)                               = pq_ij.y - qtemp;
}






/*
   Construct Newton's method LHS matrix
   This contains the definitions of
   the jacobians dF(u1)/du2

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     dxi        - spacing of xi in computational grid
                     deta       - spacing of eta in computational grid
                     pq         - array containing values of control
                                  function (p, q)
                     dgrid      - array containing values of first
                                  derivatives of (x, y) with respect to (xi, eta)
                     coeffs     - array containing values of first
                                  derivative coefficients
                     d2grid     - array containing values of second-order
                                  derivatives of (x, y) with respect to (xi, eta)
                     d2pq       - array containing values of second-order
                                  derivatives of (p, q) with respect to (xi, eta)
*/
void bh_construct_newton_matrix (int nx, int ny, long double dxi, long double deta,
                                 point_2D **pq, grid_der_2D **dgrid, coeffs_1 **coeffs,
                                 grid_dder_2D **d2grid, grid_dder_2D **d2pq, long double cxi,
                                 long double ceta, bh_jacobian_2D ***dFdu)
{
    /* Local variables */
    int                     i1, i2, j1, j2, k1, k2;
    point_2D                Dx_Dy_ij, nablas_ij, **nx_bounds, **ny_bounds;


    /* Read outward normal vectors for all boundaries */
    nx_bounds                           = read_xbound_normals_2D (nx);
    ny_bounds                           = read_ybound_normals_2D (ny);

    /* k1 is used to traverse the rows of
       the newton method matrix */
    for (k1 = 0; k1 < nx * ny; k1++)
    {
        i1                              = k1%nx;
        j1                              = k1/nx;

        /* For interior points */
        if ((i1 >= 1 && i1 <= nx - 2) &&
            (j1 >= 1 && j1 <= ny - 2))
        {
            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = interior_point_ders (i1, j1, i2, j2, dxi, deta, pq[i1][j1],
                                                               dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1],
                                                               d2pq[i1][j1]);
            }
        }

        /* xi = 0 boundary points */
        if (i1 == 0 && (j1 >= 1 && j1 <= ny - 2))
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = xi0_point_ders (i1, j1, i2, j2, dxi, deta, cxi, dgrid[i1][j1],
                                                          coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                          nablas_ij, ny_bounds[0][j1]);
            }
        }

        /* xi = 1 boundary points */
        if (i1 == nx - 1 && (j1 >= 1 && j1 <= ny - 2))
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = xi1_point_ders (i1, j1, i2, j2, dxi, deta, cxi, dgrid[i1][j1],
                                                          coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                          nablas_ij, ny_bounds[1][j1]);
            }
        }

        /* eta = 0 boundary points */
        if ((i1 >= 1 && i1 <= nx - 2) && j1 == 0)
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = eta0_point_ders (i1, j1, i2, j2, dxi, deta, ceta, dgrid[i1][j1],
                                                           coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                           nablas_ij, nx_bounds[0][i1]);
            }
        }

        /* eta = 1 boundary points */
        if ((i1 >= 1 && i1 <= nx - 2) && j1 == ny - 1)
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = eta1_point_ders (i1, j1, i2, j2, dxi, deta, ceta, dgrid[i1][j1],
                                                           coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                           nablas_ij, nx_bounds[1][i1]);
            }
        }

        /* xi = 0, eta = 0 */
        if (i1 == 0 && j1 == 0)
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = xi0_eta0_point_ders (i1, j1, i2, j2, dxi, deta, ceta, dgrid[i1][j1],
                                                               coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                               nablas_ij, nx_bounds[0][i1]);
            }
        }

        /* xi = 1, eta = 0 */
        if (i1 == nx - 1 && j1 == 0)
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = xi1_eta0_point_ders (i1, j1, i2, j2, dxi, deta, ceta, dgrid[i1][j1],
                                                               coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                               nablas_ij, nx_bounds[0][i1]);
            }
        }

        /* xi = 0, eta = 1 */
        if (i1 == 0 && j1 == ny - 1)
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = xi0_eta1_point_ders (i1, j1, i2, j2, dxi, deta, ceta, dgrid[i1][j1],
                                                               coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                               nablas_ij, nx_bounds[1][i1]);
            }
        }

        /* xi = 1, eta = 1 */
        if (i1 == nx - 1 && j1 == ny - 1)
        {
            Dx_Dy_ij                    = compute_Dx_Dy_ij (coeffs[i1][j1], d2grid[i1][j1]);
            nablas_ij                   = compute_laplacians_ij (dgrid[i1][j1], coeffs[i1][j1], d2grid[i1][j1]);

            /* k2 is used to traverse columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* dFdu[k1][k2] = dF(u(i1, j1))/du(i2, j2) */
                (*dFdu)[k1][k2]         = xi1_eta1_point_ders (i1, j1, i2, j2, dxi, deta, ceta, dgrid[i1][j1],
                                                               coeffs[i1][j1], d2grid[i1][j1], Dx_Dy_ij,
                                                               nablas_ij, nx_bounds[1][i1]);
            }
        }
    }
    /* Print dFdu as a matrix */
    //print_as_2D_matrix (nx * ny, dFdu);


    /* Free outward normal vector arrays */
    free_2D_point_2D_array ("nx_bounds", 2, nx_bounds);
    free_2D_point_2D_array ("ny_bounds", 2, ny_bounds);
}






/*
   Construct RHS vector for Newton's
   method (2D)

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dxi    - xi spacing in computational grid
                     deta   - eta spacing in computational grid
                     dgrid  - array containing first derivatives of
                              (x, y) with respect to (xi, eta)
                     coeffs - array containing first derivative
                              coefficients
                     d2grid - array containing second order derivatives
                              of (x, y) with respect to (xi, eta)
                     pq     - array containing values of control
                              functions (p, q)
                     grid   - array containing (x, y) coordinates
*/
void bh_construct_newton_rhs_vector (int nx, int ny, point_2D **grid, point_2D **pq,
                                     grid_der_2D **dgrid, coeffs_1 **coeffs,
                                     grid_dder_2D **d2grid, grid_dder_2D **d2pq,
                                     long double cxi, long double ceta, bh_sol_2D **rhs)
{
    /* Local variables */
    int                     i, j, k;
    point_2D                **nx_bounds, **ny_bounds, **x_bounds, **y_bounds;


    /* Read coordinates along all boundaries */
    x_bounds                            = get_x_bounds_2D (nx, ny, grid);
    y_bounds                            = get_y_bounds_2D (nx, ny, grid);

    /* Read normal vectors for all boundaries */
    nx_bounds                           = read_xbound_normals_2D (nx);
    ny_bounds                           = read_ybound_normals_2D (ny);

    for (k = 0; k < nx * ny; k++)
    {
        i                               = k%nx;
        j                               = k/nx;

        if ((i >= 1 && i <= nx - 2) &&
            (j >= 1 && j <= ny - 2))
        {
            /* Compute f(x) and f(y) for interior
               point (i, j) */
            compute_interior_fx_fy_2D (dgrid[i][j], coeffs[i][j], d2grid[i][j], pq[i][j],
                                    &((*rhs)[k].x), &((*rhs)[k].y));

            /* Compute f(p) and f(q) for interior
               point (i, j) */
            compute_interior_fp_fq_2D (coeffs[i][j], d2pq[i][j], &((*rhs)[k].p), &((*rhs)[k].q));
        }

        /* xi = 0, 1 boundary */
        if ((i == 0 || i == nx - 1) &&
            (j >= 1 && j <= ny - 2))
        {
            /* f(x) and f(y) */
            (*rhs)[k].x                 = grid[i][j].x - y_bounds[i/(nx - 1)][j].x;
            (*rhs)[k].y                 = grid[i][j].y - y_bounds[i/(nx - 1)][j].y;

            /* f(p) and f(q) */
            if (i == 0)
            {
                /* xi = 0 */
                compute_xi_01_fp_fq_2D (cxi, coeffs[i][j], dgrid[i][j], d2grid[i][j], pq[i][j],
                                        ny_bounds[0][j], &((*rhs)[k].p), &((*rhs)[k].q));
            } else {
                /* xi = 1 */
                compute_xi_01_fp_fq_2D (cxi, coeffs[i][j], dgrid[i][j], d2grid[i][j], pq[i][j],
                                        ny_bounds[1][j], &((*rhs)[k].p), &((*rhs)[k].q));
            }
        }

        /* eta = 0, 1 boundary */
        if (j == 0 || j == ny - 1)
        {
            /* f(x) and f(y) */
            (*rhs)[k].x                 = grid[i][j].x - x_bounds[j/(ny - 1)][i].x;
            (*rhs)[k].y                 = grid[i][j].y - x_bounds[j/(ny - 1)][i].y;

            /* f(p) and f(q) */
            if (j == 0)
            {
                /* eta = 0 */
                compute_eta_01_fp_fq_2D (ceta, coeffs[i][j], dgrid[i][j], d2grid[i][j], pq[i][j],
                                         nx_bounds[0][i], &((*rhs)[k].p), &((*rhs)[k].q));
            } else {
                /* eta = 1 */
                compute_eta_01_fp_fq_2D (ceta, coeffs[i][j], dgrid[i][j], d2grid[i][j], pq[i][j],
                                         nx_bounds[1][i], &((*rhs)[k].p), &((*rhs)[k].q));
            }
        }

        /* Compute -f(u) since this is the RHS */
        (*rhs)[k].x                     *= -ONE;
        (*rhs)[k].y                     *= -ONE;
        (*rhs)[k].p                     *= -ONE;
        (*rhs)[k].q                     *= -ONE;
    }


    /* Free outward normal arrays */
    free_2D_point_2D_array ("x_bounds", 2, x_bounds);
    free_2D_point_2D_array ("y_bounds", 2, y_bounds);
    free_2D_point_2D_array ("nx_bounds", 2, nx_bounds);
    free_2D_point_2D_array ("ny_bounds", 2, ny_bounds);
}






/*
   Update the 2D point_2D (x, y) and (p, q)
   arrays after a Newton iteration using
   the bh_sol_2D update array

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     update - 1D array of type bh_sol_2D
                              containing (x, y, p, q) update values
*/
void bh_update_solution (int nx, int ny, long double damping,
                         bh_sol_2D *update, point_2D ***grid,
                         point_2D ***pq)
{
    for (int k = 0; k < nx * ny; k++)
    {
        int i                           = k%nx;
        int j                           = k/nx;

        (*grid)[i][j].x                 += damping * update[k].x;
        (*grid)[i][j].y                 += damping * update[k].y;
        (*pq)[i][j].x                   += damping * update[k].p;
        (*pq)[i][j].y                   += damping * update[k].q;
    }
}






/* ----------------------------------------------------------- */
/*                                                             */
/* Functions needed for defining Newton's method for           */
/* elliptic grid generation                                    */
/*                                                             */
/* ----------------------------------------------------------- */

/*
   Compute F(x(i, j)) and F(y(i, j)) of the
   right hand side vector for Newton's method
   in the interior of the grid generation
   domain (2D)

   Input parameters: coeffs_ij  - coefficients derived from first
                                  derivatives at point (i, j)
                     d2grid_ij  - second order derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
*/
point_2D ell_compute_interior_fx_fy_2D (coeffs_1 coeff_ij, grid_dder_2D d2grid_ij)
{
    /* Return F */
    point_2D                F;


    /* Compute F(x(i, j)) and F(y(i, j)) */
    F.x                                 = (coeff_ij.alpha * d2grid_ij.x_der.x) -
                                          (TWO * coeff_ij.beta * d2grid_ij.xy_der.x) +
                                          (coeff_ij.gamma * d2grid_ij.y_der.x);
    F.y                                 = (coeff_ij.alpha * d2grid_ij.x_der.y) -
                                          (TWO * coeff_ij.beta * d2grid_ij.xy_der.y) +
                                          (coeff_ij.gamma * d2grid_ij.y_der.y);


    return F;
}






/*
   Construct Newton's method LHS matrix
   for elliptic grid generation. This
   contains the definitions of the jacobians
   dF(u1)/du2

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dxi    - spacing of xi in computational grid
                     deta   - spacing of eta in computational grid
                     dgrid  - array containing values of first derivatives
                              of (x, y) with respect to (xi, eta)
                     coeffs - array containing values of first derivative
                              coefficients
                     d2grid - array containing values of second-order
                              derivatives of (x, y) with respect to (xi, eta)
*/
void ell_construct_newton_matrix (int nx, int ny, long double dxi, long double deta,
                                  grid_der_2D **dgrid, coeffs_1 **coeffs, grid_dder_2D **d2grid,
                                  ell_jacobian_2D ***dFdu)
{
    /* Local variables */
    int                     i1, j1, i2, j2, k1, k2;


    /* k1 is used to traverse the rows of
       the newton method matrix */
    for (k1 = 0; k1 < nx * ny; k1++)
    {
        i1                              = k1%nx;
        j1                              = k1/nx;

        /* For interior points */
        if ((i1 >= 1 && i1 <= nx - 2) &&
            (j1 >= 1 && j1 <= ny - 2))
        {
            /* k2 is used to traverse the columns of
               the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                (*dFdu)[k1][k2]         = ell_interior_point_ders (i1, j1, i2, j2, dxi, deta, dgrid[i1][j1],
                                                                   coeffs[i1][j1], d2grid[i1][j1]);
            }
        }
        else /* For boundary points */
        {
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* The derivatives are non-zero only
                   when i2 == i1 and j2 == j1 */
                if (i2 == i1 && j2 == j1)
                {
                    (*dFdu)[k1][k2]     = identity_ders_ell ();
                }
                else
                {
                    (*dFdu)[k1][k2]     = zero (ell_jacobian_2D);
                }
            }
        }
    }
}






/*
   Construct RHS vector for Newton's
   method for elliptic grid generation (2D).
   We are actually computing -F(u) here

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     grid   - array containing (x, y) coordinates
                     coeffs - array containing first derivative coefficients
                     d2grid - array containing second derivative values
*/
void ell_construct_newton_rhs_vector (int nx, int ny, point_2D **grid, coeffs_1 **coeffs,
                                      grid_dder_2D **d2grid, point_2D **rhs, int print_from)
{
    /* Local variables */
    int                     i, j, k;
    point_2D                **x_bounds, **y_bounds;


    /* Read coordinates along all boundaries */
    x_bounds                            = read_x_bounds_2D (nx);
    y_bounds                            = read_y_bounds_2D (ny);


    /* Compute F(x(i, j)) and F(y(i, j)) */
    for (k = 0; k < nx * ny; k++)
    {
        i                               = k%nx;
        j                               = k/nx;

        /* Domain interior */
        if ((i >= 1 && i <= nx - 2) &&
            (j >= 1 && j <= ny - 2))
        {
            (*rhs)[k]                   = ell_compute_interior_fx_fy_2D (coeffs[i][j], d2grid[i][j]);
        }

        /* xi = 0, 1 boundaries */
        if ((i == 0 || i == nx - 1) &&
            (j >= 1 && j <= ny - 2))
        {
            (*rhs)[k].x                 = -grid[i][j].x + y_bounds[i/(nx - 1)][j].x;
            (*rhs)[k].y                 = -grid[i][j].y + y_bounds[i/(nx - 1)][j].y;
        }

        /* eta = 0, 1 boundaries */
        if (j == 0 || j == ny - 1)
        {
            (*rhs)[k].x                 = -grid[i][j].x + x_bounds[j/(ny - 1)][i].x;
            (*rhs)[k].y                 = -grid[i][j].y + x_bounds[j/(ny - 1)][i].y;
        }

        /* Compute -F(u) since this is the RHS */
        (*rhs)[k].x                 *= -ONE;
        (*rhs)[k].y                 *= -ONE;
    }


    /* Free array memory */
    free_2D_point_2D_array ("x_bounds", 2, x_bounds);
    free_2D_point_2D_array ("y_bounds", 2, y_bounds);
}






/*
   Update the 2D point_2D (x, y) and (p, q)
   arrays after a Newton iteration using
   the bh_sol_2D update array

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     update - 1D array of type bh_sol_2D
                              containing (x, y, p, q) update values
*/
void ell_update_solution (int nx, int ny, point_2D *update, point_2D ***grid)
{
    for (int k = 0; k < nx * ny; k++)
    {
        int i                           = k%nx;
        int j                           = k/nx;

        (*grid)[i][j].x                 += update[k].x;
        (*grid)[i][j].y                 += update[k].y;
    }
}






/* ----------------------------------------------------------- */
/*                                                             */
/* Functions needed for defining Newton's method for           */
/* elliptic grid generation with control functions             */
/*                                                             */
/* ----------------------------------------------------------- */

/*
   Compute F(x(i, j)) and F(y(i, j)) of the
   right hand side vector for Newton's method
   in the interior of the grid generation
   domain

   Input parameters: dgrid_ij   - derivatives of (x, y) with respect to
                                  (xi, eta) at point (i, j)
                     coeff_ij   - first derivative coefficients at point (i, j)
                     d2grid_ij  - second-order derivatives of (x, y) with
                                  respect to (xi, eta) at point (i, j)
                     pq_ij      - control function values at point (i, j)
*/
point_2D psn_compute_interior_fx_fy_2D (grid_der_2D dgrid_ij, coeffs_1 coeff_ij,
                                        grid_dder_2D d2grid_ij, point_2D pq_ij)
{
    /* Return ders */
    point_2D                F;


    /* Compute F(x(i, j)) and F(y(i, j)) */
    F.x                                 = (coeff_ij.alpha * d2grid_ij.x_der.x) -
                                          (TWO * coeff_ij.beta * d2grid_ij.xy_der.x) +
                                          (coeff_ij.gamma * d2grid_ij.y_der.x) +
                                          (coeff_ij.J * coeff_ij.J * ((pq_ij.x * dgrid_ij.x_der.x) +
                                           (pq_ij.y * dgrid_ij.y_der.x)));
    F.y                                 = (coeff_ij.alpha * d2grid_ij.x_der.y) -
                                          (TWO * coeff_ij.beta * d2grid_ij.xy_der.y) +
                                          (coeff_ij.gamma * d2grid_ij.y_der.y) +
                                          (coeff_ij.J * coeff_ij.J * ((pq_ij.x * dgrid_ij.x_der.y) +
                                           (pq_ij.y * dgrid_ij.y_der.y)));


    return F;
}






/*
   Construct Newton's method LHS matrix for
   elliptic grid generation with non-zero control
   functions. This contains the definition of the
   jacobians dF(u1)/du2

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dxi    - spacing of xi in computational grid
                     deta   - spacing of eta in computational grid
                     dgrid  - arry containing first derivative values
                     coeffs - array containing first derivative coefficients
                     d2grid - array containing second derivative values
                     pq     - array containing control function values
*/
void psn_construct_newton_matrix (int nx, int ny, long double dxi, long double deta,
                                  grid_der_2D **dgrid, coeffs_1 **coeffs, grid_dder_2D **d2grid,
                                  point_2D **pq, ell_jacobian_2D ***dFdu)
{
    /* Local variables */
    int                     i1, j1, i2, j2, k1, k2;


    /* k1 is used to traverse the rows of
       the newton method matrix */
    for (k1 = 0; k1 < nx * ny; k1++)
    {
        i1                              = k1%nx;
        j1                              = k1/nx;

        /* For interior points */
        if ((i1 >= 1 && i1 <= nx - 2) &&
            (j1 >= 1 && j1 <= ny - 2))
        {
            /* k2 is used to traverse the columns of
            the newton method matrix */
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                (*dFdu)[k1][k2]         = psn_interior_point_ders (i1, j1, i2, j2, dxi, deta, dgrid[i1][j1],
                                                                   coeffs[i1][j1], d2grid[i1][j1], pq[i1][j1]);
            }
        }
        else /* For boundary points */
        {
            for (k2 = 0; k2 < nx * ny; k2++)
            {
                i2                      = k2%nx;
                j2                      = k2/nx;

                /* The derivatives are non-zero only
                   when i2 == i1 and j2 == j1 */
                if (i2 == i1 && j2 == j1)
                {
                    (*dFdu)[k1][k2]     = identity_ders_ell ();
                }
                else
                {
                    (*dFdu)[k1][k2]     = zero (ell_jacobian_2D);
                }
            }
        }
    }
}






/*
   Construct RHS vector for Newton's method
   for elliptic grid generation (2D) with
   non-zero control functions. We are actually
   computing -F(u) here

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     grid   - array containing (x, y) coordinates
                     dgrid  - array containing first derivative values
                     coeffs - array containing first derivative coefficients
                     d2grid - array containg second derivative values
                     pq     - array containing control function values
*/
void psn_construct_newton_rhs_vector (int nx, int ny, point_2D **grid, grid_der_2D **dgrid,
                                      coeffs_1 **coeffs, grid_dder_2D **d2grid, point_2D **pq,
                                      point_2D **rhs)
{
    /* Local variables */
    int                     i, j, k;
    point_2D                **x_bounds, **y_bounds;


    /* Read coordinates along all boundaries */
    x_bounds                            = read_x_bounds_2D (nx);
    y_bounds                            = read_y_bounds_2D (ny);


    /* Compute F(x(i, j)) and F(y(i, j)) */
    for (k = 0; k < nx * ny; k++)
    {
        i                               = k%nx;
        j                               = k/nx;

        /* Domain interior */
        if ((i >= 1 && i <= nx - 2) &&
            (j >= 1 && j <= ny - 2))
        {
            (*rhs)[k]                   = psn_compute_interior_fx_fy_2D (dgrid[i][j], coeffs[i][j], d2grid[i][j],
                                                                         pq[i][j]);
        }

        /* xi = 0, 1 boundaries */
        if ((i == 0 || i == nx - 1) &&
            (j >= 1 && j <= ny - 2))
        {
            (*rhs)[k].x                 = -grid[i][j].x + y_bounds[i/(nx - 1)][j].x;
            (*rhs)[k].y                 = -grid[i][j].y + y_bounds[i/(nx - 1)][j].y;
        }

        /* eta = 0, 1 boundaries */
        if (j == 0 || j == ny - 1)
        {
            (*rhs)[k].x                 = -grid[i][j].x + x_bounds[j/(ny - 1)][i].x;
            (*rhs)[k].y                 = -grid[i][j].y + x_bounds[j/(ny - 1)][i].y;
        }

        /* Compute -F(u) since this is the RHS */
        (*rhs)[k].x                     *= -ONE;
        (*rhs)[k].y                     *= -ONE;
    }


    /* Free memory */
    free_2D_point_2D_array ("x_bounds", 2, x_bounds);
    free_2D_point_2D_array ("y_bounds", 2, y_bounds);
}
