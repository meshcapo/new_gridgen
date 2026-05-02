#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <criterion/logging.h>
#include <math.h>
#include "types.h"
#include "aux_functions.h"
#include "newton.h"
#include "grid_solvers.h"
#include "geometry.h"




/*
    Integrated FD-vs-analytic check for ell_construct_newton_matrix
    on the clustered quadratic test grid

        x(i, j) = xi  + a xi^2  + c eta^2
        y(i, j) = eta + b eta^2 + c xi^2

    with a = 0.3, b = -a, c = 0.1, nx = ny = 5,
    xi = i/(nx - 1), eta = j/(ny - 1).

    The analytic Jacobian column dFdu[k1][k2] (k1 = i1 + j1*nx,
    k2 = i2 + j2*nx) is the exact derivative of F at (i1, j1) with
    respect to u at (i2, j2), where F is computed from FD-derived
    coeffs and second derivatives via ell_compute_interior_fx_fy_2D.
    The FD column is computed by perturbing u at (i2, j2),
    recomputing dgrid / coeffs / d2grid on the perturbed grid, and
    forming a centered difference of F at (i1, j1). The two columns
    must agree to FD-truncation precision for every neighbor in the
    9-point stencil around (i1, j1).

    Layout of the analytic 2x2 block (per the comments in
    jacobian_ders.c):
        ddx.x = dF_x/dx_neighbor   ddy.x = dF_x/dy_neighbor
        ddx.y = dF_y/dx_neighbor   ddy.y = dF_y/dy_neighbor

    Test point: (i1, j1) = (2, 2) (center of 5x5 grid; 9-pt stencil
    is fully interior, so all dispatched stencil shapes are
    exercised).

    Tolerance: eps = 1e-7L gives FD truncation O(eps^2 * F''') and
    roundoff O(machine_eps * F/eps), summing to ~1e-11 for F
    and derivatives of order 1; tolerance 1e-9L allows two-three
    digits of safety while still being a meaningful check.

    What this catches: dispatcher bugs in ell_interior_point_ders,
    loop-bound errors, sign errors that would make the column
    disagree with d/du F as the code actually computes F. (Catches
    rare misses of the per-stencil unit tests in test_jacobian_ders.c
    by sweeping all 9 stencil shapes simultaneously through one
    integrated check.)
*/
Test(newton, ell_construct_newton_matrix_fd)
{
    // Grid and discretization parameters
    int                     nx = 5, ny = 5, i, j;
    long double             a = 0.3L, b, c = 0.1L;
    long double             dxi, deta, eps, two_eps;
    long double             xi, eta;

    // Grid arrays and analytic Jacobian
    point_2D                **grid;
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    ell_jacobian_2D         **dFdu;

    // FD test point and stencil neighbor indices
    int                     i1, j1, k1;
    int                     di, dj, i2, j2, k2;

    // Storage during perturbation
    long double             orig_x, orig_y;
    point_2D                F_plus, F_minus;
    long double             FD_dx_x, FD_dx_y, FD_dy_x, FD_dy_y;


    // Setup
    b = -a;
    dxi = 1.0L/(long double) (nx - 1);
    deta = 1.0L/(long double) (ny - 1);
    eps = 1.0e-7L;
    two_eps = 2.0L * eps;


    // Build the clustered quadratic grid
    grid = allocate_2D_point_2D_array ("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i/(long double) (nx - 1);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j/(long double) (ny - 1);
            grid[i][j].x = xi  + (a * xi  * xi)  + (c * eta * eta);
            grid[i][j].y = eta + (b * eta * eta) + (c * xi  * xi);
        }
    }


    // Compute first derivatives, coefficients, second derivatives
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    coeffs = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    first_der_coefficients (nx, ny, dgrid, &coeffs);

    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);


    // Compute analytic Jacobian dFdu (interior-only Newton system:
    // boundary points are Dirichlet constants and not part of the matrix).
    // Index map for interior (i, j) with i in [1, nx-2], j in [1, ny-2]:
    //     k = (i - 1) + (j - 1) * (nx - 2)
    int n_int = (nx - 2) * (ny - 2);
    dFdu = allocate_2D_ell_jacobian_2D_array ("dFdu", n_int, n_int);
    ell_construct_newton_matrix (nx, ny, dxi, deta, dgrid, coeffs, d2grid, &dFdu);


    // Test point — center of a 5x5 grid is i1=2, j1=2 (an interior point).
    // The 9-point stencil around it stays entirely in the interior, so all
    // stencil neighbors have columns in the interior-only matrix.
    i1 = 2;
    j1 = 2;
    k1 = (i1 - 1) + ((j1 - 1) * (nx - 2));


    // Sweep the 9-point stencil around (i1, j1)
    for (dj = -1; dj <= 1; dj++)
    {
        for (di = -1; di <= 1; di++)
        {
            i2 = i1 + di;
            j2 = j1 + dj;
            k2 = (i2 - 1) + ((j2 - 1) * (nx - 2));

            orig_x = grid[i2][j2].x;
            orig_y = grid[i2][j2].y;


            // Centered FD with respect to x at (i2, j2)
            grid[i2][j2].x = orig_x + eps;
            grid_first_ders_2D (nx, ny, grid, &dgrid);
            first_der_coefficients (nx, ny, dgrid, &coeffs);
            grid_second_ders_2D (nx, ny, grid, &d2grid);
            F_plus = ell_compute_interior_fx_fy_2D (coeffs[i1][j1], d2grid[i1][j1]);

            grid[i2][j2].x = orig_x - eps;
            grid_first_ders_2D (nx, ny, grid, &dgrid);
            first_der_coefficients (nx, ny, dgrid, &coeffs);
            grid_second_ders_2D (nx, ny, grid, &d2grid);
            F_minus = ell_compute_interior_fx_fy_2D (coeffs[i1][j1], d2grid[i1][j1]);

            grid[i2][j2].x = orig_x;

            FD_dx_x = (F_plus.x - F_minus.x)/two_eps;
            FD_dx_y = (F_plus.y - F_minus.y)/two_eps;


            // Centered FD with respect to y at (i2, j2)
            grid[i2][j2].y = orig_y + eps;
            grid_first_ders_2D (nx, ny, grid, &dgrid);
            first_der_coefficients (nx, ny, dgrid, &coeffs);
            grid_second_ders_2D (nx, ny, grid, &d2grid);
            F_plus = ell_compute_interior_fx_fy_2D (coeffs[i1][j1], d2grid[i1][j1]);

            grid[i2][j2].y = orig_y - eps;
            grid_first_ders_2D (nx, ny, grid, &dgrid);
            first_der_coefficients (nx, ny, dgrid, &coeffs);
            grid_second_ders_2D (nx, ny, grid, &d2grid);
            F_minus = ell_compute_interior_fx_fy_2D (coeffs[i1][j1], d2grid[i1][j1]);

            grid[i2][j2].y = orig_y;

            FD_dy_x = (F_plus.x - F_minus.x)/two_eps;
            FD_dy_y = (F_plus.y - F_minus.y)/two_eps;


            // Compare the analytic 2x2 block to the FD 2x2 block
            cr_assert_float_eq (dFdu[k1][k2].ddx.x, FD_dx_x, 1e-9L,
                                "dF_x/dx at (i1=%d,j1=%d) wrt (i2=%d,j2=%d): analytic=%.15Le, FD=%.15Le",
                                i1, j1, i2, j2, dFdu[k1][k2].ddx.x, FD_dx_x);
            cr_assert_float_eq (dFdu[k1][k2].ddx.y, FD_dx_y, 1e-9L,
                                "dF_y/dx at (i1=%d,j1=%d) wrt (i2=%d,j2=%d): analytic=%.15Le, FD=%.15Le",
                                i1, j1, i2, j2, dFdu[k1][k2].ddx.y, FD_dx_y);
            cr_assert_float_eq (dFdu[k1][k2].ddy.x, FD_dy_x, 1e-9L,
                                "dF_x/dy at (i1=%d,j1=%d) wrt (i2=%d,j2=%d): analytic=%.15Le, FD=%.15Le",
                                i1, j1, i2, j2, dFdu[k1][k2].ddy.x, FD_dy_x);
            cr_assert_float_eq (dFdu[k1][k2].ddy.y, FD_dy_y, 1e-9L,
                                "dF_y/dy at (i1=%d,j1=%d) wrt (i2=%d,j2=%d): analytic=%.15Le, FD=%.15Le",
                                i1, j1, i2, j2, dFdu[k1][k2].ddy.y, FD_dy_y);
        }
    }


    // Cleanup
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_ell_jacobian_2D_array ("dFdu", n_int, dFdu);
}



