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


    // Compute analytic Jacobian dFdu
    dFdu = allocate_2D_ell_jacobian_2D_array ("dFdu", nx * ny, nx * ny);
    ell_construct_newton_matrix (nx, ny, dxi, deta, dgrid, coeffs, d2grid, &dFdu);


    // Test point and its flat index
    i1 = 2;
    j1 = 2;
    k1 = i1 + (j1 * nx);


    // Sweep the 9-point stencil around (i1, j1)
    for (dj = -1; dj <= 1; dj++)
    {
        for (di = -1; di <= 1; di++)
        {
            i2 = i1 + di;
            j2 = j1 + dj;
            k2 = i2 + (j2 * nx);

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
    free_2D_ell_jacobian_2D_array ("dFdu", nx * ny, dFdu);
}




/*
    Boundary-row check for ell_construct_newton_matrix on the same
    clustered quadratic test grid.

    Boundary points have Dirichlet conditions in this Newton system.
    The corresponding rows of dFdu encode the constraint
    "u(boundary) is fixed", which gives:

        dFdu[k1][k1] = identity_ders_ell ()    (= [[1, 0], [0, 1]])
        dFdu[k1][k2] = zero ell_jacobian_2D    for all k2 != k1

    where k1 corresponds to a boundary index (i1 = 0, nx-1, or
    j1 = 0, ny-1).

    This test asserts that pattern at every boundary row across all
    columns. Catches errors in the boundary branch of the matrix
    constructor (loop bounds, mis-set diagonal, missed off-diagonals).

    Tolerance is tight (1e-15L) since the values are exact constants
    coming from identity_ders_ell () and zero (ell_jacobian_2D), not
    derived from any FD computation.
*/
Test(newton, ell_construct_newton_matrix_boundary)
{
    // Grid and discretization parameters
    int                     nx = 5, ny = 5, i, j;
    long double             a = 0.3L, b, c = 0.1L;
    long double             dxi, deta;
    long double             xi, eta;

    // Grid arrays and analytic Jacobian
    point_2D                **grid;
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    ell_jacobian_2D         **dFdu;

    // Loop indices
    int                     k1, k2, i1, j1;


    // Setup
    b = -a;
    dxi = 1.0L/(long double) (nx - 1);
    deta = 1.0L/(long double) (ny - 1);


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


    // Compute analytic Jacobian dFdu
    dFdu = allocate_2D_ell_jacobian_2D_array ("dFdu", nx * ny, nx * ny);
    ell_construct_newton_matrix (nx, ny, dxi, deta, dgrid, coeffs, d2grid, &dFdu);


    // Sweep all rows; check boundary ones
    for (k1 = 0; k1 < nx * ny; k1++)
    {
        i1 = k1 % nx;
        j1 = k1/nx;

        // Skip strictly interior rows
        if (i1 >= 1 && i1 <= nx - 2 && j1 >= 1 && j1 <= ny - 2)
        {
            continue;
        }

        // Boundary row: every column k2 must be either identity (k2 == k1)
        // or zero (k2 != k1)
        for (k2 = 0; k2 < nx * ny; k2++)
        {
            if (k2 == k1)
            {
                cr_assert_float_eq (dFdu[k1][k2].ddx.x, 1.0L, 1e-15L,
                                    "Boundary diagonal at k1=%d (i1=%d, j1=%d) ddx.x not 1",
                                    k1, i1, j1);
                cr_assert_float_eq (dFdu[k1][k2].ddy.y, 1.0L, 1e-15L,
                                    "Boundary diagonal at k1=%d (i1=%d, j1=%d) ddy.y not 1",
                                    k1, i1, j1);
                cr_assert_float_eq (dFdu[k1][k2].ddx.y, 0.0L, 1e-15L,
                                    "Boundary diagonal at k1=%d (i1=%d, j1=%d) ddx.y not 0",
                                    k1, i1, j1);
                cr_assert_float_eq (dFdu[k1][k2].ddy.x, 0.0L, 1e-15L,
                                    "Boundary diagonal at k1=%d (i1=%d, j1=%d) ddy.x not 0",
                                    k1, i1, j1);
            }
            else
            {
                cr_assert_float_eq (dFdu[k1][k2].ddx.x, 0.0L, 1e-15L,
                                    "Boundary off-diagonal k1=%d, k2=%d ddx.x not 0",
                                    k1, k2);
                cr_assert_float_eq (dFdu[k1][k2].ddx.y, 0.0L, 1e-15L,
                                    "Boundary off-diagonal k1=%d, k2=%d ddx.y not 0",
                                    k1, k2);
                cr_assert_float_eq (dFdu[k1][k2].ddy.x, 0.0L, 1e-15L,
                                    "Boundary off-diagonal k1=%d, k2=%d ddy.x not 0",
                                    k1, k2);
                cr_assert_float_eq (dFdu[k1][k2].ddy.y, 0.0L, 1e-15L,
                                    "Boundary off-diagonal k1=%d, k2=%d ddy.y not 0",
                                    k1, k2);
            }
        }
    }


    // Cleanup
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_ell_jacobian_2D_array ("dFdu", nx * ny, dFdu);
}
