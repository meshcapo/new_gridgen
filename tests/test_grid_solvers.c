#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <criterion/logging.h>
#include <math.h>
#include "types.h"
#include "aux_functions.h"
#include "grid_solvers.h"
#include "geometry.h"


/* 
    Test that computational_grid with unit = 0 
    produces a normalized uniform grid on 
    [0, 1] x [0, 1].

    Using a 5x3 grid:
        xi values: 0, 0.25, 0.5, 0.75, 1.0 (i = 0, 1, 2, 3, 4)
        eta values: 0, 0.5, 1.0            (j = 0, 1, 2)

    Checks all grid points
*/
Test(grid_solvers, computational_grid_normalized)
{
    // Grid dimensions 
    int                     nx = 5, ny = 3, i, j;
    point_2D                **grid;

    // Call function 
    grid = computational_grid (nx, ny, 0);

    // Check all points 
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            cr_assert_float_eq (grid[i][j].x, (long double) i/(long double) (nx - 1), 1e-15L,
                                "Computational grid failure at grid[%d][%d].x", i, j);
            cr_assert_float_eq (grid[i][j].y, (long double) j/(long double) (ny - 1), 1e-15L,
                                "Computational grid failure at grid[%d][%d].y", i, j);
        }
    }

    // Cleanup
    free_2D_point_2D_array ("grid", nx, grid);
}

/*
    Test that computational_grid with unit = 1 
    produces the correct uniform grid 

    Using a 5x3 grid: grid[i][j] should be (i, j) as long double 
*/
Test(grid_solvers, computational_grid_unit)
{
    // Grid dimensions 
    int                     nx = 5, ny = 3, i, j;
    point_2D                **grid;

    // Call function 
    grid = computational_grid (nx, ny, 1);

    // Check all points 
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            cr_assert_float_eq (grid[i][j].x, (long double) i, 1e-15L,
                                "Computational grid failure at grid[%d][%d].x", i, j);
            cr_assert_float_eq (grid[i][j].y, (long double) j, 1e-15L,
                                "Computational grid failure at grid[%d][%d].y", i, j);
        }
    }

    // Cleanup
    free_2D_point_2D_array ("grid", nx, grid);
}


/*
    Test algebraic_grid_2D on a parallelogram domain with 
    corners (0, 0), (1, 0), (1.5, 1), (0.5, 1) - a unit 
    square sheared by 0.5 in x

    Boundaries:
        eta = 0 (bottom): x = xi, y = 0 
        eta = 1 (top): x = xi + 0.5, y = 1 
        xi = 0 (left): x = 0.5 * eta, y = eta 
        xi = 1 (right): x = 1 + (0.5 * eta), y = eta 

    Analytical TFI solution: x(xi, eta) = xi + (0.5 * eta), y(xci, eta) = eta 
*/
Test(grid_solvers, algebraic_grid_2D_pgram)
{
    int                     nx = 5, ny = 3, i, j;
    long double             xi, eta;
    point_2D                **x_bounds, **y_bounds, **grid;

    // Build boundary arrays 
    x_bounds = allocate_2D_point_2D_array ("x_bounds", 2, nx);
    y_bounds = allocate_2D_point_2D_array ("y_bounds", 2, ny);

    for (i = 0; i < nx; i++)    // Bottom boundary
    {
        x_bounds[0][i].x = (long double) i/(long double) (nx - 1);
        x_bounds[0][i].y = 0.0L;
    }
    for (i = 0; i < nx; i++)    // Top boundary
    {
        x_bounds[1][i].x = ((long double) i/(long double) (nx - 1)) + 0.5L;
        x_bounds[1][i].y = 1.0L;
    }
    for (j = 0; j < ny; j++)    // Left boundary
    {
        y_bounds[0][j].x = ((long double) j/(long double) (ny - 1)) * 0.5L;
        y_bounds[0][j].y = (long double) j/(long double) (ny - 1);

    }
    for (j = 0; j < ny; j++)    // Right boundary 
    {
        y_bounds[1][j].x = 1.0L + (((long double) j/(long double) (ny - 1)) * 0.5L);
        y_bounds[1][j].y = (long double) j/(long double) (ny - 1);

    }

    // Call function
    grid = algebraic_grid_2D (nx, ny, x_bounds, y_bounds);

    // Check all points against analytical solution 
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            xi = (long double) i/(long double) (nx - 1);
            eta = (long double) j/(long double) (ny - 1);
            cr_assert_float_eq (grid[i][j].x, xi + (0.5L * eta), 1e-15L, 
                                "Algebraic grid generation failed at grid[%d][%d].x", i, j);
            cr_assert_float_eq (grid[i][j].y, eta, 1e-15L, 
                                "Algebraic grid generation failed at grid[%d][%d].y", i, j);
        }
    }

    // Cleanup 
    free_2D_point_2D_array ("x_bounds", 2, x_bounds);
    free_2D_point_2D_array ("y_bounds", 2, y_bounds);
    free_2D_point_2D_array ("grid", nx, grid);
}

/*
    Test algebraic_grid_2D on a parallelogram domain with 
    exponential clustering applied to the bottom and top 
    boundaries (eta = 0 and eta = 1). Left and right 
    boundaries remain uniformly spaced 

    Boundaries:
        eta = 0 (bottom): x = f(xi_i), y = 0 
        eta = 1 (top): x = f(xi_i) + 0.5, y = 1 
        xi = 0 (left): x = 0.5 * eta_j, y = eta_j
        xi = 1 (right): x = 1 + (0.5 * eta_j), y = eta_j

    where f = exp_clustering (nx, 1, 2.0, 0.0, 1.0)

    Analytical TFI solution: x(i, j) = f(xi_i) + (0.5 * eta_j), y(i, j) = eta_j
*/
Test(grid_solvers, algebraic_grid_2D_pgram_clustered)
{
    int                     nx = 5, ny = 3, i, j;
    long double             eta, *f;
    point_2D                **x_bounds, **y_bounds, **grid;

    // Generate clustered xi distribution
    f = exp_clustering (nx, 1, 2.0L, 0.0L, 1.0L);

    // Build boundary arrays
    x_bounds = allocate_2D_point_2D_array ("x_bounds", 2, nx);
    y_bounds = allocate_2D_point_2D_array ("y_bounds", 2, ny);

    for (i = 0; i < nx; i++)    // Bottom boundary
    {
        x_bounds[0][i].x = f[i];
        x_bounds[0][i].y = 0.0L;
    }
    for (i = 0; i < nx; i++)    // Top boundary
    {
        x_bounds[1][i].x = f[i] + 0.5L;
        x_bounds[1][i].y = 1.0L;
    }
    for (j = 0; j < ny; j++)    // Left boundary
    {
        y_bounds[0][j].x = ((long double) j/(long double) (ny - 1)) * 0.5L;
        y_bounds[0][j].y = (long double) j/(long double) (ny - 1);

    }
    for (j = 0; j < ny; j++)    // Right boundary 
    {
        y_bounds[1][j].x = 1.0L + (((long double) j/(long double) (ny - 1)) * 0.5L);
        y_bounds[1][j].y = (long double) j/(long double) (ny - 1);

    }

    // Call function
    grid = algebraic_grid_2D (nx, ny, x_bounds, y_bounds);

    // Check all points against analytical solution 
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j/(long double) (ny - 1);
            cr_assert_float_eq (grid[i][j].x, f[i] + (0.5L * eta), 1e-12L, 
                                "Algebraic grid generation failed at grid[%d][%d].x", i, j);
            cr_assert_float_eq (grid[i][j].y, eta, 1e-15L, 
                                "Algebraic grid generation failed at grid[%d][%d].y", i, j);
        }
    }

    // Cleanup 
    free_1D_long_double_array ("clustered_pts", f);
    free_2D_point_2D_array ("x_bounds", 2, x_bounds);
    free_2D_point_2D_array ("y_bounds", 2, y_bounds);
    free_2D_point_2D_array ("grid", nx, grid);
}


/*
    Test grid_first_ders_2D on a parallelogram domain 
    with corners (0, 0), (1, 0), (1.5, 1), (0.5, 1)

    Boundaries:
        eta = 0 (bottom): x = xi, y = 0 
        eta = 1 (top): x = xi + 0.5, y = 1 
        xi = 0 (left): x = 0.5 * eta, y = eta 
        xi = 1 (right): x = 1 + (0.5 * eta), y = eta 

    Since x and y are linear in xi and eta, all finite 
    difference schemes give exact results. Expected 
    results:
        x_der = (dx/dxi, dy/dxi) = (1, 0)
        y_der = (dx/deta, dy/deta) = (0.5, 1)
*/
Test(grid_solvers, grid_first_ders_2D_pgram)
{
    int                     nx = 5, ny = 3, i, j;
    long double             xi, eta;
    point_2D                **grid;
    grid_der_2D             **dgrid;


    // Build parallelogram grid with x = xi + (0.5 * eta), 
    // y = eta
    grid = allocate_2D_point_2D_array ("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            xi = (long double) i/(long double) (nx - 1);
            eta = (long double) j/(long double) (ny - 1);
            grid[i][j].x = xi + (0.5L * eta);
            grid[i][j].y = eta;
        }
    }

    // Allocate and compute derivatives 
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    // Check all points 
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            cr_assert_float_eq (dgrid[i][j].x_der.x, 1.0L, 1e-15L, 
                                "Incorrect grid derivative dx/dxi at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].x_der.y, 0.0L, 1e-15L, 
                                "Incorrect grid derivative dy/dxi at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].y_der.x, 0.5L, 1e-15L, 
                                "Incorrect grid derivative dx/deta at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].y_der.y, 1.0L, 1e-15L, 
                                "Incorrect grid derivative dy/deta at (%d, %d)", i, j);
        }
    }

    // Cleanup
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
}

/*
    Test grid_first_ders_2D on the smoothbump domain using a
    bilinear TFI grid:
        x(xi, eta) = -1.5 + (3 * xi)
        y(xi, eta) = ((1 - eta) * f_b(x)) + (0.8 * eta)
    where f_b(x) = 0.0625 * exp(-25 * x^2)

    Analytical derivatives:
        dx/dxi  = 3 (exact, x linear in xi)
        dy/dxi  = 3 * (1 - eta) * f_b'(x) (O(h^2) truncation error)
        dx/deta = 0 (exact, x independent of eta)
        dy/deta = 0.8 - f_b(x) (exact, y linear in eta)
    where f_b'(x) = -3.125 * x * exp(-25 * x^2)
*/
Test(grid_solvers, grid_first_ders_2D_smoothbump)
{
    int             nx = 101, ny = 11, i, j;
    long double     xi, eta, x, f_b, f_b_prime;
    point_2D        **grid;
    grid_der_2D     **dgrid;

    /* Build smoothbump TFI grid */
    grid = allocate_2D_point_2D_array ("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i / (long double)(nx - 1);
        x = -1.5L + (3.0L * xi);
        f_b = 0.0625L * expl (-25.0L * x * x);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j / (long double)(ny - 1);
            grid[i][j].x = x;
            grid[i][j].y = ((1.0L - eta) * f_b) + (0.8L * eta);
        }
    }

    /* Allocate and compute derivatives */
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    /* Check all points against analytical derivatives */
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i / (long double)(nx - 1);
        x = -1.5L + (3.0L * xi);
        f_b = 0.0625L * expl (-25.0L * x * x);
        f_b_prime = -3.125L * x * expl (-25.0L * x * x);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j / (long double)(ny - 1);

            cr_assert_float_eq (dgrid[i][j].x_der.x, 3.0L, 1e-14L,
                                "Incorrect grid derivative dx/dxi at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].x_der.y, 3.0L * (1.0L - eta) * f_b_prime, 2e-2L,
                                "Incorrect grid derivative dy/dxi at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].y_der.x, 0.0L, 1e-14L,
                                "Incorrect grid derivative dx/deta at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].y_der.y, 0.8L - f_b, 1e-14L,
                                "Incorrect grid derivative dy/deta at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
}

/*
    Test grid_first_ders_2D on a quadratically clustered grid
    with the same parallelogram corners (0,0), (1,0), (1.5,1), (0.5,1).

    Grid: x(i,j) = xi_i^2 + (0.5 * eta_j^2), y(i,j) = eta_j^2
    where xi_i = i/(nx-1), eta_j = j/(ny-1)

    Analytical derivatives (exact for all FD schemes since
    x and y are quadratic in xi and eta):
        dx/dxi  = 2 * xi_i
        dy/dxi  = 0
        dx/deta = eta_j
        dy/deta = 2 * eta_j
*/
Test(grid_solvers, grid_first_ders_2D_pgram_clustered)
{
    /* Grid dimensions */
    int             nx = 5, ny = 3, i, j;
    long double     xi, eta;
    point_2D        **grid;
    grid_der_2D     **dgrid;

    /* Build quadratically clustered grid */
    grid = allocate_2D_point_2D_array ("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i / (long double)(nx - 1);
        for (j = 0; j < ny; j++)
        {
            eta             = (long double) j / (long double)(ny - 1);
            grid[i][j].x   = (xi * xi) + (0.5L * eta * eta);
            grid[i][j].y   = eta * eta;
        }
    }

    /* Allocate and compute derivatives */
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    /* Check all points — exact for quadratic functions */
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i / (long double)(nx - 1);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j / (long double)(ny - 1);

            cr_assert_float_eq (dgrid[i][j].x_der.x, 2.0L * xi, 1e-15L,
                                "Incorrect grid derivative dx/dxi at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].x_der.y, 0.0L, 1e-15L,
                                "Incorrect grid derivative dy/dxi at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].y_der.x, eta, 1e-15L,
                                "Incorrect grid derivative dx/deta at (%d, %d)", i, j);
            cr_assert_float_eq (dgrid[i][j].y_der.y, 2.0L * eta, 1e-15L,
                                "Incorrect grid derivative dy/deta at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
}


/*
    Test grid_second_ders_2D on the uniform parallelogram domain
    with corners (0, 0), (1, 0), (1.5, 1), (0.5, 1)

    Grid: x(i,j) = xi_i + (0.5 * eta_j), y(i,j) = eta_j

    Since x and y are linear in xi and eta, all six second
    derivatives are zero everywhere
*/
Test(grid_solvers, grid_second_ders_2D_pgram)
{
    int             nx = 5, ny = 5, i, j;
    long double     xi, eta;
    point_2D        **grid;
    grid_dder_2D    **d2grid;

    /* Build parallelogram grid with x = xi + (0.5 * eta), y = eta */
    grid = allocate_2D_point_2D_array ("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i / (long double)(nx - 1);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j / (long double)(ny - 1);
            grid[i][j].x = xi + (0.5L * eta);
            grid[i][j].y = eta;
        }
    }

    /* Allocate and compute second derivatives */
    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);

    /* Check all points — all second derivatives are zero for linear x and y */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            cr_assert_float_eq (d2grid[i][j].x_der.x, 0.0L, 1e-15L,
                                "Incorrect second derivative d2x/dxi2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].x_der.y, 0.0L, 1e-15L,
                                "Incorrect second derivative d2y/dxi2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].xy_der.x, 0.0L, 1e-15L,
                                "Incorrect second derivative d2x/dxideta at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].xy_der.y, 0.0L, 1e-15L,
                                "Incorrect second derivative d2y/dxideta at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].y_der.x, 0.0L, 1e-15L,
                                "Incorrect second derivative d2x/deta2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].y_der.y, 0.0L, 1e-15L,
                                "Incorrect second derivative d2y/deta2 at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
}

/*
    Test grid_second_ders_2D on the smoothbump domain using a
    bilinear TFI grid:
        x(xi, eta) = -1.5 + (3 * xi)
        y(xi, eta) = ((1 - eta) * f_b(x)) + (0.8 * eta)
    where f_b(x) = 0.0625 * exp(-25 * x^2)

    Analytical second derivatives:
        d2x/dxi2    = 0 (exact, x linear in xi)
        d2y/dxi2    = 9 * (1 - eta) * f_b''(x) (O(h^2) truncation error)
        d2x/dxideta = 0 (exact, dx/dxi = 3 is constant)
        d2y/dxideta = -3 * f_b'(x) (O(h^2) truncation error)
        d2x/deta2   = 0 (exact, x independent of eta)
        d2y/deta2   = 0 (exact, dy/deta = 0.8 - f_b(x) independent of eta)
    where f_b'(x)  = -3.125 * x * exp(-25 * x^2)
          f_b''(x) = ((156.25 * x^2) - 3.125) * exp(-25 * x^2)
*/
Test(grid_solvers, grid_second_ders_2D_smoothbump)
{
    int             nx = 10001, ny = 11, i, j;
    long double     xi, eta, x, f_b_prime, f_b_double_prime;
    point_2D        **grid;
    grid_dder_2D    **d2grid;

    /* Build smoothbump TFI grid */
    grid = allocate_2D_point_2D_array ("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i / (long double)(nx - 1);
        x  = -1.5L + (3.0L * xi);
        for (j = 0; j < ny; j++)
        {
            eta          = (long double) j / (long double)(ny - 1);
            grid[i][j].x = x;
            grid[i][j].y = ((1.0L - eta) * 0.0625L * expl (-25.0L * x * x)) + (0.8L * eta);
        }
    }

    /* Allocate and compute second derivatives */
    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);

    /* Check all points against analytical second derivatives */
    for (i = 0; i < nx; i++)
    {
        xi               = (long double) i / (long double)(nx - 1);
        x                = -1.5L + (3.0L * xi);
        f_b_prime        = -3.125L * x * expl (-25.0L * x * x);
        f_b_double_prime = ((156.25L * (x * x)) - 3.125L) * expl (-25.0L * x * x);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j/(long double)(ny - 1);

            cr_assert_float_eq (d2grid[i][j].x_der.x, 0.0L, 1e-9L,
                                "Incorrect second derivative d2x/dxi2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].x_der.y, 9.0L * (1.0L - eta) * f_b_double_prime, 1e-4L,
                                "Incorrect second derivative d2y/dxi2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].xy_der.x, 0.0L, 1e-9L,
                                "Incorrect second derivative d2x/dxideta at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].xy_der.y, -3.0L * f_b_prime, 1e-4L,
                                "Incorrect second derivative d2y/dxideta at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].y_der.x, 0.0L, 1e-9L,
                                "Incorrect second derivative d2x/deta2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].y_der.y, 0.0L, 1e-9L,
                                "Incorrect second derivative d2y/deta2 at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
}

/*
    Test grid_second_ders_2D on a quadratically clustered grid
    with the same parallelogram corners (0,0), (1,0), (1.5,1), (0.5,1)

    Grid: x(i,j) = xi_i^2 + (0.5 * eta_j^2), y(i,j) = eta_j^2
    where xi_i = i/(nx-1), eta_j = j/(ny-1)

    Analytical second derivatives (exact for all FD schemes since
    x and y are quadratic in xi and eta):
        d2x/dxi2    = 2
        d2y/dxi2    = 0
        d2x/dxideta = 0
        d2y/dxideta = 0
        d2x/deta2   = 1
        d2y/deta2   = 2
*/
Test(grid_solvers, grid_second_ders_2D_pgram_clustered)
{
    int             nx = 5, ny = 5, i, j;
    long double     xi, eta;
    point_2D        **grid;
    grid_dder_2D    **d2grid;

    /* Build quadratically clustered grid */
    grid = allocate_2D_point_2D_array ("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i / (long double)(nx - 1);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j / (long double)(ny - 1);
            grid[i][j].x = (xi * xi) + (0.5L * (eta * eta));
            grid[i][j].y = eta * eta;
        }
    }

    /* Allocate and compute second derivatives */
    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);

    /* Check all points — exact for quadratic functions */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            cr_assert_float_eq (d2grid[i][j].x_der.x, 2.0L, 1e-15L,
                                "Incorrect second derivative d2x/dxi2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].x_der.y, 0.0L, 1e-15L,
                                "Incorrect second derivative d2y/dxi2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].xy_der.x, 0.0L, 1e-15L,
                                "Incorrect second derivative d2x/dxideta at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].xy_der.y, 0.0L, 1e-15L,
                                "Incorrect second derivative d2y/dxideta at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].y_der.x, 1.0L, 1e-15L,
                                "Incorrect second derivative d2x/deta2 at (%d, %d)", i, j);
            cr_assert_float_eq (d2grid[i][j].y_der.y, 2.0L, 1e-15L,
                                "Incorrect second derivative d2y/deta2 at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
}


/*
    Test compute_thomas_middlecoff (P, Q) on the clustered
    quadratic test grid

        x(i, j) = xi_i  + a * xi_i^2  + c * eta_j^2
        y(i, j) = eta_j + b * eta_j^2 + c * xi_i^2

    with a = 0.3, b = -a, c = 0.1, nx = ny = 5,
    xi_i = i/(nx - 1), eta_j = j/(ny - 1)

    Choosing b = -a keeps the SS control functions additively
    separable in (xi, eta) (used by the SS tests below); for TM
    no such restriction is needed but the same grid is used
    here for consistency

    Grid metric quantities (with b = -a):
        r_xi  = (1 + 2a*xi, 2c*xi)
        r_eta = (2c*eta, 1 - 2a*eta)
        r_xixi   = (2a, 2c)
        r_etaeta = (2c, -2a)
        gamma(xi)  = (1 + 2a*xi)^2  + 4c^2*xi^2
        alpha(eta) = 4c^2*eta^2     + (1 - 2a*eta)^2

    Thomas-Middlecoff (P, Q):
        r_xi  . r_xixi   =  2a + 4(a^2 + c^2)*xi
        r_eta . r_etaeta = -2a + 4(a^2 + c^2)*eta
        P(xi)  = -[ 2a + 4(a^2 + c^2)*xi ]/gamma(xi)^2
        Q(eta) = -[-2a + 4(a^2 + c^2)*eta]/alpha(eta)^2
*/
Test(grid_solvers, thomas_middlecoff_pq)
{
    /* Grid dimensions and constants */
    int             nx = 5, ny = 5, i, j;
    long double     a = 0.3L, b, c = 0.1L;
    long double     xi, eta, gamma_ij, alpha_ij, ac2;
    long double     g1, a1, P_exact, Q_exact;
    point_2D        **grid, **pq;
    grid_der_2D     **dgrid;
    grid_dder_2D    **d2grid;

    b = -a;
    ac2 = (a * a) + (c * c);

    /* Build test grid */
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

    /* Compute first derivatives */
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    /* Compute second derivatives */
    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);

    /* Compute Thomas-Middlecoff control functions */
    pq = allocate_2D_point_2D_array ("pq", nx, ny);
    compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &pq);

    /* Check all points against analytic formula */
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i/(long double) (nx - 1);
        g1 = 1.0L + (2.0L * a * xi);
        gamma_ij = (g1 * g1) + (4.0L * c * c * xi * xi);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j/(long double) (ny - 1);
            a1 = 1.0L + (2.0L * b * eta);
            alpha_ij = (4.0L * c * c * eta * eta) + (a1 * a1);
            P_exact = -((2.0L * a) + (4.0L * ac2 * xi))/(gamma_ij * gamma_ij);
            Q_exact = -((2.0L * b) + (4.0L * ac2 * eta))/(alpha_ij * alpha_ij);

            cr_assert_float_eq (pq[i][j].x, P_exact, 1e-14L,
                                "Incorrect Thomas-Middlecoff P at (%d, %d)", i, j);
            cr_assert_float_eq (pq[i][j].y, Q_exact, 1e-14L,
                                "Incorrect Thomas-Middlecoff Q at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_point_2D_array ("pq", nx, pq);
}


/*
    Test compute_thomas_middlecoff_alt (alpha-phi form) on the
    clustered quadratic test grid

        x(i, j) = xi_i  + a * xi_i^2  + c * eta_j^2
        y(i, j) = eta_j + b * eta_j^2 + c * xi_i^2

    with a = 0.3, b = -a, c = 0.1, nx = ny = 5,
    xi_i = i/(nx - 1), eta_j = j/(ny - 1)

    Grid metric quantities (with b = -a):
        gamma(xi)  = (1 + 2a*xi)^2  + 4c^2*xi^2
        alpha(eta) = 4c^2*eta^2     + (1 - 2a*eta)^2

    Thomas-Middlecoff formulas in (phi, psi) form (assumes beta = 0):
        r_xi  . r_xixi   =  2a + 4(a^2 + c^2)*xi
        r_eta . r_etaeta = -2a + 4(a^2 + c^2)*eta
        phi(xi)  = -[ 2a + 4(a^2 + c^2)*xi ]/gamma(xi)
        psi(eta) = -[-2a + 4(a^2 + c^2)*eta]/alpha(eta)

    phi depends only on xi and psi only on eta, so TFI in the
    interior recovers the analytic formula exactly. FD stencils are
    also exact since the grid coordinates are quadratic in (xi, eta)
*/
Test(grid_solvers, thomas_middlecoff_phispi)
{
    /* Grid dimensions and constants */
    int             nx = 5, ny = 5, i, j;
    long double     a = 0.3L, b, c = 0.1L;
    long double     xi, eta, gamma_ij, alpha_ij, ac2;
    long double     g1, a1, phi_exact, psi_exact;
    point_2D        **grid, **phipsi;
    grid_der_2D     **dgrid;
    grid_dder_2D    **d2grid;

    b = -a;
    ac2 = (a * a) + (c * c);

    /* Build test grid */
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

    /* Compute first derivatives */
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    /* Compute second derivatives */
    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);

    /* Compute Thomas-Middlecoff control functions (alpha-phi form) */
    phipsi = allocate_2D_point_2D_array ("phipsi", nx, ny);
    compute_thomas_middlecoff_alt (nx, ny, dgrid, d2grid, &phipsi);

    /* Check all points against analytic formula */
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i/(long double) (nx - 1);
        g1 = 1.0L + (2.0L * a * xi);
        gamma_ij = (g1 * g1) + (4.0L * c * c * xi * xi);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j/(long double) (ny - 1);
            a1 = 1.0L + (2.0L * b * eta);
            alpha_ij = (4.0L * c * c * eta * eta) + (a1 * a1);
            phi_exact = -((2.0L * a) + (4.0L * ac2 * xi))/gamma_ij;
            psi_exact = -((2.0L * b) + (4.0L * ac2 * eta))/alpha_ij;

            cr_assert_float_eq (phipsi[i][j].x, phi_exact, 1e-14L,
                                "Incorrect Thomas-Middlecoff phi at (%d, %d)", i, j);
            cr_assert_float_eq (phipsi[i][j].y, psi_exact, 1e-14L,
                                "Incorrect Thomas-Middlecoff psi at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_point_2D_array ("phipsi", nx, phipsi);
}


/*
    Test compute_steger_sorenson with mode = 0 (P, Q form) on the
    clustered quadratic test grid

        x(i, j) = xi_i  + a * xi_i^2  + c * eta_j^2
        y(i, j) = eta_j + b * eta_j^2 + c * xi_i^2

    with a = 0.3, b = -a, c = 0.1, nx = ny = 5,
    xi_i = i/(nx - 1), eta_j = j/(ny - 1)

    Choosing b = -a kills the cross term 4c(a + b)*xi (and its
    eta counterpart) in r_xi.r_etaeta and r_eta.r_xixi, leaving
    (phi, psi) additively separable as f(xi) + g(eta) so TFI in
    the interior recovers the analytic formula exactly

    Grid metric quantities (with b = -a):
        gamma(xi)  = (1 + 2a*xi)^2  + 4c^2*xi^2
        alpha(eta) = 4c^2*eta^2     + (1 - 2a*eta)^2
        r_xi  . r_xixi   =  2a + 4(a^2 + c^2)*xi
        r_eta . r_etaeta = -2a + 4(a^2 + c^2)*eta
        r_xi  . r_etaeta = 2c
        r_eta . r_xixi   = 2c

    Steger-Sorenson formulas in (phi, psi) form (assumes beta = 0):
        phi(xi, eta) = -[ 2a + 4(a^2 + c^2)*xi ]/gamma(xi)  - 2c/alpha(eta)
        psi(xi, eta) = -[-2a + 4(a^2 + c^2)*eta]/alpha(eta) - 2c/gamma(xi)

    With mode = 0 the function divides phi by gamma and psi by alpha
    to produce (P, Q):
        P(xi, eta) = phi(xi, eta)/gamma(xi)
        Q(xi, eta) = psi(xi, eta)/alpha(eta)
*/
Test(grid_solvers, steger_sorenson_pq)
{
    /* Grid dimensions and constants */
    int             nx = 5, ny = 5, i, j;
    long double     a = 0.3L, b, c = 0.1L;
    long double     xi, eta, gamma_ij, alpha_ij, ac2;
    long double     g1, a1, phi_val, psi_val, P_exact, Q_exact;
    point_2D        **grid, **pq;
    grid_der_2D     **dgrid;
    grid_dder_2D    **d2grid;
    coeffs_1        **coeffs;

    b = -a;
    ac2 = (a * a) + (c * c);

    /* Build test grid */
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

    /* Compute first derivatives */
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    /* Compute first-derivative coefficients (alpha, beta, gamma, J) */
    coeffs = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    first_der_coefficients (nx, ny, dgrid, &coeffs);

    /* Compute second derivatives */
    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);

    /* Compute Steger-Sorenson control functions in (P, Q) form */
    pq = allocate_2D_point_2D_array ("pq", nx, ny);
    compute_steger_sorenson (nx, ny, 0, coeffs, dgrid, d2grid, &pq);

    /* Check all points against analytic formula */
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i/(long double) (nx - 1);
        g1 = 1.0L + (2.0L * a * xi);
        gamma_ij = (g1 * g1) + (4.0L * c * c * xi * xi);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j/(long double) (ny - 1);
            a1 = 1.0L + (2.0L * b * eta);
            alpha_ij = (4.0L * c * c * eta * eta) + (a1 * a1);
            phi_val = (-((2.0L * a) + (4.0L * ac2 * xi))/gamma_ij) - ((2.0L * c)/alpha_ij);
            psi_val = (-((2.0L * b) + (4.0L * ac2 * eta))/alpha_ij) - ((2.0L * c)/gamma_ij);
            P_exact = phi_val/gamma_ij;
            Q_exact = psi_val/alpha_ij;

            cr_assert_float_eq (pq[i][j].x, P_exact, 1e-14L,
                                "Incorrect Steger-Sorenson P at (%d, %d)", i, j);
            cr_assert_float_eq (pq[i][j].y, Q_exact, 1e-14L,
                                "Incorrect Steger-Sorenson Q at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_point_2D_array ("pq", nx, pq);
}


/*
    Test compute_steger_sorenson with mode = 1 (phi, psi form) on the
    clustered quadratic test grid

        x(i, j) = xi_i  + a * xi_i^2  + c * eta_j^2
        y(i, j) = eta_j + b * eta_j^2 + c * xi_i^2

    with a = 0.3, b = -a, c = 0.1, nx = ny = 5,
    xi_i = i/(nx - 1), eta_j = j/(ny - 1)

    Grid metric quantities (with b = -a):
        gamma(xi)  = (1 + 2a*xi)^2  + 4c^2*xi^2
        alpha(eta) = 4c^2*eta^2     + (1 - 2a*eta)^2
        r_xi  . r_xixi   =  2a + 4(a^2 + c^2)*xi
        r_eta . r_etaeta = -2a + 4(a^2 + c^2)*eta
        r_xi  . r_etaeta = 2c
        r_eta . r_xixi   = 2c

    Steger-Sorenson formulas in (phi, psi) form (assumes beta = 0):
        phi(xi, eta) = -[ 2a + 4(a^2 + c^2)*xi ]/gamma(xi)  - 2c/alpha(eta)
        psi(xi, eta) = -[-2a + 4(a^2 + c^2)*eta]/alpha(eta) - 2c/gamma(xi)

    With b = -a, (phi, psi) are additively separable in (xi, eta),
    so TFI in the interior recovers the analytic formula exactly.
    Mode = 1 leaves (phi, psi) unscaled, so the coeffs argument is
    not used inside the function but a valid array is still passed
*/
Test(grid_solvers, control_steger_sorenson_phipsi)
{
    /* Grid dimensions and constants */
    int             nx = 5, ny = 5, i, j;
    long double     a = 0.3L, b, c = 0.1L;
    long double     xi, eta, gamma_ij, alpha_ij, ac2;
    long double     g1, a1, phi_exact, psi_exact;
    point_2D        **grid, **phipsi;
    grid_der_2D     **dgrid;
    grid_dder_2D    **d2grid;
    coeffs_1        **coeffs;

    b = -a;
    ac2 = (a * a) + (c * c);

    /* Build test grid */
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

    /* Compute first derivatives */
    dgrid = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    grid_first_ders_2D (nx, ny, grid, &dgrid);

    /* Compute first-derivative coefficients (alpha, beta, gamma, J) */
    coeffs = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    first_der_coefficients (nx, ny, dgrid, &coeffs);

    /* Compute second derivatives */
    d2grid = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    grid_second_ders_2D (nx, ny, grid, &d2grid);

    /* Compute Steger-Sorenson control functions in (phi, psi) form */
    phipsi = allocate_2D_point_2D_array ("phipsi", nx, ny);
    compute_steger_sorenson (nx, ny, 1, coeffs, dgrid, d2grid, &phipsi);

    /* Check all points against analytic formula */
    for (i = 0; i < nx; i++)
    {
        xi = (long double) i/(long double) (nx - 1);
        g1 = 1.0L + (2.0L * a * xi);
        gamma_ij = (g1 * g1) + (4.0L * c * c * xi * xi);
        for (j = 0; j < ny; j++)
        {
            eta = (long double) j/(long double) (ny - 1);
            a1 = 1.0L + (2.0L * b * eta);
            alpha_ij = (4.0L * c * c * eta * eta) + (a1 * a1);
            phi_exact = (-((2.0L * a) + (4.0L * ac2 * xi))/gamma_ij) - ((2.0L * c)/alpha_ij);
            psi_exact = (-((2.0L * b) + (4.0L * ac2 * eta))/alpha_ij) - ((2.0L * c)/gamma_ij);

            cr_assert_float_eq (phipsi[i][j].x, phi_exact, 1e-14L,
                                "Incorrect Steger-Sorenson phi at (%d, %d)", i, j);
            cr_assert_float_eq (phipsi[i][j].y, psi_exact, 1e-14L,
                                "Incorrect Steger-Sorenson psi at (%d, %d)", i, j);
        }
    }

    /* Cleanup */
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_point_2D_array ("phipsi", nx, phipsi);
}
