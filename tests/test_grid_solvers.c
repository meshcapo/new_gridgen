#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <criterion/logging.h>
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
