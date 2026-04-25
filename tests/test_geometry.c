#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <criterion/logging.h>
#include "aux_functions.h"
#include "geometry.h"
#include "types.h"



/*
   Test that get_x_bounds_2D correctly extracts the
   eta = 0 and eta = 1 boundaries from a manufactured
   grid

   Grid is filled with points where x = i, y = j,
   so expected bounds are:
     bounds[0][i] = {i, 0}      (eta = 0, j = 0)
     bounds[1][i] = {i, 30}     (eta = 1, j = ny-1)
*/
Test(geometry, get_x_bounds_2D)
{
    // Grid dimensions 
    int                 nx = 21, ny = 31, i, j;

    // Build a known grid
    point_2D **grid = allocate_2D_point_2D_array("grid", nx, ny);
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            grid[i][j].x = (long double) i;
            grid[i][j].y = (long double) j;
        }
    }

    // Call the function under test
    point_2D **bounds = get_x_bounds_2D(nx, ny, grid);

    // Check eta=0 boundary (j = 0) 
    for (i = 0; i < nx; i++)
    {
        cr_assert_float_eq (bounds[0][i].x, (long double) i, 1e-15L, "bounds[0][%d].x mismatch", i);
        cr_assert_float_eq (bounds[0][i].y, (long double) 0, 1e-15L, "bounds[0][%d].y mismatch", i);
    }

    // Check eta=1 boundary (j = ny - 1)
    for (i = 0; i < nx; i++)
    {
        cr_assert_float_eq (bounds[1][i].x, (long double) i, 1e-15L, "bounds[1][%d].x mismatch", i);
        cr_assert_float_eq (bounds[1][i].y, (long double) (ny - 1), 1e-15L, "bounds[1][%d].y mismatch", i);
    }

    // Free memory
    free_2D_point_2D_array("grid",   nx, grid);
    free_2D_point_2D_array("bounds", 2,  bounds);
}



/*
   Test that get_y_bounds_2D correctly extracts the
   xi = 0 and xi = 1 boundaries from a manufactured
   grid

   Grid is filled with points where x = i, y = j,
   so expected bounds are:
     bounds[0][j] = {0, j}      (xi = 0, i = 0)
     bounds[1][j] = {20, j}     (xi = 1, i = nx-1)
*/
Test(geometry, get_y_bounds_2D)
{
    /* Grid dimensions */ 
    int                 nx = 21, ny = 31, i, j;

    // Manufacture grid 
    point_2D **grid = allocate_2D_point_2D_array("grid", nx, ny); 
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            grid[i][j].x = (long double) i; 
            grid[i][j].y = (long double) j;
        }
    }

    // Call function 
    point_2D **bounds = get_y_bounds_2D(nx, ny, grid);

    // Check xi = 0 boundary (i = 0)
    for (j = 0; j < ny; j++)
    {
        cr_assert_float_eq (bounds[0][j].x, (long double) 0, 1e-15L, "bounds[0][%d].x mismatch", j);
        cr_assert_float_eq (bounds[0][j].y, (long double) j, 1e-15L, "bounds[0][%d].y mismatch", j);
    }

    // Check xi = 1 boundary (i = nx - 1) 
    for (j = 0; j < ny; j++)
    {
        cr_assert_float_eq (bounds[1][j].x, (long double) (nx - 1), 1e-15L, "bounds[1][%d].x mismatch", j);
        cr_assert_float_eq (bounds[1][j].y, (long double) j, 1e-15L, "bounds[1][%d].y mismatch", j);
    }

    // Free memory 
    free_2D_point_2D_array ("grid", nx, grid);
    free_2D_point_2D_array ("bounds", 2, bounds);
}


// Tests for exp_clustering
// Check if endpoints lie at x = 0.0 and x = 1.0
Test(geometry, exp_clustering_endpoints)
{
    int                 n = 11, dir;
    long double         *pts; 

    for (dir = 1; dir <= 4; dir++)
    {
        // Call function 
        pts = exp_clustering (11, dir, 2.0L, 0.0L, 1.0L);

        // Assert that the first and last points are equal 
        // to min and max 
        cr_assert_float_eq (pts[0], 0.0L, 1e-15L, "Exponential clustering failed: direction %d, pts[0] != 0.0", dir);
        cr_assert_float_eq (pts[10], 1.0L, 1e-15L, "Exponential clustering failed: direction %d, pts[n - 1] != 1.0", dir);

        // Free memory
        free_1D_long_double_array ("clustered_pts", pts);
    }
}

// Check monotonicity of clustered points
Test(geometry, exp_clustering_monotonicity)
{
    int                 dir, i; 
    long double         *pts;

    for (dir = 1; dir <= 4; dir++)
    {
        pts = exp_clustering (11, dir, 2.0L, 0.0L, 1.0L);

        // Check monotonicity
        for (i = 0; i < 10; i++)
        {
            cr_assert (pts[i + 1] > pts[i], "Exponential clustering failed: direction %d, pts[%d] >= pts[%d]", dir, i + 1, i);
        }

        // Free memory 
        free_1D_long_double_array ("clustered_pts", pts);
    }
}

// Check symmetry of direction 1 and 2
Test(geometry, exp_clustering_symmetry)
{
    int                 i; 
    long double         *pts1, *pts2; 

    // Generate directions 1 and 2 
    pts1 = exp_clustering (11, 1, 2.0L, 0.0L, 1.0L);
    pts2 = exp_clustering (11, 2, 2.0L, 0.0L, 1.0L);

    // Check symmetry 
    for (i = 0; i < 11; i++)
    {
        cr_assert_float_eq (pts1[i] + pts2[10 - i], 1.0L, 1e-10L, "Exponential clustering failed: symmetry violated at i = %d", i);
    }

    // Free memory 
    free_1D_long_double_array ("clustered_pts 1", pts1);
    free_1D_long_double_array ("clustered_pts 2", pts2);
}


// Tests for tanh clustering 
// Check if endpoints lie at x = 0.0 and x = 1.0
Test(geometry, tanh_clustering_endpoints)
{
    int                 dir;
    long double         *pts; 

    for (dir = 1; dir <= 5; dir++)
    {
        // Call function
        pts = tanh_clustering (11, dir, 2.0L, 0.0L, 1.0L);

        // Assert that the first and last points are equal 
        // to min and max
        cr_assert_float_eq (pts[0], 0.0L, 1e-15L, "tanh clustering failed: direction %d, pts[0] != 0.0", dir);
        cr_assert_float_eq (pts[10], 1.0L, 1e-15L, "tanh clustering failed: direction %d, pts[10] != 1.0", dir);

        // Free memory 
        free_1D_long_double_array ("clustered_pts", pts);
    }
}

// Check monotonicity of clustered points 
Test(geometry, tanh_clustering_monotonicity)
{
    int                 dir, i;
    long double         *pts; 

    for (dir = 1; dir <= 5; dir++)
    {
        // Call function
        pts = tanh_clustering (11, dir, 2.0L, 0.0L, 1.0L);

        // Check monotonicity
        for (i = 0; i < 10; i++)
        {
            cr_assert (pts[i + 1] > pts[i], "tanh clustering failed: direction %d, pts[%d] <= pts[%d]", dir, i + 1, i);
        }

        // Free memory 
        free_1D_long_double_array ("clustered_pts", pts);
    }
}

// Check uniformity of clustered points 
Test(geometry, tanh_clustering_uniform)
{
    int                 i;
    long double         *pts; 

    // Call function
    pts = tanh_clustering (11, 5, 2.0L, 0.0L, 1.0L);

    // Check for uniform spacing
    for (i = 0; i < 11; i++)
    {
        cr_assert_float_eq (pts[i], ((long double) i/10.0L), 1e-10L, 
                            "tanh clustering failed: uniformity violated at i = %d", i);
    }

    // Free memory 
    free_1D_long_double_array ("clustered_pts", pts);
}

// Check symmetry of direction 1 and 2
Test(geometry, tanh_clustering_symmetry)
{
    int                 i; 
    long double         *pts1, *pts2; 

    // Generate directions 1 and 2 
    pts1 = tanh_clustering (11, 1, 2.0L, 0.0L, 1.0L);
    pts2 = tanh_clustering (11, 2, 2.0L, 0.0L, 1.0L);

    // Check symmetry 
    for (i = 0; i < 11; i++)
    {
        cr_assert_float_eq (pts1[i] + pts2[10 - i], 1.0L, 1e-10L, "tanh clustering failed: symmetry violated at i = %d", i);
    }

    // Free memory 
    free_1D_long_double_array ("clustered_pts 1", pts1);
    free_1D_long_double_array ("clustered_pts 2", pts2);
}
