#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "globvar.h"
#include "config.h"
#include "aux_functions.h"
#include "file_operations.h"
#include "geometry.h"
#include "grid_solvers.h"
#include "errors.h"



/*
   Main driver function of the grid generator
*/
int main (int argc, char *argv[])
{
    /* Run configuration parsed from input file */
    run_args                cfg;
    solver_params           params;

    /* Local variables */
    int                     i, j;
    grid_2D                 init_grid, final_grid;
    grid_der_2D             **dgrid_q;
    coeffs_1                **coeffs_q;
    long double             **J_field;
    grid_quality_2D         q;


    /* Require a single config-file argument */
    if (argc != 2)
    {
        fprintf (stderr, "usage: %s <case.in>\n", argv[0]);
        return 1;
    }


    /* Parse config file and apply defaults */
    read_config_file (argv[1], &cfg);
    apply_solver_defaults (&cfg);

    /* If an initial grid file was provided, read its
       dimensions from the PLOT3D header into cfg */
    if (strlen (cfg.init_grid_file) > 0)
    {
        cfg.nblocks = read_nblocks_from_plot3D (cfg.init_grid_file);
        cfg.npts_x  = read_npts_x_from_plot3D (cfg.mode, cfg.init_grid_file, cfg.nblocks);
        cfg.npts_y  = read_npts_y_from_plot3D (cfg.mode, cfg.init_grid_file, cfg.nblocks);
    }


    /* Set initial grid parameters */
    init_grid.nx = cfg.npts_x[0];
    init_grid.ny = cfg.npts_y[0];

    if (strlen (cfg.init_grid_file) > 0)    /* Initial grid provided */
    {
        init_grid.grid     = read_2D_singleblock_plot3D (cfg.init_grid_file);
        init_grid.x_bounds = get_x_bounds_2D ((int) init_grid.nx, (int) init_grid.ny, init_grid.grid);
        init_grid.y_bounds = get_y_bounds_2D ((int) init_grid.nx, (int) init_grid.ny, init_grid.grid);
        write_bounds_2D ("x_boundaries.dat", (int) init_grid.nx, init_grid.x_bounds);
        write_bounds_2D ("y_boundaries.dat", (int) init_grid.ny, init_grid.y_bounds);
    }
    else    /* No initial grid provided */
    {
        init_grid.x_bounds = read_x_bounds_2D ((int) init_grid.nx);
        init_grid.y_bounds = read_y_bounds_2D ((int) init_grid.ny);
        init_grid.grid     = algebraic_grid_2D ((int) init_grid.nx, (int) init_grid.ny,
                                                init_grid.x_bounds, init_grid.y_bounds);
    }
    write_2D_singleblock_vts ("initial_grid.vts", (int) init_grid.nx, (int) init_grid.ny, init_grid.grid, 1);


    /* Final grid parameters are the same as the initial grid */
    final_grid.nx       = init_grid.nx;
    final_grid.ny       = init_grid.ny;
    final_grid.x_bounds = init_grid.x_bounds;
    final_grid.y_bounds = init_grid.y_bounds;


    /* Pack runtime tunables into the solver-facing struct. Each
       solver uses only the subset applicable to it */
    params.omega = cfg.omega;
    params.lambda_max = cfg.lambda_max;
    params.ramp_iters = cfg.ramp_iters;
    params.inner_niter = cfg.inner_niter;
    params.tol_step = cfg.tol_step;
    params.tol_resid = cfg.tol_resid;

    /* Compute final grid by dispatching on the chosen solver */
    switch (cfg.solver)
    {
        case SOLVER_ELL:
            final_grid.grid = elliptic_grid_2D_point ((int) init_grid.nx, (int) init_grid.ny,
                                                      init_grid.grid, (int) cfg.niter, params);
            break;
        case SOLVER_ELL_NEWTON:
            final_grid.grid = elliptic_grid_2D ((int) init_grid.nx, (int) init_grid.ny,
                                                init_grid.grid, (int) cfg.niter, params);
            break;
        case SOLVER_PSN_TM:
            final_grid.grid = poisson_grid_2D_point ((int) init_grid.nx, (int) init_grid.ny,
                                                     init_grid.grid, (int) cfg.niter, params);
            break;
        case SOLVER_PSN_LIM:
            final_grid.grid = poisson_grid_2D_point_lim ((int) init_grid.nx, (int) init_grid.ny,
                                                         init_grid.grid, (int) cfg.niter, params);
            break;
        case SOLVER_NEWTON:
            final_grid.grid = poisson_grid_2D ((int) init_grid.nx, (int) init_grid.ny,
                                               init_grid.grid, (int) cfg.niter, params);
            break;
        default:
            fatal_error ("main: unhandled solver type in dispatch");
    }
    write_2D_singleblock_vts ("final_grid.vts", (int) final_grid.nx, (int) final_grid.ny, final_grid.grid, 0);


    /* Grid quality diagnostics */
    // Allocate arrays
    dgrid_q  = allocate_2D_grid_der_2D_array ("dgrid_q", (int) final_grid.nx, (int) final_grid.ny);
    coeffs_q = allocate_2D_coeffs_1_array ("coeffs_q", (int) final_grid.nx, (int) final_grid.ny);

    // Compute grid quality metrics
    grid_first_ders_2D ((int) final_grid.nx, (int) final_grid.ny, final_grid.grid, &dgrid_q);
    first_der_coefficients ((int) final_grid.nx, (int) final_grid.ny, dgrid_q, &coeffs_q);
    q = compute_grid_quality_2D ((int) final_grid.nx, (int) final_grid.ny, coeffs_q);

    // Append orthogonality field to final_grid.vts and close the file
    write_long_double_to_vts ("final_grid.vts", (int) final_grid.nx, (int) final_grid.ny,
                              "orthogonality", q.ortho, 0);

    // Append aspect ratio field to final_grid.vts and close the file
    write_long_double_to_vts ("final_grid.vts", (int) final_grid.nx, (int) final_grid.ny,
                              "aspect_ratio", q.aspect_ratio, 0);

    // Store Jacobian in a temporary long double array
    J_field = allocate_2D_long_double_array("J_field", (int) final_grid.nx, (int) final_grid.ny);
    for (i = 0; i < (int) final_grid.nx; i++)
    {
        for (j = 0; j < (int) final_grid.ny; j++)
        {
            J_field[i][j] = coeffs_q[i][j].J;
        }
    }

    // Append Jacobian field to final_grid.vts and close the file
    write_long_double_to_vts ("final_grid.vts", (int) final_grid.nx, (int) final_grid.ny,
                              "jacobian", J_field, 1);

    // Print diagnostics
    fprintf (stdout, "Orthogonality: max = %10.4Lf, mean = %10.4Lf, degenerate = %d\n",
             q.ortho_max, q.ortho_mean, q.n_degenerate);
    fprintf (stdout, "Aspect ratio: max = %10.4Lf, min = %10.4Lf\n",
             q.aspect_ratio_max, q.aspect_ratio_min);
    fprintf (stdout, "Jacobian: max = %10.4Lf, min = %10.4Lf, folded = %d\n",
             q.J_max, q.J_min, q.n_folded);

    // Print warnings
    if (q.ortho_max > 0.95L)
    {
        fprintf (stderr, "WARNING: Peak orthogonality %10.4Le exceeds 0.95 (near-folding)\n", q.ortho_max);
    }
    if (q.n_folded > 0)
    {
        fprintf (stderr, "WARNING: %d folded cells detected (J <= 0): grid is unusable\n", q.n_folded);
    }


    /* Free memory from grid quality arrays */
    free_2D_grid_der_2D_array ("dgrid_q", (int) final_grid.nx, dgrid_q);
    free_2D_coeffs_1_array ("coeffs_q", (int) final_grid.nx, coeffs_q);
    free_2D_long_double_array ("q.ortho", (int) final_grid.nx, q.ortho);
    free_2D_long_double_array ("J_field", (int) final_grid.nx, J_field);

    /* Free memory from final grid */
    free_2D_point_2D_array ("final_grid", final_grid.nx, final_grid.grid);

    /* Free memory from initial grid */
    free_2D_point_2D_array ("init_grid_x_bounds", 2, init_grid.x_bounds);
    free_2D_point_2D_array ("init_grid_y_bounds", 2, init_grid.y_bounds);
    free_2D_point_2D_array ("init_grid", cfg.npts_x[0], init_grid.grid);

    /* Release run config (frees npts_x/y/z) */
    free_run_args (&cfg);


    return 0;
}
