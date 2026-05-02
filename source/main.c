#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "globvar.h"
#include "aux_functions.h"
#include "file_operations.h"
#include "geometry.h"
#include "grid_solvers.h"



/*
   Main driver function of the grid generator
*/
int main (int argc, char *argv[] )
{
    char                    *mode, *init_gridname;
    long                    nblocks, *npts_x, *npts_y, *npts_z, niter;

    /* Local variables */
    grid_2D                 init_grid, final_grid;


    /* Initialize grid generator inputs */
    mode                                = "";
    nblocks                             = 0;
    npts_x                              = NULL;
    npts_y                              = NULL;
    npts_z                              = NULL;
    init_gridname                       = "";
    niter                               = 0;

    /* Read command line arguments into global variables
       read_inputs is defined in aux_functions.c */
    read_inputs (argc, argv, &mode, &nblocks, &npts_x, &npts_y, &npts_z, &init_gridname, &niter);


    /* Set initial grid parameters */
    init_grid.nx                    = npts_x[0];
    init_grid.ny                    = npts_y[0];

    if (strcmp(init_gridname, "") != 0) /* Initial grid provided */
    {
        init_grid.grid                  = read_2D_singleblock_plot3D (init_gridname);
        init_grid.x_bounds              = get_x_bounds_2D ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid);
        init_grid.y_bounds              = get_y_bounds_2D ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid);
        write_bounds_2D ("x_boundaries.dat", (int)init_grid.nx, init_grid.x_bounds);
        write_bounds_2D ("y_boundaries.dat", (int)init_grid.ny, init_grid.y_bounds);
    }
    else    /* No initial grid provided */
    {
        //init_grid.x_bounds              = create_x_bounds_2D ((int)init_grid.nx);
        //init_grid.y_bounds              = create_y_bounds_2D ((int)init_grid.ny);
        //write_bounds_2D ("x_boundaries.dat", (int)init_grid.nx, init_grid.x_bounds);
        //write_bounds_2D ("y_boundaries.dat", (int)init_grid.ny, init_grid.y_bounds);
        init_grid.x_bounds              = read_x_bounds_2D ((int)init_grid.nx);
        init_grid.y_bounds              = read_y_bounds_2D ((int)init_grid.ny);
        init_grid.grid                  = algebraic_grid_2D ((int)init_grid.nx, (int)init_grid.ny, init_grid.x_bounds,
                                                             init_grid.y_bounds);
    }
    //write_2D_singleblock_plot3D ("initial_grid.x", (int) init_grid.nx, (int) init_grid.ny, init_grid.grid);
    write_2D_singleblock_vts ("initial_grid.vts", (int) init_grid.nx, (int) init_grid.ny, init_grid.grid, 1);

    /* Final grid parameters are the same as the initial grid */
    final_grid.nx                       = init_grid.nx;
    final_grid.ny                       = init_grid.ny;
    final_grid.x_bounds                 = init_grid.x_bounds;
    final_grid.y_bounds                 = init_grid.y_bounds;

    /* Compute final biharmonic grid */
    //final_grid.grid                 = biharmonic_grid_2D ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid,
    //                                                      (int)niter);
    //final_grid.grid                     = elliptic_grid_2D ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid,
    //                                                        (int)niter);
    //final_grid.grid                     = elliptic_grid_2D_point ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid,
    //                                                              (int)niter);
    //final_grid.grid                     = poisson_grid_2D ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid,
    //                                                        (int)niter);
    //final_grid.grid                     = poisson_grid_2D_point ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid,
    //                                                             (int)niter);
    final_grid.grid                     = poisson_grid_2D_point_lim ((int)init_grid.nx, (int)init_grid.ny, init_grid.grid,
                                                                     (int)niter);
    write_2D_singleblock_plot3D ("final_grid.x", (int) init_grid.nx, (int) init_grid.ny, final_grid.grid);


    /* Free memory from final grid */
    free_2D_point_2D_array ("final_grid", final_grid.nx, final_grid.grid);

    /* Free memory from initial grid */
    free_2D_point_2D_array ("init_grid_x_bounds", 2, init_grid.x_bounds);
    free_2D_point_2D_array ("init_grid_y_bounds", 2, init_grid.y_bounds);
    free_2D_point_2D_array ("init_grid", npts_x[0], init_grid.grid);

    /* Free memory from arrays containing number of points */
    free_1D_long_array ("npts_x", npts_x);
    free_1D_long_array ("npts_y", npts_y);
    if (strcmp(mode, "3D") == 0) free_1D_long_array ("npts_z", npts_z);


    return 0;
}
