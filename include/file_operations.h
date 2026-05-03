#ifndef FILE_OPERATIONS_H_
#define FILE_OPERATIONS_H_
#include "types.h"



/*
   Determine the extension of a file - used for
   reading the initial grid, if provided
*/
char *get_file_extension (char *);

/*
   Read number of blocks from a plot3D file
*/
long read_nblocks_from_plot3D (char *);

/*
   Read number of x, y and z points in each block
   of a structured grid from a plot3D file
*/
long *read_npts_x_from_plot3D (char *, char *, long);
long *read_npts_y_from_plot3D (char *, char *, long);
long *read_npts_z_from_plot3D (char *, long);

/*
   Read coordinates along the x boundaries (which
   map to eta = 0 and eta = 1) and y boundaries
   (which map to xi = 0 and xi = 1) of a 2D grid
*/
point_2D **read_x_bounds_2D (int);
point_2D **read_y_bounds_2D (int);

void write_bounds_2D (char *, int, point_2D **);

/*
   Read the normals to the x and y boundaries
   of a 2D grid
*/
point_2D **read_xbound_normals_2D (int);
point_2D **read_ybound_normals_2D (int);

/*
   Read a formatted 2D plot3D (.x) file containing a
   single block grid
*/
point_2D **read_2D_singleblock_plot3D (char *);

/*
   Write a single block structured 2D grid to a
   formatted 2D plot3D (.x)/XML structured VTK
   (.vts) file
*/
void write_2D_singleblock_plot3D (char *, int, int, point_2D **);
void write_2D_singleblock_vts (char *, int, int, point_2D **, int);

/*
   Append per-point scalar / vector field DataArrays to a .vts file
   in the open state (file already started by write_2D_singleblock_vts
   with write_data = 0). Pass write_data = 1 on the final call to
   close the file structure, or use close_vts explicitly
*/
void write_long_double_to_vts (char *, int, int, char *, long double **, int);
void write_point_2D_to_vts (char *, int, int, char *, point_2D **, int);
void write_grid_der_2D_to_vts (char *, int, int, grid_der_2D **, int);
void write_grid_dder_2D_to_vts (char *, int, int, grid_dder_2D **, int);

/*
   Close a .vts file (writes </PointData> + closing structure tags)
*/
void close_vts (char *);

/*
   Functions to write 2D arrays to .dat files
*/
void write_long_double_1D (char *, int, long double *);
void write_long_double_2D (char *, int, int, long double **);
void write_1D_point_2D (char *, int, point_2D *);
void write_2D_point_2D (char *, char *, int, int, point_2D **);
void write_grid_der_2D (char *, char *, char *, char *, int, int, grid_der_2D **);
void write_grid_dder_2D (int, int, grid_dder_2D **);
void write_x_to_file_2D (char *, int, int, point_2D **, point_2D **);
void write_dgrid_to_file_2D (char *, int, int, grid_der_2D **);
void write_coeffs_to_file_2D (char *, int, int, coeffs_1 **);
void write_d2grid_to_file_2D (char *, int, int, grid_dder_2D **);
void write_1D_bh_sol_2D (char *, int, int, bh_sol_2D *);



#endif // FILE_OPERATIONS_H_
