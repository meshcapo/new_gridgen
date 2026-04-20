#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aux_functions.h"
#include "errors.h"
#include "file_operations.h"
#include "types.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

#define _swap(x, y) { typeof(x) tmp = x; x = y; y = tmp; }


/*
   Dynamically allocate memory to a 1D int array

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
*/
int *allocate_1D_int_array (char *array_name, int nrows)
{
    /* Return allocated array */
    int                     *array;

    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(int));
    if (array == NULL) array_allocation_msg (array_name, 0);


    return array;
}






/*
   Dynamically allocate memory to a 1D long int array

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
*/
long *allocate_1D_long_array (char *array_name, int nrows)
{
    /* Return allocated array */
    long                     *array;

    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(long));
    if (array == NULL) array_allocation_msg (array_name, 0);


    return array;
}






/*
   Dynamically allocate memory to a 1D long double array

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
*/
long double *allocate_1D_long_double_array (char *array_name, int nrows)
{
    /* Return allocated array */
    long double             *array;

    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(long double));
    if (array == NULL) array_allocation_msg (array_name, 0);


    return array;
}






/*
   Dynamically allocate memory to a 1D point_2D array

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
*/
point_2D *allocate_1D_point_2D_array (char *array_name, int nrows)
{
    /* Return allocated array */
    point_2D                *array;

    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(point_2D));
    if (array == NULL) array_allocation_msg (array_name, 0);


    return array;
}






/*
   Dynamically allocate memory to a 1D
   bh_sol_2D array

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
*/
bh_sol_2D *allocate_1D_bh_sol_2D_array (char * array_name, int nrows)
{
    /* Return allocated array */
    bh_sol_2D               *array;

    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(bh_sol_2D));
    if (array == NULL) array_allocation_msg (array_name, 0);

    /* Initialize to zero */
    for (int i = 0; i < nrows; i++)
    {
        array[i].x                      = (long double) 0.0;
        array[i].y                      = (long double) 0.0;
        array[i].p                      = (long double) 0.0;
        array[i].q                      = (long double) 0.0;
    }


    return array;
}






/*
   Free memory from a dynamically allocated 1D int array

   Input parameters: array_name - array descriptor
                     array      - the array being deallocated
*/
void free_1D_int_array (char *array_name, int *array)
{
    if (array == NULL)
    {
        printf ("Unable to free memory from array %s.\n", array_name);
        return;
    }
    else
    {
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 1D long int array

   Input parameters: array_name - array descriptor
                     array      - the array being deallocated
*/
void free_1D_long_array (char *array_name, long *array)
{
    if (array == NULL)
    {
        printf ("Unable to free memory from array %s.\n", array_name);
        return;
    }
    else
    {
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 1D long double array

   Input parameters: array_name - array descriptor
                     array      - the array being deallocated
*/
void free_1D_long_double_array (char *array_name, long double *array)
{
    if (array == NULL)
    {
        printf ("Unable to free memory from array %s.\n", array_name);
        return;
    }
    else
    {
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 1D point_2D int array

   Input parameters: array_name - array descriptor
                     array      - the array being deallocated
*/
void free_1D_point_2D_array (char *array_name, point_2D *array)
{
    if (array == NULL)
    {
        printf ("Unable to free memory from array %s.\n", array_name);
        return;
    }
    else
    {
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 1D
   bh_sol_2D array

   Input parameters: array_name - array descriptor
                     array      - the array being deallocated
*/
void free_1D_bh_sol_2D_array (char *array_name, bh_sol_2D *array)
{
    if (array == NULL)
    {
        printf ("Unable to free memory from array %s.\n", array_name);
        return;
    }
    else
    {
        free (array);
    }
}






/*
   Dynamically allocate memory to a 2D array of type
   int

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
int **allocate_2D_int_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    int                     **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(int *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(int));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Dynamically allocate memory to a 2D array of type
   long double

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
long double **allocate_2D_long_double_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    long double             **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(long double *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(long double));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Dynamically allocate memory to a 2D array of type
   point_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
point_2D **allocate_2D_point_2D_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    point_2D                **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(point_2D *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(point_2D));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Dynamically allocate memory to a 2D array of type
   grid_der_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
grid_der_2D **allocate_2D_grid_der_2D_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    grid_der_2D             **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(grid_der_2D *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(grid_der_2D));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Dynamically allocate memory to a 2D array of type
   coeffs_1

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
coeffs_1 **allocate_2D_coeffs_1_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    coeffs_1                **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(coeffs_1 *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(coeffs_1));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Dynamically allocate memory to a 2D array of type
   grid_dder_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
grid_dder_2D **allocate_2D_grid_dder_2D_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    grid_dder_2D            **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(grid_dder_2D *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(grid_dder_2D));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Dynamically allocate memory to a 2D array of type
   bh_jacobian_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
ell_jacobian_2D **allocate_2D_ell_jacobian_2D_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    ell_jacobian_2D          **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(ell_jacobian_2D *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(ell_jacobian_2D));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Dynamically allocate memory to a 2D array of type
   bh_jacobian_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows
                     ncols      - number of columns
*/
bh_jacobian_2D **allocate_2D_bh_jacobian_2D_array (char *array_name, int nrows, int ncols)
{
    /* Return allocated array */
    bh_jacobian_2D          **array;

    /* Local variables */
    int                     i;


    /* Initialization */
    array                               = NULL;

    /* Allocation of memory */
    array                               = malloc (nrows * sizeof(bh_jacobian_2D *));
    if (array == NULL) array_allocation_msg (array_name, 0);

    for (i = 0; i < nrows; i++)
    {
        array[i]                        = malloc (ncols * sizeof(bh_jacobian_2D));
        if (array[i] == NULL) array_allocation_msg (array_name, i);
    }


    return array;
}






/*
   Free memory from a dynamically allocated 2D
   array of type int

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_int_array (char *array_name, int nrows, int **array)
{
    /* Local variables */
    int                     i;

    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 2D
   array of type long_double

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_long_double_array (char *array_name, int nrows, long double **array)
{
    /* Local variables */
    int                     i;

    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 2D
   array of type point_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_point_2D_array (char *array_name, int nrows, point_2D **array)
{
    /* Local variables */
    int                     i;

    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 2D
   array of type grid_der_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_grid_der_2D_array (char *array_name, int nrows, grid_der_2D **array)
{
    /* Local variables */
    int                     i;


    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 2D
   array of type coeffs_1

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_coeffs_1_array (char *array_name, int nrows, coeffs_1 **array)
{
    /* Local variables */
    int                     i;


    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 2D
   array of type grid_dder_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_grid_dder_2D_array (char *array_name, int nrows, grid_dder_2D **array)
{
    /* Local variables */
    int                     i;


    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 2D
   array of type ell_jacobian_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_ell_jacobian_2D_array (char *array_name, int nrows, ell_jacobian_2D **array)
{
    /* Local variables */
    int                     i;


    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Free memory from a dynamically allocated 2D
   array of type bh_jacobian_2D

   Input parameters: array_name - array descriptor
                     nrows      - number of rows in 2D array
                     array      - the array being deallocated
*/
void free_2D_bh_jacobian_2D_array (char *array_name, int nrows, bh_jacobian_2D **array)
{
    /* Local variables */
    int                     i;


    if (array == NULL)
    {
        printf("Trying to deallocate NULL array %s.\n", array_name);
        return;
    } else
    {
        for (i = 0; i < nrows; i++) {
            free (array[i]);
        }
        free (array);
    }
}






/*
   Read the number of x coordinates in each block
   when provided as command line arguments for 2D
   grids

   Input parameters: arg_list   - all command line arguments
                     nblks      - number of blocks in structured grid
*/
long *get_npts_x_arg_2D (char *arg_list[], long nblks)
{
    /* Return nx */
    long                    *nx;

    /* Local variables */
    int                     i;
    char                    *ptr;


    /* Allocate memory to array */
    nx                                  = allocate_1D_long_array ("npts_x", (int)nblks);

    for (i = 0; i < (int)nblks; i++)
    {
        nx[i]                           = strtol (arg_list[(2*(i + 2)) - 1], &ptr, 10);
    }


    return nx;
}






/*
   Read the number of y coordinates in each block
   when provided as command line arguments for 2D
   grids

   Input parameters: arg_list   - all command line arguments
                     nblks      - number of blocks in structured grid
*/
long *get_npts_y_arg_2D (char *arg_list[], long nblks)
{
    /* Return ny */
    long                    *ny;

    /* Local variables */
    int                     i;
    char                    *ptr;


    /* Allocate memory to array */
    ny                                  = allocate_1D_long_array ("npts_y", (int)nblks);

    for (i = 0; i < (int)nblks; i++)
    {
        ny[i]                           = strtol (arg_list[2*(i + 2)], &ptr, 10);
    }


    return ny;
}






/*
   Read the number of x coordinates in each block
   when provided as command line arguments for 3D
   grids

   Input parameters: arg_list   - all command line arguments
                     nblks      - number of blocks in structured grid
*/
long *get_npts_x_arg_3D (char *arg_list[], long nblks)
{
    /* Return nx */
    long                    *nx;

    /* Local variables */
    int                     i;
    char                    *ptr;


    /* Allocate memory to array */
    nx                                  = allocate_1D_long_array ("npts_x", (int)nblks);

    for (i = 0; i < (int)nblks; i++)
    {
        nx[i]                           = strtol (arg_list[3*(i + 1)], &ptr, 10);
    }

    return nx;
}






/*
   Read the number of y coordinates in each block
   when provided as command line arguments for 3D
   grids

   Input parameters: arg_list   - all command line arguments
                     nblks      - number of blocks in structured grid
*/
long *get_npts_y_arg_3D (char *arg_list[], long nblks)
{
    /* Return ny */
    long                    *ny;

    /* Local variables */
    int                     i;
    char                    *ptr;


    /* Allocate memory to array */
    ny                                  = allocate_1D_long_array ("npts_y", (int)nblks);

    for (i = 0; i < (int)nblks; i++)
    {
        ny[i]                           = strtol (arg_list[(3*(i + 1)) + 1], &ptr, 10);
    }

    return ny;
}






/*
   Read the number of z coordinates in each block
   when provided as command line arguments for 3D
   grids

   Input parameters: arg_list   - all command line arguments
                     nblks      - number of blocks in structured grid
*/
long *get_npts_z_arg_3D (char *arg_list[], long nblks)
{
    /* Return nz */
    long                    *nz;

    /* Local variables */
    int                     i;
    char                    *ptr;


    /* Allocate memory to array */
    nz                                  = allocate_1D_long_array ("npts_z", (int)nblks);

    for (i = 0; i < (int)nblks; i++)
    {
        nz[i]                           = strtol (arg_list[(3*(i + 1)) + 2], &ptr, 10);
    }

    return nz;
}






/*
   Function to read command line arguments into global
   variables

   Input parameters: narg       - number of command line arguments
                     arg_list   - all command line arguments
                     run_mode   - grid dimensionality (2D/3D)
                     nx         - number of points in x
                     ny         - number of points in y
                     nz         - number of points in z
                     filename   - name of the file containing the initial grid, if any
*/
void read_inputs (int narg, char *arg_list[], char **run_mode, long *nblks,
                  long **nx, long **ny, long **nz, char **filename, long *niter)
{
    /* Local variables */
    char                    *ptr_nblks, *ext, *ptr_niter;


    /* Raise a fatal error if there are no command
       line arguments */
    if (narg == 1)
    {
        fatal_error ("No inputs specified.");
    }
    else
    {
        /* Get the grid dimensionality as a string */
        *run_mode                       = arg_list[1];

        /* Get the number of blocks */
        *nblks                          = strtol (arg_list[2], &ptr_nblks, 10);

        /* Read necessary inputs in 2D mode */
        if (strcmp(*run_mode, "2D") == 0)
        {
            /* Get the number of (x, y) points for each block
               from the command line arguments if an initial
               grid hasn't been provided */
            if (strcmp(ptr_nblks, "") == 0)
            {
                *nx                     = get_npts_x_arg_2D (arg_list, *nblks);
                *ny                     = get_npts_y_arg_2D (arg_list, *nblks);
                *nz                     = NULL;
            }
            else    // Initial grid provided
            {
                *filename               = arg_list[2];

                /* Determine extension of initial grid file */
                ext                     = get_file_extension (*filename);
                file_extension_error (ext);

                /* Read inputs from a plot3D file */
                if (strcmp(ext, ".x") == 0)
                {
                    *nblks              = read_nblocks_from_plot3D (*filename);
                    *nx                 = read_npts_x_from_plot3D (*run_mode, *filename, *nblks);
                    *ny                 = read_npts_y_from_plot3D (*run_mode, *filename, *nblks);
                    *nz                 = NULL;
                }
            }
        }

        /* Read necessary inputs in 3D mode */
        else if (strcmp(*run_mode, "3D") == 0)
        {
            /* Get the number of (x,y,z) points for each block
               from command line arguments if an initial grid
               hasn't been provided */
            if (strcmp(ptr_nblks, "") == 0)
            {
                *nx                     = get_npts_x_arg_3D (arg_list, *nblks);
                *ny                     = get_npts_y_arg_3D (arg_list, *nblks);
                *nz                     = get_npts_z_arg_3D (arg_list, *nblks);
            }
            if (strcmp(ptr_nblks, "") != 0)
            {
                *filename                   = arg_list[2];

                /* Determine extension of initial grid file */
                ext                     = get_file_extension (*filename);
                file_extension_error (ext);

                /* Read inputs from a plot3D file */
                if (strcmp(ext, ".x") == 0)
                {
                    *nblks              = read_nblocks_from_plot3D (*filename);
                    *nx                 = read_npts_x_from_plot3D (*run_mode, *filename, *nblks);
                    *ny                 = read_npts_y_from_plot3D (*run_mode, *filename, *nblks);
                    *nz                 = read_npts_z_from_plot3D (*filename, *nblks);
                }
            }
        }

        /* Raise fatal error if code isn't run in 2D or 3D mode */
        else
        {
            fatal_error ("Wrong mode provided.");
        }

        /* Get the number of iterations */
        *niter                          = strtol (arg_list[narg - 1], &ptr_niter, 10);
    }
}






/*
   Hash functions for converting 2D array
   indices to 1D array index

   Input parameters: i  - 1st index of 2D array
                     j  - 2nd index of 2D array
                     nx - number of rows in 2D array
*/
int hash (int i, int j, int nx)
{
    /* Return key value */
    int                     key;


    /* Key value */
    key                                 = i + ((j - 1) * (nx - 2)) - 1;
    //key                                 = i + (nx * j);


    return key;
}






/*
   Convert a 2D point_2D array to a
   1D array

   Input parameters: nx     - number of rows in 2D array
                     ny     - number of columns in 2D array
                     array  - the 2D array to be converted
*/
point_2D *convert_array_to_column (int nx, int ny, point_2D **array)
{
    /* Return col */
    point_2D                *col;

    /* Local variables */
    int                     i, j, k;


    /* Allocate memory */
    col                                 = allocate_1D_point_2D_array ("column",
                                          (nx - 2) * (ny - 2));

    /* Store 2D array elements into a
       1D column array */
    for (j = 1; j < ny - 1; j++)
    {
        for (i = 1; i < nx - 1; i++)
        {
            k                           = hash (i, j, nx);
            col[k].x                    = array[i][j].x;
            col[k].y                    = array[i][j].y;
        }
    }


    return col;
}






/*
   Convert a 1D point_2D array to a
   2D array

   Input parameter: nx      - number of rows in 2D array
                    ny      - number of columns in 2D array
                    col     - the 1D array to be converted
                    array   - resulting 2D array
*/
void convert_column_to_array (int nx, int ny, point_2D *col, point_2D ***array)
{
    /* Local variables */
    int                     i, j, k;


    for (j = 1; j < ny - 1; j++)
    {
        for (i = 1; i < nx - 1; i++)
        {
            k                           = hash (i, j, nx);
            (*array)[i][j].x            = col[k].x;
            (*array)[i][j].y            = col[k].y;
        }
    }
}






/***************************************************************************/
/*                                                                         */
/*              OPERATOR OVERLOADING FOR USER DEFINED DATA TYPES           */
/*                                                                         */
/***************************************************************************/



/*                          SET VARIABLES AS ZERO                          */

/* Set point_2D variable to zero */
point_2D zero_point_2D (void)
{
    /* Return a */
    point_2D                    a;


    /* a = 0 */
    a.x                                 = 0.0L;
    a.y                                 = 0.0L;


    return a;
}



/* Set grid_der_2D variable to zero */
grid_der_2D zero_grid_der_2D (void)
{
    /* Return a */
    grid_der_2D                 a;


    /* a = 0 */
    a.x_der                             = zero (point_2D);
    a.y_der                             = zero (point_2D);


    return a;
}



/* Set coeffs_1 variable to zero */
coeffs_1 zero_coeffs_1 (void)
{
    /* Return a */
    coeffs_1                    a;


    /* a = 0 */
    a.alpha                             = 0.0L;
    a.beta                              = 0.0L;
    a.gamma                             = 0.0L;
    a.J                                 = 0.0L;


    return a;
}



/* Set grid_dder_2D variable to zero */
grid_dder_2D zero_grid_dder_2D (void)
{
    /* Return a */
    grid_dder_2D                a;


    /* a = 0 */
    a.x_der                             = zero (point_2D);
    a.y_der                             = zero (point_2D);
    a.xy_der                            = zero (point_2D);


    return a;
}



/* Set bh_sol_2D variable to zero */
bh_sol_2D zero_bh_sol_2D (void)
{
    /* Return a */
    bh_sol_2D                   a;


    /* a = 0 */
    a.x                                 = 0.0L;
    a.y                                 = 0.0L;
    a.p                                 = 0.0L;
    a.q                                 = 0.0L;


    return a;
}



/* Set ell_jacobian_2D variable to zero */
ell_jacobian_2D zero_ell_jacobian_2D (void)
{
    /* Return a */
    ell_jacobian_2D             a;


    /* a = 0 */
    a.ddx                               = zero (point_2D);
    a.ddy                               = zero (point_2D);


    return a;
}



/* Set bh_jacobian_2D variable to zero */
bh_jacobian_2D zero_bh_jacobian_2D (void)
{
    /* Return a */
    bh_jacobian_2D             a;


    /* a = 0 */
    a.ddx                               = zero (bh_sol_2D);
    a.ddy                               = zero (bh_sol_2D);
    a.ddp                               = zero (bh_sol_2D);
    a.ddq                               = zero (bh_sol_2D);

    return a;
}



/* --------------------- CHECK FOR VARIABLE EQUALITY --------------------- */

/* Check equality of two point_2D variables */
int equal_point_2D (point_2D a, point_2D b)
{
    return (fabsl (a.x - b.x) < 1e-15L) &&
           (fabsl (a.y - b.y) < 1e-15L);
}



/* Check equality of two grid_der_2D variables */
int equal_grid_der_2D (grid_der_2D a, grid_der_2D b)
{
    return (equal (a.x_der, b.x_der) &&
            equal (a.y_der, b.y_der));
}



/* Check equality of two coeffs_1 variables */
int equal_coeffs_1 (coeffs_1 a, coeffs_1 b)
{
    return (fabsl (a.alpha - b.alpha) < 1e-15L) &&
           (fabsl (a.beta - b.beta) < 1e-15L) &&
           (fabsl (a.gamma - b.gamma) < 1e-15L) &&
           (fabsl (a.J - b.J) < 1e-15L);
}



/* Check equality of two grid_dder_2D variables */
int equal_grid_dder_2D (grid_dder_2D a, grid_dder_2D b)
{
    return (equal (a.x_der, b.x_der)) &&
           (equal (a.y_der, b.y_der)) &&
           (equal (a.xy_der, b.xy_der));
}



/* Check equality of two bh_sol_2D variables */
int equal_bh_sol_2D(bh_sol_2D a, bh_sol_2D b)
{
    return (fabsl (a.x - b.x) < 1e-15L) &&
           (fabsl (a.y - b.y) < 1e-15L) &&
           (fabsl (a.p - b.p) < 1e-15L) &&
           (fabsl (a.q - b.q) < 1e-15L);
}



/* Check equality of two ell_jacobian_2D variables */
int equal_ell_jacobian_2D (ell_jacobian_2D a, ell_jacobian_2D b)
{
    return (equal (a.ddx, b.ddx)) &&
           (equal (a.ddy, b.ddy));
}



/* Check equality of two bh_jacobian_2D variables */
int equal_bh_jacobian_2D (bh_jacobian_2D a, bh_jacobian_2D b)
{
    return (equal (a.ddx, b.ddx)) &&
           (equal (a.ddy, b.ddy)) &&
           (equal (a.ddp, b.ddp)) &&
           (equal (a.ddq, b.ddq));
}



/* Set point_2D variable equal to another */
point_2D set_as_equal_point_2D (point_2D a)
{
    /* Return b */
    point_2D                b;


    /* b = a */
    b.x                                 = a.x;
    b.y                                 = a.y;


    return b;
}



/* Multiply a point_2D variable with -1 */
void negative_point_2D (point_2D *a)
{
    (*a).x                              *= (long double) -1.0;
    (*a).y                              *= (long double) -1.0;
}



/* Divide a point_2D variable by a number */
point_2D divide_point_2D (point_2D a, long double d)
{
    /* Return b */
    point_2D                b;


    /* Check if d == 0 */ 
    if (fabsl(d) < 1e-15L) fatal_error ("Error in divide_point_2D: Attempting to divide by zero");

    /* b = a/d */
    b.x                                 = a.x/d;
    b.y                                 = a.y/d;

    return b;
}



/* In a ell_jacobian_2D variable, set df(x)/dx
   and df(y)/dy components as one */
ell_jacobian_2D identity_ders_ell ()
{
    /* return ders */
    ell_jacobian_2D         ders;

    ders.ddx.x                          = (long double) 1.0;
    ders.ddx.y                          = (long double) 0.0;
    ders.ddy.x                          = (long double) 0.0;
    ders.ddy.y                          = (long double) 1.0;

    return ders;
}



/*
   Check if a variable of type
   bh_jacobian_2D is zero or not.
   If zero, return 1 otherwise return 0
*/
int is_zero (bh_jacobian_2D ders)
{
    /* Return is_zero */
    int                     is_zero = 1;

    /* Local variables */
    bh_jacobian_2D          local_ders;


    /* Populate local variable with zeros */
    local_ders = zero (bh_jacobian_2D);

    /* Check if incoming variable is zero or not */ 
    is_zero = equal (ders, local_ders);
    

    return is_zero;
}






/*
   Print a bh_jacobian_2D object as a
   4x4 2Daarray of type long double
*/
void print_as_2D_matrix (int n, bh_jacobian_2D **ders)
{
    /* Local variables */
    int                     k1, k2;


    /* Print as 2D matrix */
    for (k1 = 0; k1 < n; k1++)
    {
        for (k2 = 0; k2 < n - 1; k2++)
        {
            printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF    ",
                    ders[k1][k2].ddx.x, ders[k1][k2].ddy.x, ders[k1][k2].ddp.x, ders[k1][k2].ddq.x);
        }
        printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF",
                ders[k1][n - 1].ddx.x, ders[k1][n - 1].ddy.x, ders[k1][n - 1].ddp.x, ders[k1][n - 1].ddq.x);
        printf ("\n");
        for (k2 = 0; k2 < n - 1; k2++)
        {
            printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF    ",
                    ders[k1][k2].ddx.y, ders[k1][k2].ddy.y, ders[k1][k2].ddp.y, ders[k1][k2].ddq.y);
        }
        printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF",
                ders[k1][n - 1].ddx.y, ders[k1][n - 1].ddy.y, ders[k1][n - 1].ddp.y, ders[k1][n - 1].ddq.y);
        printf ("\n");
        for (k2 = 0; k2 < n - 1; k2++)
        {
            printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF    ",
                    ders[k1][k2].ddx.p, ders[k1][k2].ddy.p, ders[k1][k2].ddp.p, ders[k1][k2].ddq.p);
        }
        printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF",
                ders[k1][n - 1].ddx.p, ders[k1][n - 1].ddy.p, ders[k1][n - 1].ddp.p, ders[k1][n - 1].ddq.p);
        printf ("\n");
        for (k2 = 0; k2 < n - 1; k2++)
        {
            printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF    ",
                    ders[k1][k2].ddx.q, ders[k1][k2].ddy.q, ders[k1][k2].ddp.q, ders[k1][k2].ddq.q);
        }
        printf ("%10.6LF    %10.6LF    %10.6LF    %10.6LF",
                ders[k1][n - 1].ddx.q, ders[k1][n - 1].ddy.q, ders[k1][n - 1].ddp.q, ders[k1][n - 1].ddq.q);
    }
}






/***************************************************************************/
/*                                                                         */
/*                          ARRAY ARITHMETIC FUNCTIONS                     */
/*                                                                         */
/***************************************************************************/



/*
   Allocate memory to a 1D point_2D array
   and initialize it to zero

   Input parameters: array_name - array descriptor string
                     nrows      - number of rows
*/
point_2D *zero_1D_point_2D_array (char *array_name, int nrows)
{
    /* Return array */
    point_2D                *array;

    /* Local variables */
    int                     i;


    /* Allocate memory */
    array                               = allocate_1D_point_2D_array (array_name, nrows);

    /* Set all elements to 0 */
    for (i = 0; i < nrows; i++)
    {
        array[i]                        = zero(point_2D);
    }


    return array;
}



/*
   Allocate memory to a 2D point_2D array
   and initialize it to zero

   Input parameters: array_name - array descriptor string
                     nrows      - number of rows
                     ncols      - number of columns
*/
point_2D **zero_2D_point_2D_array (char *array_name, int nrows, int ncols)
{
    /* Return array */
    point_2D                **array;

    /* Local variables */
    int                     i, j;


    /* Allocate memory */
    array                               = allocate_2D_point_2D_array (array_name, nrows, ncols);

    /* Set all elements to 0.0 */
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            array[i][j]                 = zero(point_2D);
        }
    }


    return array;
}



/*
   Set 1D array_2 equal to array_2

   Input parameters: nelem   - number of elements in incoming array
                     array_1 - array to be copied
*/
point_2D *equals_1D_point_2D_array (int nelem, point_2D *array_1)
{
    /* Return array_2 */
    point_2D                *array_2;


    /* Set array_2 = array_1 */
    array_2                             = allocate_1D_point_2D_array ("array_2", nelem);
    for (int i = 0; i < nelem; i++)
    {
        array_2[i]                      = set_as_equal_point_2D (array_1[i]);
    }


    return array_2;
}



/*
   Set 2D array_2 equal to array_1

   Input parameters: nrows      - number of rows
                     ncols      - number of columns
*/
void equals_2D_point_2D_array (int nrows, int ncols, point_2D **array_1, point_2D ***array_2)
{
    /* Local variables */
    int                     i, j;


    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            (*array_2)[i][j]            = set_as_equal_point_2D (array_1[i][j]);
        }
    }
}



/*
   Element-wise subtraction between two
   point_2D 2D arrays

   Input parameters: nrows      - number of rows
                     ncols      - number of columns
*/
point_2D **sub_2D_point_2D_arrays (int nrows, int ncols, point_2D **array_1, point_2D **array_2)
{
    /* Return array */
    point_2D                **array;

    /* Local variables */
    int                     i, j;


    /* Initialization */
    array                               = allocate_2D_point_2D_array ("sub_array", nrows, ncols);

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            array[i][j].x               = array_2[i][j].x - array_1[i][j].x;
            array[i][j].y               = array_2[i][j].y - array_1[i][j].y;
        }
    }


    return array;
}



/*
   Find maximum element in a 2D long double
   array

   Input parameters: nrows  - number of rows
                     ncols  - number of columns
                     array  - 2D array
*/
long double array_max_2D_long_double (int nrows, int ncols, long double **array)
{
    /* Return max */
    long double             max;

    /* Local variables */
    int                     i, j;


    /* Initialize max */
    max                                 = array[0][0];

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            /* Find maximum x */
            if (array[i][j] > max)
                max                     = array[i][j];
        }
    }


    return max;
}



/*
   Find maximum element in a 2D point_2D
   array. The maximums for both (x, y)
   componenets are computed individually

   Input parameters: nrows  - number of rows
                     ncols  - number of columns
                     array  - incoming 2D point_2D array
*/
point_2D array_max_2D_point_2D (int nrows, int ncols, point_2D **array)
{
    /* Return max */
    point_2D                max;

    /* Local variables */
    int                     i, j;


    /* Initialize max */
    max                                 = array[0][0];

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            /* Find maximum x */
            if (array[i][j].x > max.x)
                max.x                   = array[i][j].x;

            /* Find maximum y */
            if (array[i][j].y > max.y)
                max.y                   = array[i][j].y;
        }
    }


    return max;
}



/*
   Allocate memory to a 2D long double array
   and initialize it to zero

   Input parameters: array_name - array descriptor string
                     nrows      - number of rows
                     ncols      - number of columns
*/
long double **zero_2D_long_double_array (char *array_name, int nrows, int ncols)
{
    /* Return array */
    long double             **array;

    /* Local variables */
    int                     i, j;


    /* Allocate memory */
    array                               = allocate_2D_long_double_array (array_name, nrows, ncols);

    /* Set all elements to 0.0 */
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            array[i][j]                 = (long double)0.0;
        }
    }


    return array;
}



/*
   Set all of the elements of an existing
   2D long double array to zero

   Input parameters: nrows - number of rows
                     ncols - number of columns
                     array - array to be set to zero
*/
void equals_zero_2D_long_double_array (int nrows, int ncols, long double ***array)
{
    /* Local variables */
    int                     i, j;

    /* Set all elements to 0 */
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            (*array)[i][j]              = (long double) 0.0;
        }
    }
}



/*
   Set all elements of an existing
   1D point_2D array to zero

   Input parameters: nelem - number of array elements
                     array - array to be set to zero
*/
void equals_zero_1D_point_2D_array (int nelem, point_2D **array)
{
    /* Local variables */
    int                     i;

    /* Set all elements to 0 */
    for (i = 0; i < nelem; i++)
    {
        (*array)[i].x                   = (long double) 0.0;
        (*array)[i].y                   = (long double) 0.0;
    }
}






/*
   Convert a 2D array of type ell_jacobian_2D
   to an array of type long double. Used for
   solving linear system using LU decomposition

   TODO: Replace with a direct LU decomp
         routine using point_2D and
         ell_jacobian_2D
*/
void convert_ell_jacobian_2D (int n, ell_jacobian_2D **mat, long double ***array)
{
    for (int k1 = 0; k1 < n; k1++)
    {
        int l1                          = 2 * k1;
        for (int k2 = 0; k2 < n; k2++)
        {
            int l2                      = 2 * k2;
            (*array)[l1][l2]            = mat[k1][k2].ddx.x;
            (*array)[l1][l2 + 1]        = mat[k1][k2].ddy.x;
            (*array)[l1 + 1][l2]        = mat[k1][k2].ddx.y;
            (*array)[l1 + 1][l2 + 1]    = mat[k1][k2].ddy.y;
        }
    }
}






/*
   Convert a 2D array of type ell_jacobian_2D
   to an array of type gsl_spmatrix. Used for
   solving linear system using GMRES
*/
void convert_to_gsl_spmatrix_2D (int n, ell_jacobian_2D **mat, gsl_spmatrix **array)
{
    for (int k1 = 0; k1 < n; k1++)
    {
        int l1                          = 2 * k1;
        for (int k2 = 0; k2 < n; k2++)
        {
            int l2                      = 2 * k2;
            gsl_spmatrix_set (*array, l1, l2, (double) mat[k1][k2].ddx.x);
            gsl_spmatrix_set (*array, l1, l2 + 1, (double) mat[k1][k2].ddy.x);
            gsl_spmatrix_set (*array, l1 + 1, l2, (double) mat[k1][k2].ddx.y);
            gsl_spmatrix_set (*array, l1 + 1, l2 + 1, (double) mat[k1][k2].ddy.y);
        }
    }
}






/*
   Convert an array of type point_2D
   to an array of type long double. Used
   for solving linear system using LU
   decomposition

   TODO: Replace with a direct LU decomp
         routine using point_2D and
         ell_jacobian_2D
*/
void convert_point_2D (int n, point_2D *vec, long double **array)
{
    for (int k = 0; k < n; k++)
    {
        int l                           = 2 * k;
        (*array)[l]                     = vec[k].x;
        (*array)[l + 1]                 = vec[k].y;
    }
}






/*
   Convert an array of type long double
   to an array of type point_2D. Used for
   solving linear system using LU decomp

   TODO: Replace with a direct LU decomp
         routine using point_2D and
         ell_jacobian_2D
*/
void convert_to_point_2D (int n, long double *array, point_2D **vec)
{
    for (int k = 0; k < n; k++)
    {
        int l                           = 2 * k;
        (*vec)[k].x                     = array[l];
        (*vec)[k].y                     = array[l + 1];
    }
}






/*
   Convert an array of type point_2D
   to an array of type gsl_vector. Used
   for solving linear system using GMRES
*/
void convert_to_gsl_vector_2D (int n, point_2D *vec, gsl_vector **array)
{
    for (int k = 0; k < n; k++)
    {
        int l                           = 2 * k;
        gsl_vector_set (*array, l, (double) vec[k].x);
        gsl_vector_set (*array, l + 1, (double) vec[k].y);
    }
}






/*
   Convert an array of type gsl_vector to
   an array of type point_2D. Used when
   solving linear system using GMRES
*/
void convert_gsl_vector_2D (int n, gsl_vector *array, point_2D **vec)
{
    for (int k = 0; k < n; k++)
    {
        int l                           = 2 * k;
        (*vec)[k].x                     = (long double) gsl_vector_get (array, l);
        (*vec)[k].y                     = (long double) gsl_vector_get (array, l + 1);
    }
}






/*
   Convert a 2D array of type bh_jacobian_2D
   to an array of type long double. Used for
   solving linear system using LU decomposition

   TODO: Replace with a direct LU decomp
         routine using bh_sol_2D and
         bh_jacobian_2D
*/
void convert_bh_jacobian_2D (int n, bh_jacobian_2D **mat, long double ***array)
{
    for (int k1 = 0; k1 < n; k1++)
    {
        int l1                          = 4 * k1;
        for (int k2 = 0; k2 < n; k2++)
        {
            int l2                      = 4 * k2;
            (*array)[l1][l2]            = mat[k1][k2].ddx.x;
            (*array)[l1][l2 + 1]        = mat[k1][k2].ddy.x;
            (*array)[l1][l2 + 2]        = mat[k1][k2].ddp.x;
            (*array)[l1][l2 + 3]        = mat[k1][k2].ddq.x;
            (*array)[l1 + 1][l2]        = mat[k1][k2].ddx.y;
            (*array)[l1 + 1][l2 + 1]    = mat[k1][k2].ddy.y;
            (*array)[l1 + 1][l2 + 2]    = mat[k1][k2].ddp.y;
            (*array)[l1 + 1][l2 + 3]    = mat[k1][k2].ddq.y;
            (*array)[l1 + 2][l2]        = mat[k1][k2].ddx.p;
            (*array)[l1 + 2][l2 + 1]    = mat[k1][k2].ddy.p;
            (*array)[l1 + 2][l2 + 2]    = mat[k1][k2].ddp.p;
            (*array)[l1 + 2][l2 + 3]    = mat[k1][k2].ddq.p;
            (*array)[l1 + 3][l2]        = mat[k1][k2].ddx.q;
            (*array)[l1 + 3][l2 + 1]    = mat[k1][k2].ddy.q;
            (*array)[l1 + 3][l2 + 2]    = mat[k1][k2].ddp.q;
            (*array)[l1 + 3][l2 + 3]    = mat[k1][k2].ddq.q;
        }
    }
}






/*
   Convert an array of type bh_sol_2D
   to an array of type long double. Used
   for solving linear system using LU
   decomposition.

   TODO: Replace with a direct LU decomp
         routine using bh_sol_2D and
         bh_jacobian_2D
*/
void convert_bh_sol_2D (int n, bh_sol_2D *vec, long double **array)
{
    for (int k = 0; k < n; k++)
    {
        int l                           = 4 * k;
        (*array)[l]                     = vec[k].x;
        (*array)[l + 1]                 = vec[k].y;
        (*array)[l + 2]                 = vec[k].p;
        (*array)[l + 3]                 = vec[k].q;
    }
}






/*
   Convert an array of type long double
   to an array of type bh_sol_2D. Used
   for solving linear system using LU
   decomposition

   TODO: Replace with a direct LU decomp
         routine using bh_sol_2D and
         bh_jacobian_2D
*/
void convert_to_bh_sol_2D (int n, long double *array, bh_sol_2D **vec)
{
    for (int k = 0; k < n; k++)
    {
        int l                           = 4 * k;
        (*vec)[k].x                     = array[l];
        (*vec)[k].y                     = array[l + 1];
        (*vec)[k].p                     = array[l + 2];
        (*vec)[k].q                     = array[l + 3];
    }
}






/*
   Compute LU factorization of a square
   matrix M

   Input parameters: n - number of rows and columns in matrix
*/
void std_matrix_LU (int n, long double **M, long double ***L, long double ***U)
{
    /* Local variables */
    int                     i, j, k;
    long double             sum;


    /* Set L and U to zero */
    equals_zero_2D_long_double_array (n, n, L);
    equals_zero_2D_long_double_array (n, n, U);

    /* Compute LU factorization */
    for (i = 0; i < n; i++)
    {
        for (k = i; k < n; k++)
        {
            sum                         = 0;
            for (j = 0; j < i; j++)
            {
                sum                     += ((*L)[i][j] * (*U)[j][k]);
            }
            (*U)[i][k]                  = M[i][k] - sum;
        }

        for (k = i; k < n; k++)
        {
            if (i == k)
                (*L)[i][i]              = 1.0;
            else {
                sum                     = 0;
                for (j = 0; j < i; j++)
                {
                    sum                 += ((*L)[k][j] * (*U)[j][i]);
                }
                (*L)[k][i]              = (M[k][i] - sum)/(*U)[i][i];
            }
        }
    }
}






/*
   Solve a linear system Mx = b using LU
   decomposition of matrix M. The function
   first solves a system Lx_mid = b using
   forward substitution and then solves
   the system Ux = x_mid using backward
   substitution

   Input parameters: n - order of the linear system
                     M - nxn matrix
                     b - rhs vector
*/
void LU_linear_system_solve (int n, long double **M, long double *b, long double **x)
{
    /* Local variables */
    int                     i, j;
    long double             **L, **U, *x_mid, sum;


    /* Get LU factorization of incoming matrix M */
    L                                   = allocate_2D_long_double_array ("LUdecomp_L", n, n);
    U                                   = allocate_2D_long_double_array ("LUdecomp_U", n, n);
    std_matrix_LU (n, M, &L, &U);

    /* Allocate memory for x_mid */
    x_mid                               = allocate_1D_long_double_array ("x_mid", n);

    /* Solve the system Lx_mid = b using forward
       substitution */
    x_mid[0]                            = b[0];
    for (i = 1; i < n; i++)
    {
        sum                             = 0.0;
        for (j = 0; j <= i - 1; j++)
        {
            sum                         += L[i][j] * x_mid[j];
        }
        x_mid[i]                        = b[i] - sum;
    }

    /* Now solve the system Ux = x_mid using
       backward substitution */
    (*x)[n - 1]                         = x_mid[n - 1]/U[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        sum                             = 0.0;
        for (j = i + 1; j <= n - 1; j++)
        {
            sum                         += U[i][j] * (*x)[j];
        }
        (*x)[i]                         = (x_mid[i] - sum)/U[i][i];
    }


    /* Free memory */
    free_2D_long_double_array ("LUdecomp_L", n, L);
    free_2D_long_double_array ("LUdecomp_U", n, U);
    free_1D_long_double_array ("x_mid", x_mid);
}






/*
   Solve the linear system Mx = b created when
   using Newton's method for grid generation
   with GMRES

   TODO: Currently set up for elliptic grid
         generation only

   Input parameters: debug  - flag for printing optional output
                     nx     - number of x points
                     ny     - number of y points
                     dFdu   - 2D array of type ell_jacobian_2D containing
                              the Jacobian derivatives which constitute the
                              RHS matrix M
                     fu     - 1D array of type point_2D containing the
                              formulations of F(u) which constitute the
                              LHS vector b
*/
void gmres_linear_system_solve (int debug, int nx, int ny, ell_jacobian_2D **dFdu, point_2D *fu,
                                point_2D **del_u)
{
    /* Local variables */
    const int                           dim = 2 * nx * ny;
    gsl_spmatrix                        *M, *C;
    gsl_vector                          *b, *x;
    const double                        tol = 10E-6;
    const size_t                        max_iter = 100;
    const gsl_splinalg_itersolve_type   *T;
    gsl_splinalg_itersolve              *work;
    size_t                              iter = 0;
    double                              residual;
    int                                 status;


    /* Allocate GSL matrices */
    M                                   = gsl_spmatrix_alloc (dim, dim);

    /* Allocate GSL vectors */
    b                                   = gsl_vector_alloc (dim);
    x                                   = gsl_vector_alloc (dim);


    /* Convert incoming ell_jacobian_2D lhs
       array to a sparse matrix needed by GSL */
    convert_to_gsl_spmatrix_2D (nx * ny, dFdu, &M);

    /* Convert incoming point_2D rhs array to
       a vector needed by GSL GMRES routine */
    convert_to_gsl_vector_2D (nx * ny, fu, &b);

    /* Convert to compressed column format */
    C                                   = gsl_spmatrix_ccs (M);


    /* Set up for GMRES */
    T                                   = gsl_splinalg_itersolve_gmres;
    work                                = gsl_splinalg_itersolve_alloc (T, dim, 0);
    gsl_vector_set_zero (x);

    /* Solve the system Mx = b */
    do
    {
        status                          = gsl_splinalg_itersolve_iterate (C, b, tol, x, work);
        residual                        = gsl_splinalg_itersolve_normr (work);

        /* Print residual value if needed */
        //if (debug) printf ("%zu    %22.16LF\n", iter, (long double) residual);
    } while (status == GSL_CONTINUE && ++iter < max_iter);

    /* Print number of GMRES iterations needed for convergence */
    if (debug) printf ("GMRES iterations = %zu    %22.16LF\n", iter, (long double) residual);

    /* Convert solution to a point_2D array */
    convert_gsl_vector_2D (nx * ny, x, del_u);


    /* Free memory */
    gsl_splinalg_itersolve_free (work);
    gsl_spmatrix_free (M);
    gsl_spmatrix_free (C);
    gsl_vector_free (b);
    gsl_vector_free (x);
}
