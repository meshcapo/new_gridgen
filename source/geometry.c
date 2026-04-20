#include <stdio.h>
#include <math.h>
#include "types.h"
#include "aux_functions.h"



/*
   Get the points along the x boundaries (which map to
   xi = 0 and xi = 1) from a given 2D grid

   Input parameters: nx     - number of x points in the grid
                     ny     - number of y points in the grid
                     grid   - the 2D grid
*/
point_2D **get_x_bounds_2D (int nx, int ny, point_2D **grid)
{
    /* Return bounds */
    point_2D                **bounds;

    /* Local variables */
    int                     i;


    /* Initialization */
    bounds                              = allocate_2D_point_2D_array ("x_bounds", 2, nx);

    for (i = 0; i < nx; i++)
    {
        bounds[0][i]                    = set_as_equal_point_2D (grid[i][0]);
        bounds[1][i]                    = set_as_equal_point_2D (grid[i][ny - 1]);
    }


    return bounds;
}






/*
   Get the points along the y boundaries (which map to
   eta = 0 and eta = 1) from a given 2D grid

   Input parameters: nx     - number of x points in the grid
                     ny     - number of y points in the grid
                     grid   - the 2D grid
*/
point_2D **get_y_bounds_2D (int nx, int ny, point_2D **grid)
{
    /* Return bounds */
    point_2D                **bounds;

    /* Local variables */
    int                     j;


    /* Initialization */
    bounds                              = allocate_2D_point_2D_array ("y_bounds", 2, ny);

    for (j = 0; j < ny; j++)
    {
        bounds[0][j]                    = set_as_equal_point_2D (grid[0][j]);
        bounds[1][j]                    = set_as_equal_point_2D (grid[nx - 1][j]);
    }


    return bounds;
}






/*
   Use an exponential equation to cluster
   points given a minimum and maximum. Used
   for constructing either straight line
   geometries, the independent variable in a
   y = f(x) geometry or a parameter used in
   parametric definition of a curve.

   Input parameters: n          - number of points
                     direction  - index to determine where to cluster points
                     alpha      - clustering parameter
                     min        - minimum value
                     max        - maximum value
*/
long double *exp_clustering (int n, int direction, long double alpha, long double min, long double max)
{
    /* Return clustered_pts */
    long double             *clustered_pts;

    /* Local variables */
    const long double       HALF = 0.5, ONE = 1.0, TWO = 2.0;
    long double             *xi, xihalf[(n + 1)/2], temp3, temp4[(n + 1)/2], temp5[(n + 1)/2];
    int                     i;
    double                  temp1, temp2;


    /* Initialization */
    clustered_pts                       = allocate_1D_long_double_array ("clustered_pts", n);
    xi                                  = allocate_1D_long_double_array ("uniform_xi", n);


    /* Generate clustered points */
    switch (direction)
    {
        case 1: // cluster towards min

            for (i = 0; i < n; i++)
            {
                xi[i]                   = (long double) i/((long double) (n - 1));
                temp1                   = (double) (xi[i] * alpha);
                temp2                   = (double) alpha;
                temp3                   = ((long double) exp(temp1) - ONE)/((long double) exp(temp2) - ONE);
                clustered_pts[i]        = min + ((max - min) * temp3);
            }
            break;

        case 2: // cluster towards max

            for (i = 0; i < n; i++)
            {
                xi[i]                   = (long double) i/((long double) (n - 1));
                temp1                   = (double) (xi[i] * alpha);
                temp2                   = (double) alpha;
                temp3                   = ((long double) exp(temp1) - ONE)/((long double) exp(temp2) - ONE);
                clustered_pts[n - 1 - i]= min + ((max - min) * (ONE - temp3));
            }
            break;

        case 3: // cluster in the middle

            for (i = 0; i < (n + 1)/2; i++)
            {
                xihalf[i]               = (long double) i/((long double) ((n - 1)/2));
                temp1                   = (double) (xihalf[i] * alpha);
                temp2                   = (double) alpha;
                temp3                   = ((long double) exp(temp1) - ONE)/((long double) exp(temp2) - ONE);
                temp4[((n - 1)/2) - i]  = HALF * (ONE - temp3); // cluster towards 0.5
                temp5[i]                = HALF * (ONE + temp3); // cluster away from 0.5
            }

            /* Clustered points */
            for (i = 0; i < n; i++)
            {
                if (i < (n + 1)/2) {
                    clustered_pts[i]    = min + ((max - min) * temp4[i]);
                } else {
                    clustered_pts[i]    = min + ((max - min) * (temp5[i - ((n - 1)/2)]));
                }
            }
            break;

        case 4: // cluster towards both ends

            for (i = 0; i < (n + 1)/2; i++)
            {
                xihalf[i]               = (long double) i/((long double) ((n - 1)/2));
                temp1                   = (double) (xihalf[i] * alpha);
                temp2                   = (double) alpha;
                temp3                   = ((long double) exp(temp1) - ONE)/((long double) exp(temp2) - ONE);
                temp4[i]                = HALF * temp3; // cluster towards 0
                temp5[((n - 1)/2) - i]  = HALF * (TWO - temp3); // cluster towards 1
            }

            /* Clustered points */
            for (i = 0; i < n; i++)
            {
                if (i < (n + 1)/2) {
                    clustered_pts[i]    = min + ((max - min) * temp4[i]);
                } else {
                    clustered_pts[i]    = min + ((max - min) * (temp5[i - ((n - 1)/2)]));
                }
            }
            break;
    }


    /* Free memory */
    free_1D_long_double_array ("uniform_xi", xi);


    return clustered_pts;
}






/*
   Use a hyperolic tangent equation to cluster
   points given a minimum and maximum. Used
   for constructing either straight line
   geometries, the independent variable in a
   y = f(x) geometry or a parameter used in
   parametric definition of a curve.

   Input parameters: n          - number of points
                     direction  - index to determine where to cluster points
                     alpha      - clustering parameter
                     min        - minimum value
                     max        - maximum value
*/
long double *tanh_clustering (int n, int direction, long double alpha, long double min, long double max)
{
    /* Return clustered_pts */
    long double             *clustered_pts;

    /* Local variables */
    const int               nhalf = (n + 1)/2;
    const long double       HALF = 0.5, ONE = 1.0, TWO = 2.0;
    long double             *xi, temp3, xihalf[nhalf], temp4[nhalf], temp5[nhalf];
    double                  temp1, temp2;
    int                     i;


    /* Initialization */
    clustered_pts                       = allocate_1D_long_double_array ("clustered_pts", n);
    xi                                  = allocate_1D_long_double_array ("uniform_xi", n);


    /* Cluster points */
    switch (direction)
    {
        case 1: // cluster towards min

            for (i = 0; i < n; i++)
            {
                xi[i]                   = (long double) i/((long double) (n - 1));
                temp1                   = (double) (alpha * (xi[i] - ONE));
                temp2                   = (double) alpha;
                temp3                   = ONE + (((long double) tanh(temp1))/((long double) tanh(temp2)));
                clustered_pts[i]        = min + ((max - min) * temp3);
            }
            break;

        case 2: // cluster towards max

            for (i = 0; i < n; i++)
            {
                xi[i]                   = (long double) i/((long double) (n - 1));
                temp1                   = (double) (alpha * (xi[i] - ONE));
                temp2                   = (double) alpha;
                temp3                   = -(((long double) tanh(temp1))/((long double) tanh(temp2)));
                clustered_pts[n - i - 1]= min + ((max - min) * temp3);
            }
            break;

        case 3: // cluster towards middle

            for (i = 0; i < nhalf; i++)
            {
                xihalf[i]               = (long double) i/((long double) (nhalf - 1));
                temp1                   = (double) (alpha * (xihalf[i] - ONE));
                temp2                   = (double) alpha;
                temp3                   = (((long double) tanh(temp1))/((long double) tanh(temp2)));
                temp4[nhalf - i - 1]    = -HALF * temp3;
                temp5[i]                = HALF * (TWO + temp3);
            }

            /* Clustered points */
            for (i = 0; i < n; i++)
            {
                if (i < nhalf) {
                    clustered_pts[i]    = min + ((max - min) * temp4[i]);
                } else {
                    clustered_pts[i]    = min + ((max - min) * (temp5[i - nhalf + 1]));
                }
            }
            break;

        case 4: // cluster towards both ends

            for (i = 0; i < n; i++)
            {
                xi[i]                   = (long double) i/((long double) (n - 1));
                temp1                   = (double) (TWO * alpha * (xi[i] - HALF));
                temp2                   = (double) (alpha);
                temp3                   = ((long double) tanh(temp1))/((long double) tanh(temp2));
                clustered_pts[i]        = min + ((max - min) * (HALF * (ONE + temp3)));
            }
            break;

        case 5: // no clustering

            for (i = 0; i < n; i++)
            {
                xi[i]                   = (long double) i/((long double) (n - 1));
                clustered_pts[i]        = min + ((max - min) * xi[i]);
            }
    }


    /* Free memory */
    free_1D_long_double_array ("uniform_xi", xi);


    return clustered_pts;
}






/*
   Create boundaries mapping to eta = 0 and eta = 1

   Input parameters: n - number of points
*/
point_2D **create_x_bounds_2D (int n)
{
    /* Return x_bounds */
    point_2D                **x_bounds;

    /* Local variables */
    FILE                    *fptr1, *fptr2;


    /* Initialization */
    x_bounds                            = allocate_2D_point_2D_array ("x_bounds", 2, n);


    /* Read original blade section along with
       offset curve */
    fptr1                               = fopen ("blade.1.2.9S", "r");
    fptr2                               = fopen ("2D_inflated_hub_section.dat", "r");
    for (int i = 0; i < n; i++)
    {
        fscanf (fptr1, "%LF    %LF\n", &x_bounds[0][i].x, &x_bounds[0][i].y);
        fscanf (fptr2, "%LF    %LF\n", &x_bounds[1][i].x, &x_bounds[1][i].y);
    }
    fclose (fptr1);
    fclose (fptr2);


    return x_bounds;
}






/*
   Create boundaries mapping xi = 0 and xi = 1

   Input parameters: n - number of points
*/
point_2D **create_y_bounds_2D (int n)
{
    /* Return y_bounds */
    point_2D                **y_bounds;

    /* Local variables */
    point_2D                blade_TE, offset_TE;
    FILE                    *fptr1, *fptr2;
    long double             *pts, m, c;


    /* Initialization */
    y_bounds                            = allocate_2D_point_2D_array ("y_bounds", 2, n);


    /* Read TE points for original blade section
       and offset curve */
    fptr1                               = fopen ("blade.1.2.9S", "r");
    fptr2                               = fopen ("2D_inflated_hub_section.dat", "r");
    fscanf (fptr1, "%LF    %LF\n", &blade_TE.x, &blade_TE.y);
    fscanf (fptr2, "%LF    %LF\n", &offset_TE.x, &offset_TE.y);
    fclose (fptr1);
    fclose (fptr2);


    /* Cluster points between x coordinates of
       the TE points */
    pts                                 = tanh_clustering (n, 1, (long double) 1.15, blade_TE.x, offset_TE.x);

    /* Slope and intercept of the line segment
       connecting the two TE points */
    m                                   = (offset_TE.y - blade_TE.y)/(offset_TE.x - blade_TE.x);
    c                                   = blade_TE.y - (m * blade_TE.x);

    /* Generate the boundaries */
    for (int i = 0; i < n; i++)
    {
        y_bounds[0][i].x                = pts[i];
        y_bounds[0][i].y                = (m * y_bounds[0][i].x) + c;
        y_bounds[1][i]                  = set_as_equal_point_2D (y_bounds[0][i]);
    }


    /* Free memory */
    free_1D_long_double_array ("clustered_pts", pts);


    return y_bounds;
}
