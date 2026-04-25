#include <stdio.h>
#include <math.h>
#include "globvar.h"
#include "types.h"
#include "aux_functions.h"
#include "file_operations.h"
#include "newton.h"


#define sqr(x) (x * x)

/* Global variables */
//const long double           ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0,
//                            FOUR = 4.0, FIVE = 5.0;



/*
   Create a 2D computational grid (xi, eta)
   which is a unit square

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     unit   - index to decide if we are in a unit square
                              or not
*/
point_2D **computational_grid (int nx, int ny, int unit)
{
    /* Return comp_grid */
    point_2D                **comp_grid;

    /* Local variables */
    int                     i, j;


    /* Initialization */
    comp_grid                           = allocate_2D_point_2D_array ("xieta", nx, ny);
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (unit == 0)
            {
                comp_grid[i][j].x       = ((long double) i)/((long double) (nx - 1));
                comp_grid[i][j].y       = ((long double) j)/((long double) (ny - 1));
            }
            else
            {
                comp_grid[i][j].x       = (long double) i;
                comp_grid[i][j].y       = (long double) j;
            }
        }
    }


    return comp_grid;
}






/*
   Compute an interpolated 2D algebraic grid given
   boundary points - used to construct an initial grid

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     x_bounds   - points on the boundaries mapping to
                                  xi = 0 and xi = 1
                     y_bounds   - points on the boundaries mapping to
                                  eta = 0 and eta = 1
*/
point_2D **algebraic_grid_2D (int nx, int ny, point_2D **x_bounds, point_2D **y_bounds)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    int                     i, j;
    point_2D                **xieta;
    long double             temp1, temp2, temp3, temp4;


    /* Initialization */
    grid                                = allocate_2D_point_2D_array ("algebraic_grid", nx, ny);


    /* Set boundary points */
    for (i = 0; i < nx; i++)
    {
        grid[i][0]                      = set_as_equal_point_2D (x_bounds[0][i]);
        grid[i][ny - 1]                 = set_as_equal_point_2D (x_bounds[1][i]);
    }
    for (j = 1; j < ny - 1; j++)
    {
        grid[0][j]                      = set_as_equal_point_2D (y_bounds[0][j]);
        grid[nx - 1][j]                 = set_as_equal_point_2D (y_bounds[1][j]);
    }


    /* Compute numerical grid */
    xieta                               = computational_grid (nx, ny, 0);

    /* Generate interpolated algebraic grid */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp1                       = ONE - xieta[i][j].x;  // 1 - xi
            temp2                       = ONE - xieta[i][j].y;  // 1 - eta

            temp3                       = (temp1 * grid[0][j].x) + (xieta[i][j].x * grid[nx - 1][j].x) +
                                          (temp2 * grid[i][0].x) + (xieta[i][j].y * grid[i][ny - 1].x);
            temp4                       = (temp1 * temp2 * grid[0][0].x) + (xieta[i][j].x * temp2 * grid[nx - 1][0].x) +
                                          (temp1 * xieta[i][j].y * grid[0][ny - 1].x) + (xieta[i][j].x * xieta[i][j].y *
                                           grid[nx - 1][ny - 1].x);
            grid[i][j].x                = temp3 - temp4;

            temp3                       = (temp1 * grid[0][j].y) + (xieta[i][j].x * grid[nx - 1][j].y) +
                                          (temp2 * grid[i][0].y) + (xieta[i][j].y * grid[i][ny - 1].y);
            temp4                       = (temp1 * temp2 * grid[0][0].y) + (xieta[i][j].x * temp2 * grid[nx - 1][0].y) +
                                          (temp1 * xieta[i][j].y * grid[0][ny - 1].y) + (xieta[i][j].x * xieta[i][j].y *
                                           grid[nx - 1][ny - 1].y);
            grid[i][j].y                = temp3 - temp4;
        }
    }


    /* Free memory from numerical grid array */
    free_2D_point_2D_array ("xieta", nx, xieta);


    return grid;
}






/*
   Compute the first derivatives of all (x, y) grid
   points with respect to the computational grid
   (xi, eta) using finite differences

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     grid   - (x, y) grid coordinates
                     dgrid  - array to store derivative
*/
void grid_first_ders_2D (int nx, int ny, point_2D **grid, grid_der_2D ***dgrid)
{
    /* Local variables */
    long double             dxi, deta;
    int                     i, j;


    /* Spacing of xi and eta in the uniform
       computational grid */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Interior: Use a O(h^2) central difference */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            /* Derivatives of (x, y) wrt xi */
            (*dgrid)[i][j].x_der.x      = (grid[i + 1][j].x - grid[i - 1][j].x)/(TWO * dxi);
            (*dgrid)[i][j].x_der.y      = (grid[i + 1][j].y - grid[i - 1][j].y)/(TWO * dxi);

            /* Derivatives of (x, y) wrt eta */
            (*dgrid)[i][j].y_der.x      = (grid[i][j + 1].x - grid[i][j - 1].x)/(TWO * deta);
            (*dgrid)[i][j].y_der.y      = (grid[i][j + 1].y - grid[i][j - 1].y)/(TWO * deta);
        }
    }

    /* eta = 0 and eta = 1 boundaries */
    for (i = 1; i < nx - 1; i++)
    {
        /* eta = 0 (j = 0): O(h^2) central differences
           for xi derivatives */
        (*dgrid)[i][0].x_der.x          = (grid[i + 1][0].x - grid[i - 1][0].x)/(TWO * dxi);
        (*dgrid)[i][0].x_der.y          = (grid[i + 1][0].y - grid[i - 1][0].y)/(TWO * dxi);

        /* O(h^2) forward differences for eta derivatives */
        (*dgrid)[i][0].y_der.x          = ((-THREE * grid[i][0].x) + (FOUR * grid[i][1].x) -
                                           grid[i][2].x)/(TWO * deta);
        (*dgrid)[i][0].y_der.y          = ((-THREE * grid[i][0].y) + (FOUR * grid[i][1].y) -
                                           grid[i][2].y)/(TWO * deta);


        /* eta = 1 (j = ny - 1): O(h^2) central differences
           for xi derivatives */
        (*dgrid)[i][ny - 1].x_der.x     = (grid[i + 1][ny - 1].x - grid[i - 1][ny - 1].x)/(TWO * dxi);
        (*dgrid)[i][ny - 1].x_der.y     = (grid[i + 1][ny - 1].y - grid[i - 1][ny - 1].y)/(TWO * dxi);

        /* O(h^2) backward differences for eta derivatives */
        (*dgrid)[i][ny - 1].y_der.x     = ((THREE * grid[i][ny - 1].x) - (FOUR * grid[i][ny - 2].x) +
                                           grid[i][ny - 3].x)/(TWO * deta);
        (*dgrid)[i][ny - 1].y_der.y     = ((THREE * grid[i][ny - 1].y) - (FOUR * grid[i][ny - 2].y) +
                                           grid[i][ny - 3].y)/(TWO * deta);
    }

    /* xi = 0 and xi = 1 boundaries */
    for (j = 1; j < ny - 1; j++)
    {
        /* xi = 0 (i = 0): O(h^2) forward differences
           for xi derivatives */
        (*dgrid)[0][j].x_der.x          = ((-THREE * grid[0][j].x) + (FOUR * grid[1][j].x) -
                                           grid[2][j].x)/(TWO * dxi);
        (*dgrid)[0][j].x_der.y          = ((-THREE * grid[0][j].y) + (FOUR * grid[1][j].y) -
                                           grid[2][j].y)/(TWO * dxi);

        /* O(h^2) central differences for eta derivatives */
        (*dgrid)[0][j].y_der.x          = (grid[0][j + 1].x - grid[0][j - 1].x)/(TWO * deta);
        (*dgrid)[0][j].y_der.y          = (grid[0][j + 1].y - grid[0][j - 1].y)/(TWO * deta);


        /* xi = 1 (i = nx - 1): O(h^2) backward differences
           for xi derivatives */
        (*dgrid)[nx - 1][j].x_der.x     = ((THREE * grid[nx - 1][j].x) - (FOUR * grid[nx - 2][j].x) +
                                           grid[nx - 3][j].x)/(TWO * dxi);
        (*dgrid)[nx - 1][j].x_der.y     = ((THREE * grid[nx - 1][j].y) - (FOUR * grid[nx - 2][j].y) +
                                           grid[nx - 3][j].y)/(TWO * dxi);

        /* O(h^2) central differences for eta derivatives */
        (*dgrid)[nx - 1][j].y_der.x     = (grid[nx - 1][j + 1].x - grid[nx - 1][j - 1].x)/(TWO * deta);
        (*dgrid)[nx - 1][j].y_der.y     = (grid[nx - 1][j + 1].y - grid[nx - 1][j - 1].y)/(TWO * deta);
    }


    /* xi = 0, eta = 0: O(h^2) forward differences for
       xi and eta derivatives */
    (*dgrid)[0][0].x_der.x              = ((-THREE * grid[0][0].x) + (FOUR * grid[1][0].x) -
                                           grid[2][0].x)/(TWO * dxi);
    (*dgrid)[0][0].x_der.y              = ((-THREE * grid[0][0].y) + (FOUR * grid[1][0].y) -
                                           grid[2][0].y)/(TWO * dxi);
    (*dgrid)[0][0].y_der.x              = ((-THREE * grid[0][0].x) + (FOUR * grid[0][1].x) -
                                           grid[0][2].x)/(TWO * deta);
    (*dgrid)[0][0].y_der.y              = ((-THREE * grid[0][0].y) + (FOUR * grid[0][1].y) -
                                           grid[0][2].y)/(TWO * deta);

    /* xi = 1, eta = 0: O(h^2) backward difference for
       xi derivatives and forward difference for eta */
    (*dgrid)[nx - 1][0].x_der.x         = ((THREE * grid[nx - 1][0].x) - (FOUR * grid[nx - 2][0].x) +
                                           grid[nx - 3][0].x)/(TWO * dxi);
    (*dgrid)[nx - 1][0].x_der.y         = ((THREE * grid[nx - 1][0].y) - (FOUR * grid[nx - 2][0].y) +
                                           grid[nx - 3][0].y)/(TWO * dxi);
    (*dgrid)[nx - 1][0].y_der.x         = ((-THREE * grid[nx - 1][0].x) + (FOUR * grid[nx - 1][1].x) -
                                           grid[nx - 1][2].x)/(TWO * deta);
    (*dgrid)[nx - 1][0].y_der.y         = ((-THREE * grid[nx - 1][0].y) + (FOUR * grid[nx - 1][1].y) -
                                           grid[nx - 1][2].y)/(TWO * deta);

    /* xi = 0, eta = 1: O(h^2) forward difference for xi
       derivatives and backward difference for eta */
    (*dgrid)[0][ny - 1].x_der.x         = ((-THREE * grid[0][ny - 1].x) + (FOUR * grid[1][ny - 1].x) -
                                           grid[2][ny - 1].x)/(TWO * dxi);
    (*dgrid)[0][ny - 1].x_der.y         = ((-THREE * grid[0][ny - 1].y) + (FOUR * grid[1][ny - 1].y) -
                                           grid[2][ny - 1].y)/(TWO * dxi);
    (*dgrid)[0][ny - 1].y_der.x         = ((THREE * grid[0][ny - 1].x) - (FOUR * grid[0][ny - 2].x) +
                                           grid[0][ny - 3].x)/(TWO * deta);
    (*dgrid)[0][ny - 1].y_der.y         = ((THREE * grid[0][ny - 1].y) - (FOUR * grid[0][ny - 2].y) +
                                           grid[0][ny - 3].y)/(TWO * deta);

    /* xi = 1, eta = 1: O(h^2) backward differences for
       xi and eta derivatives */
    (*dgrid)[nx - 1][ny - 1].x_der.x    = ((THREE * grid[nx - 1][ny - 1].x) - (FOUR * grid[nx - 2][ny - 1].x) +
                                           grid[nx - 3][ny - 1].x)/(TWO * dxi);
    (*dgrid)[nx - 1][ny - 1].x_der.y    = ((THREE * grid[nx - 1][ny - 1].y) - (FOUR * grid[nx - 2][ny - 1].y) +
                                           grid[nx - 3][ny - 1].y)/(TWO * dxi);
    (*dgrid)[nx - 1][ny - 1].y_der.x    = ((THREE * grid[nx - 1][ny - 1].x) - (FOUR * grid[nx - 1][ny - 2].x) +
                                           grid[nx - 1][ny - 3].x)/(TWO * deta);
    (*dgrid)[nx - 1][ny - 1].y_der.y    = ((THREE * grid[nx - 1][ny - 1].y) - (FOUR * grid[nx - 1][ny - 2].y) +
                                           grid[nx - 1][ny - 3].y)/(TWO * deta);
}






/*
   Compute new values of grid first derivatives
   at a point (i, j). This is done when using a
   Gauss-Seidel method where the values of some
   points with indices i - 1 or j - 1 have already
   been updated

   Input parameters: dxi        - xi step size
                     deta       - eta step size
                     grid_pts   - array containing coordinates of
                                  points surrounding (i, j) in the
                                  order: (i - 1, j - 1), (i, j - 1),
                                  (i + 1, j - 1), (i - 1, j), (i, j),
                                  (i + 1, j), (i - 1, j + 1), (i, j + 1),
                                  (i + 1, j + 1)
*/
grid_der_2D compute_first_ders_2D_at_ij (long double dxi, long double deta, point_2D *grid_pts)
{
    /* Return dij */
    grid_der_2D             dij;


    /* Compute new first ders */
    dij.x_der.x                         = (grid_pts[5].x - grid_pts[3].x)/(TWO * dxi);
    dij.x_der.y                         = (grid_pts[5].y - grid_pts[3].y)/(TWO * dxi);
    dij.y_der.x                         = (grid_pts[7].x - grid_pts[1].x)/(TWO * deta);
    dij.y_der.y                         = (grid_pts[7].y - grid_pts[1].y)/(TWO * deta);


    return dij;
}






/*
   Compute coefficients and the Jacobian needed
   for solving the biharmonic equation. These are
   all computed using only the first derivatives
   of (x, y) wrt (xi, eta) and are stored in an
   array of a data type defined for these quantities

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dgrid  - the (x, y) grid 1st derivatives
                     coeffs - array to store coefficients
*/
void first_der_coefficients (int nx, int ny, grid_der_2D **dgrid, coeffs_1 ***coeffs)
{
    /* Local variables */
    int                     i, k;


    /* Compute coefficients and Jacobian */
    for (i = 0; i < nx; i++)
    {
        for (k = 0; k < ny; k++)
        {
            (*coeffs)[i][k].alpha       = (dgrid[i][k].y_der.x * dgrid[i][k].y_der.x) +
                                          (dgrid[i][k].y_der.y * dgrid[i][k].y_der.y);
            (*coeffs)[i][k].beta        = (dgrid[i][k].x_der.x * dgrid[i][k].y_der.x) +
                                          (dgrid[i][k].x_der.y * dgrid[i][k].y_der.y);
            (*coeffs)[i][k].gamma       = (dgrid[i][k].x_der.x * dgrid[i][k].x_der.x) +
                                          (dgrid[i][k].x_der.y * dgrid[i][k].x_der.y);
            (*coeffs)[i][k].J           = (dgrid[i][k].x_der.x * dgrid[i][k].y_der.y) -
                                          (dgrid[i][k].x_der.y * dgrid[i][k].y_der.x);
        }
    }
}






/*
   Compute new values of first derivative coefficients
   at a point (i, j). This is done when using a
   Gauss-Seidel method where the values of some
   points with indices i - 1 or j - 1 have already
   been updated

   Input parameters: dxi        - xi step size
                     deta       - eta step size
                     grid_pts   - array containing coordinates of
                                  points surrounding (i, j) in the
                                  order: (i - 1, j - 1), (i, j - 1),
                                  (i + 1, j - 1), (i - 1, j), (i, j),
                                  (i + 1, j), (i - 1, j + 1), (i, j + 1),
                                  (i + 1, j + 1)
*/
coeffs_1 compute_coeffs_at_ij (grid_der_2D dij)
{
    /* Return coeff_ij */
    coeffs_1                coeff_ij;


    coeff_ij.alpha                      = (dij.y_der.x * dij.y_der.x) +
                                          (dij.y_der.y * dij.y_der.y);
    coeff_ij.beta                       = (dij.x_der.x * dij.y_der.x) +
                                          (dij.x_der.y * dij.y_der.y);
    coeff_ij.gamma                      = (dij.x_der.x * dij.x_der.x) +
                                          (dij.x_der.y * dij.x_der.y);
    coeff_ij.J                          = (dij.x_der.x * dij.y_der.y) -
                                          (dij.x_der.y * dij.y_der.x);

    return coeff_ij;
}






/*
   Function to compute coefficients and indices
   required for various O(h^2) finite difference
   approximations of a second derivative. The 1st
   row contains the coefficients, the 2nd row
   contains the i coefficients and the 3rd row
   contains the j coefficients

   Input parameter: ID  - identifier for the FD approximation
                    i   - i index of the current point
                    j   - j index of the current point
*/
int **finite_difference_2 (int ID, int i, int j)
{
    /* Return fd */
    int                     **fd;


    switch (ID)
    {
        case 1: // central difference for "x" 2nd der
            fd                          = allocate_2D_int_array ("fd", 3, 3);
            fd[0][0] = 1; fd[0][1] = -2; fd[0][2] = 1.0;
            fd[1][0] = i + 1; fd[1][1] = i; fd[1][2] = i - 1;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j;
            break;
        case 2: // forward difference for "x" 2nd der
            fd                          = allocate_2D_int_array ("fd", 3, 4);
            fd[0][0] = 2; fd[0][1] = -5; fd[0][2] = 4; fd[0][3] = -1;
            fd[1][0] = i; fd[1][1] = i + 1; fd[1][2] = i + 2; fd[1][3] = i + 3;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j; fd[2][3] = j;
            break;
        case 3: // backward difference for "x" 2nd der
            fd                          = allocate_2D_int_array ("fd", 3, 4);
            fd[0][0] = 2; fd[0][1] = -5; fd[0][2] = 4; fd[0][3] = -1;
            fd[1][0] = i; fd[1][1] = i - 1; fd[1][2] = i - 2; fd[1][3] = i - 3;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j; fd[2][3] = j;
            break;
        case 4: // central difference for "y" 2nd der
            fd                          = allocate_2D_int_array ("fd", 3, 3);
            fd[0][0] = 1; fd[0][1] = -2; fd[0][2] = 1.0;
            fd[1][0] = i; fd[1][1] = i; fd[1][2] = i;
            fd[2][0] = j + 1; fd[2][1] = j; fd[2][2] = j - 1;
            break;
        case 5: // forward difference for "y" 2nd der
            fd                          = allocate_2D_int_array ("fd", 3, 4);
            fd[0][0] = 2; fd[0][1] = -5; fd[0][2] = 4; fd[0][3] = -1;
            fd[1][0] = i; fd[1][1] = i; fd[1][2] = i; fd[1][3] = i;
            fd[2][0] = j; fd[2][1] = j + 1; fd[2][2] = j + 2; fd[2][3] = j + 3;
            break;
        case 6: // forward difference for "y" 2nd der
            fd                          = allocate_2D_int_array ("fd", 3, 4);
            fd[0][0] = 2; fd[0][1] = -5; fd[0][2] = 4; fd[0][3] = -1;
            fd[1][0] = i; fd[1][1] = i; fd[1][2] = i; fd[1][3] = i;
            fd[2][0] = j; fd[2][1] = j - 1; fd[2][2] = j - 2; fd[2][3] = j - 3;
            break;
        case 7: // central difference for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 4);
            fd[0][0] = 1; fd[0][1] = -1; fd[0][2] = -1; fd[0][3] = 1;
            fd[1][0] = i + 1; fd[1][1] = i + 1; fd[1][2] = i - 1; fd[1][3] = i - 1;
            fd[2][0] = j + 1; fd[2][1] = j - 1; fd[2][2] = j + 1; fd[2][3] = j - 1;
            break;
        case 8: // central difference in i, forward difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 6);
            fd[0][0] = -3; fd[0][1] = 3; fd[0][2] = 4; fd[0][3] = -4; fd[0][4] = -1; fd[0][5] = 1;
            fd[1][0] = i + 1; fd[1][1] = i - 1; fd[1][2] = i + 1; fd[1][3] = i - 1; fd[1][4] = i + 1; fd[1][5] = i - 1;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j + 1; fd[2][3] = j + 1; fd[2][4] = j + 2; fd[2][5] = j + 2;
            break;
        case 9: // central difference in i, backward difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 6);
            fd[0][0] = 3; fd[0][1] = -3; fd[0][2] = -4; fd[0][3] = 4; fd[0][4] = 1; fd[0][5] = -1;
            fd[1][0] = i + 1; fd[1][1] = i - 1; fd[1][2] = i + 1; fd[1][3] = i - 1; fd[1][4] = i + 1; fd[1][5] = i - 1;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j - 1; fd[2][3] = j - 1; fd[2][4] = j - 2; fd[2][5] = j - 2;
            break;
        case 10: // forward difference in i, central difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 6);
            fd[0][0] = -3; fd[0][1] = 3; fd[0][2] = 4; fd[0][3] = -4; fd[0][4] = -1; fd[0][5] = 1;
            fd[1][0] = i; fd[1][1] = i; fd[1][2] = i + 1; fd[1][3] = i + 1; fd[1][4] = i + 2; fd[1][5] = i + 2;
            fd[2][0] = j + 1; fd[2][1] = j - 1; fd[2][2] = j + 1; fd[2][3] = j - 1; fd[2][4] = j + 1; fd[2][5] = j - 1;
            break;
        case 11: // backward difference in i, central difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 6);
            fd[0][0] = 3; fd[0][1] = -3; fd[0][2] = -4; fd[0][3] = 4; fd[0][4] = 1; fd[0][5] = -1;
            fd[1][0] = i; fd[1][1] = i; fd[1][2] = i - 1; fd[1][3] = i - 1; fd[1][4] = i - 2; fd[1][5] = i - 2;
            fd[2][0] = j + 1; fd[2][1] = j - 1; fd[2][2] = j + 1; fd[2][3] = j - 1; fd[2][4] = j + 1; fd[2][5] = j - 1;
            break;
        case 12: // forward difference in i, forward difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 9);
            fd[0][0] = 9; fd[0][1] = -12; fd[0][2] = 3; fd[0][3] = -12; fd[0][4] = 16; fd[0][5] = -4;
            fd[0][6] = 3; fd[0][7] = -4; fd[0][8] = 1;
            fd[1][0] = i; fd[1][1] = i + 1; fd[1][2] = i + 2; fd[1][3] = i; fd[1][4] = i + 1; fd[1][5] = i + 2;
            fd[1][6] = i; fd[1][7] = i + 1; fd[1][8] = i + 2;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j; fd[2][3] = j + 1; fd[2][4] = j + 1; fd[2][5] = j + 1;
            fd[2][6] = j + 2; fd[2][7] = j + 2; fd[2][8] = j + 2;
            break;
        case 13: // backward difference in i, forward difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 9);
            fd[0][0] = -9; fd[0][1] = 12; fd[0][2] = -3; fd[0][3] = 12; fd[0][4] = -16; fd[0][5] = 4;
            fd[0][6] = -3; fd[0][7] = 4; fd[0][8] = -1;
            fd[1][0] = i; fd[1][1] = i - 1; fd[1][2] = i - 2; fd[1][3] = i; fd[1][4] = i - 1; fd[1][5] = i - 2;
            fd[1][6] = i; fd[1][7] = i - 1; fd[1][8] = i - 2;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j; fd[2][3] = j + 1; fd[2][4] = j + 1; fd[2][5] = j + 1;
            fd[2][6] = j + 2; fd[2][7] = j + 2; fd[2][8] = j + 2;
            break;
        case 14: // forward difference in i, backward difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 9);
            fd[0][0] = -9; fd[0][1] = 12; fd[0][2] = -3; fd[0][3] = 12; fd[0][4] = -16; fd[0][5] = 4;
            fd[0][6] = -3; fd[0][7] = 4; fd[0][8] = -1;
            fd[1][0] = i; fd[1][1] = i + 1; fd[1][2] = i + 2; fd[1][3] = i; fd[1][4] = i + 1; fd[1][5] = i + 2;
            fd[1][6] = i; fd[1][7] = i + 1; fd[1][8] = i + 2;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j; fd[2][3] = j - 1; fd[2][4] = j - 1; fd[2][5] = j - 1;
            fd[2][6] = j - 2; fd[2][7] = j - 2; fd[2][8] = j - 2;
            break;
        case 15: // backward difference in i, backward difference in j for mixed der
            fd                          = allocate_2D_int_array ("fd", 3, 9);
            fd[0][0] = 9; fd[0][1] = -12; fd[0][2] = 3; fd[0][3] = -12; fd[0][4] = 16; fd[0][5] = -4;
            fd[0][6] = 3; fd[0][7] = -4; fd[0][8] = 1;
            fd[1][0] = i; fd[1][1] = i - 1; fd[1][2] = i - 2; fd[1][3] = i; fd[1][4] = i - 1; fd[1][5] = i - 2;
            fd[1][6] = i; fd[1][7] = i - 1; fd[1][8] = i - 2;
            fd[2][0] = j; fd[2][1] = j; fd[2][2] = j; fd[2][3] = j - 1; fd[2][4] = j - 1; fd[2][5] = j - 1;
            fd[2][6] = j - 2; fd[2][7] = j - 2; fd[2][8] = j - 2;
            break;
    }


    return fd;

}






/*
   Compute a O(h^2) finite difference approximation
   of a second derivative (e.g.: d2x/dxi2) using either
   a central, forward, backward difference or a combination
   of the three

   Input parameters: ID     - identifier for FD type
                     i      - i index of current point
                     j      - j index of current point
                     array  - array containing all (x, y) points

   TODO: Definitely figure out how to not pass the entire
         grid with each function call
*/
point_2D compute_finite_diff_2 (int index, int i, int j, point_2D **array, int debug)
{
    /* Return der_2 */
    point_2D                der_2;

    /* Local variables */
    int                     **fd, n, k, l;


    /* Initialization */
    der_2                               = zero (point_2D);


    /* Compute coefficients and indices needed to
       compute the FD approximation */
    fd                                  = finite_difference_2 (index, i, j);

    /* Number of columns in the array fd */
    if (index == 1 || index == 4) {
        n                               = 3;
    } else if (index == 2 || index == 3 || index == 5 || index == 6 || index == 7) {
        n                               = 4;
    } else if (index == 8 || index == 9 || index == 10 || index == 11) {
        n                               = 6;
    } else {
        n                               = 9;
    }

    if (debug == 0)
    {
        printf ("Printing for debugging from compute_finite_diff_2\n");
        for (l = 0; l < 3; l++)
        {
            for (k = 0; k < n; k++)
            {
                printf ("%d    ", fd[l][k]);
            }
            printf ("\n");
        }
    }


    /* Compute the numerator of the FD approximation */
    for (k = 0; k < n; k++)
    {
        /* fd[0][k] - coefficient
           fd[1][k] - i index
           fd[2][k] - j index */
        der_2.x                         += ((long double) fd[0][k]) * (array[fd[1][k]][fd[2][k]].x);
        der_2.y                         += ((long double) fd[0][k]) * (array[fd[1][k]][fd[2][k]].y);
    }


    /* Free memory from fd */
    free_2D_int_array ("fd", 3, fd);


    return der_2;
}






/*
   Compute second and mixed derivatives of all (x, y)
   grid points with respect to the computational grid
   (xi, eta) using finite differences

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     grid   - (x, y) grid coordinates
                     d2grid - array containing grid derivatives
*/
void grid_second_ders_2D (int nx, int ny, point_2D **grid, grid_dder_2D ***d2grid)
{
    /* Local variables */
    int                     i, j;
    long double             dxi, deta;


    /* Spacing of xi and eta in the uniform
       computational grid */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            (*d2grid)[i][j].x_der       = zero (point_2D);
            (*d2grid)[i][j].xy_der      = zero (point_2D);
            (*d2grid)[i][j].y_der       = zero (point_2D);
        }
    }


    /* Interior: O(h^2) central differences */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            (*d2grid)[i][j].x_der       = compute_finite_diff_2 (1, i, j, grid, 1);
            (*d2grid)[i][j].xy_der      = compute_finite_diff_2 (7, i, j, grid, 1);
            (*d2grid)[i][j].y_der       = compute_finite_diff_2 (4, i, j, grid, 1);
        }
    }


    /* eta = 0 and eta = 1 boundaries */
    for (i = 1; i < nx - 1; i++)
    {
        /* eta = 0 (j = 0): O(h^2) central difference in xi
           and O(h^2) forward difference in eta */
        (*d2grid)[i][0].x_der           = compute_finite_diff_2 (1, i, 0, grid, 1);
        (*d2grid)[i][0].xy_der          = compute_finite_diff_2 (8, i, 0, grid, 1);
        (*d2grid)[i][0].y_der           = compute_finite_diff_2 (5, i, 0, grid, 1);

        /* eta = 1 (j = ny - 1): O(h^2) central difference
           in xi and O(h^2) backward difference in eta */
        (*d2grid)[i][ny - 1].x_der      = compute_finite_diff_2 (1, i, ny - 1, grid, 1);
        (*d2grid)[i][ny - 1].xy_der     = compute_finite_diff_2 (9, i, ny - 1, grid, 1);
        (*d2grid)[i][ny - 1].y_der      = compute_finite_diff_2 (6, i, ny - 1, grid, 1);
    }


    /* xi = 0 and xi = 1 boundaries */
    for (j = 1; j < ny - 1; j++)
    {
        /* xi = 0 (i = 0): O(h^2) forward difference in xi
           and O(h^2) central difference in eta */
        (*d2grid)[0][j].x_der           = compute_finite_diff_2 (2, 0, j, grid, 1);
        (*d2grid)[0][j].xy_der          = compute_finite_diff_2 (10, 0, j, grid, 1);
        (*d2grid)[0][j].y_der           = compute_finite_diff_2 (4, 0, j, grid, 1);

        /* xi = 1 (i = nx - 1): O(h^2) backward difference
           in xi and O(h^2) central difference in eta */
        (*d2grid)[nx - 1][j].x_der      = compute_finite_diff_2 (3, nx - 1, j, grid, 1);
        (*d2grid)[nx - 1][j].xy_der     = compute_finite_diff_2 (11, nx - 1, j, grid, 1);
        (*d2grid)[nx - 1][j].y_der      = compute_finite_diff_2 (4, nx - 1, j, grid, 1);
    }


    /* xi = 0, eta = 0: forward difference in xi
       and forward difference in eta */
    (*d2grid)[0][0].x_der               = compute_finite_diff_2 (2, 0, 0, grid, 1);
    (*d2grid)[0][0].xy_der              = compute_finite_diff_2 (12, 0, 0, grid, 1);
    (*d2grid)[0][0].y_der               = compute_finite_diff_2 (5, 0, 0, grid, 1);

    /* xi = 1, eta = 0: backward difference in xi
       and forward difference in eta */
    (*d2grid)[nx - 1][0].x_der          = compute_finite_diff_2 (3, nx - 1, 0, grid, 1);
    (*d2grid)[nx - 1][0].xy_der         = compute_finite_diff_2 (13, nx - 1, 0, grid, 1);
    (*d2grid)[nx - 1][0].y_der          = compute_finite_diff_2 (5, nx - 1, 0, grid, 1);

    /* xi = 0, eta = 1: forward difference in xi
       and backward difference in eta */
    (*d2grid)[0][ny - 1].x_der          = compute_finite_diff_2 (2, 0, ny - 1, grid, 1);
    (*d2grid)[0][ny - 1].xy_der         = compute_finite_diff_2 (14, 0, ny - 1, grid, 1);
    (*d2grid)[0][ny - 1].y_der          = compute_finite_diff_2 (6, 0, ny - 1, grid, 1);

    /* xi = 1, eta = 1: backward difference in xi
       and backward difference in eta */
    (*d2grid)[nx - 1][ny - 1].x_der     = compute_finite_diff_2 (3, nx - 1, ny - 1, grid, 1);
    (*d2grid)[nx - 1][ny - 1].xy_der    = compute_finite_diff_2 (15, nx - 1, ny - 1, grid, 1);
    (*d2grid)[nx - 1][ny - 1].y_der     = compute_finite_diff_2 (6, nx - 1, ny - 1, grid, 1);


    /* Complete FD approximation */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            (*d2grid)[i][j].x_der       = divide_point_2D ((*d2grid)[i][j].x_der, (dxi * dxi));
            (*d2grid)[i][j].xy_der      = divide_point_2D ((*d2grid)[i][j].xy_der, (FOUR * dxi * deta));
            (*d2grid)[i][j].y_der       = divide_point_2D ((*d2grid)[i][j].y_der, (deta * deta));
        }
    }
}






/*
   Compute new values of grid second derivatives
   at a point (i, j). This is done when using a
   Gauss-Seidel method where the values of some
   points with indices i - 1 or j - 1 have already
   been updated

   Input parameters: dxi        - xi step size
                     deta       - eta step size
                     grid_pts   - array containing coordinates of
                                  points surrounding (i, j) in the
                                  order: (i - 1, j - 1), (i, j - 1),
                                  (i + 1, j - 1), (i - 1, j), (i, j),
                                  (i + 1, j), (i - 1, j + 1), (i, j + 1),
                                  (i + 1, j + 1)
*/
grid_dder_2D compute_second_ders_2D_at_ij (long double dxi, long double deta, point_2D *grid_pts)
{
    /* Return d2ij */
    grid_dder_2D            d2ij;


    /* Compute new second ders */
    d2ij.x_der.x                        = (grid_pts[5].x - (TWO * grid_pts[4].x) +
                                           grid_pts[3].x)/(dxi * dxi);
    d2ij.x_der.y                        = (grid_pts[5].y - (TWO * grid_pts[4].y) +
                                           grid_pts[3].y)/(dxi * dxi);
    d2ij.xy_der.x                       = (grid_pts[8].x - grid_pts[2].x - grid_pts[6].x +
                                           grid_pts[0].x)/(FOUR * dxi * deta);
    d2ij.xy_der.y                       = (grid_pts[8].y - grid_pts[2].y - grid_pts[6].y +
                                           grid_pts[0].y)/(FOUR * dxi * deta);
    d2ij.y_der.x                        = (grid_pts[7].x - (TWO * grid_pts[4].x) +
                                           grid_pts[1].x)/(deta * deta);
    d2ij.y_der.y                        = (grid_pts[7].y - (TWO * grid_pts[4].y) +
                                           grid_pts[1].y)/(deta * deta);


    return d2ij;
}






/*
   Compute second and mixed derivatives of all (p, q)
   grid points with respect to the computational grid
   (xi, eta) using finite differences

   Input parameters: nx   - number of x points
                     ny   - number of y points
                     pq   - (x, y) grid coordinates
                     d2pq - array containing grid derivatives
*/
void pq_second_ders_2D (int nx, int ny, point_2D **pq, grid_dder_2D ***d2pq)
{
    /* Local variables */
    int                     i, j;
    long double             dxi, deta;


    /* Spacing of xi and eta in the uniform
       computational grid */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            (*d2pq)[i][j].x_der         = zero (point_2D);
            (*d2pq)[i][j].xy_der        = zero (point_2D);
            (*d2pq)[i][j].y_der         = zero (point_2D);
        }
    }


    /* Interior: O(h^2) central differences */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            (*d2pq)[i][j].x_der         = compute_finite_diff_2 (1, i, j, pq, 1);
            (*d2pq)[i][j].xy_der        = compute_finite_diff_2 (7, i, j, pq, 1);
            (*d2pq)[i][j].y_der         = compute_finite_diff_2 (4, i, j, pq, 1);
        }
    }


    /* eta = 0 and eta = 1 boundaries */
    for (i = 1; i < nx - 1; i++)
    {
        /* eta = 0 (j = 0): O(h^2) central difference in xi
           and O(h^2) forward difference in eta */
        (*d2pq)[i][0].x_der             = compute_finite_diff_2 (1, i, 0, pq, 1);
        (*d2pq)[i][0].xy_der            = compute_finite_diff_2 (8, i, 0, pq, 1);
        (*d2pq)[i][0].y_der             = compute_finite_diff_2 (5, i, 0, pq, 1);

        /* eta = 1 (j = ny - 1): O(h^2) central difference
           in xi and O(h^2) backward difference in eta */
        (*d2pq)[i][ny - 1].x_der        = compute_finite_diff_2 (1, i, ny - 1, pq, 1);
        (*d2pq)[i][ny - 1].xy_der       = compute_finite_diff_2 (9, i, ny - 1, pq, 1);
        (*d2pq)[i][ny - 1].y_der        = compute_finite_diff_2 (6, i, ny - 1, pq, 1);
    }


    /* xi = 0 and xi = 1 boundaries */
    for (j = 1; j < ny - 1; j++)
    {
        /* xi = 0 (i = 0): O(h^2) forward difference in xi
           and O(h^2) central difference in eta */
        (*d2pq)[0][j].x_der             = compute_finite_diff_2 (2, 0, j, pq, 1);
        (*d2pq)[0][j].xy_der            = compute_finite_diff_2 (10, 0, j, pq, 1);
        (*d2pq)[0][j].y_der             = compute_finite_diff_2 (4, 0, j, pq, 1);

        /* xi = 1 (i = nx - 1): O(h^2) backward difference
           in xi and O(h^2) central difference in eta */
        (*d2pq)[nx - 1][j].x_der        = compute_finite_diff_2 (3, nx - 1, j, pq, 1);
        (*d2pq)[nx - 1][j].xy_der       = compute_finite_diff_2 (11, nx - 1, j, pq, 1);
        (*d2pq)[nx - 1][j].y_der        = compute_finite_diff_2 (4, nx - 1, j, pq, 1);
    }


    /* xi = 0, eta = 0: forward difference in xi
       and forward difference in eta */
    (*d2pq)[0][0].x_der                 = compute_finite_diff_2 (2, 0, 0, pq, 1);
    (*d2pq)[0][0].xy_der                = compute_finite_diff_2 (12, 0, 0, pq, 1);
    (*d2pq)[0][0].y_der                 = compute_finite_diff_2 (5, 0, 0, pq, 1);

    /* xi = 1, eta = 0: backward difference in xi
       and forward difference in eta */
    (*d2pq)[nx - 1][0].x_der            = compute_finite_diff_2 (3, nx - 1, 0, pq, 1);
    (*d2pq)[nx - 1][0].xy_der           = compute_finite_diff_2 (13, nx - 1, 0, pq, 1);
    (*d2pq)[nx - 1][0].y_der            = compute_finite_diff_2 (5, nx - 1, 0, pq, 1);

    /* xi = 0, eta = 1: forward difference in xi
       and backward difference in eta */
    (*d2pq)[0][ny - 1].x_der            = compute_finite_diff_2 (2, 0, ny - 1, pq, 1);
    (*d2pq)[0][ny - 1].xy_der           = compute_finite_diff_2 (14, 0, ny - 1, pq, 1);
    (*d2pq)[0][ny - 1].y_der            = compute_finite_diff_2 (6, 0, ny - 1, pq, 1);

    /* xi = 1, eta = 1: backward difference in xi
       and backward difference in eta */
    (*d2pq)[nx - 1][ny - 1].x_der       = compute_finite_diff_2 (3, nx - 1, ny - 1, pq, 1);
    (*d2pq)[nx - 1][ny - 1].xy_der      = compute_finite_diff_2 (15, nx - 1, ny - 1, pq, 1);
    (*d2pq)[nx - 1][ny - 1].y_der       = compute_finite_diff_2 (6, nx - 1, ny - 1, pq, 1);


    /* Complete FD approximation */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            (*d2pq)[i][j].x_der         = divide_point_2D ((*d2pq)[i][j].x_der, (dxi * dxi));
            (*d2pq)[i][j].xy_der        = divide_point_2D ((*d2pq)[i][j].xy_der, (FOUR * dxi * deta));
            (*d2pq)[i][j].y_der         = divide_point_2D ((*d2pq)[i][j].y_der, (deta * deta));
        }
    }
}






/*
   Compute the derivatives of the computational
   grid (xi, eta) with respect to the outward
   unit normal to the boundaries of the physical
   domain which map to eta = 0 and eta = 1

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     grid       - (x, y) grid
                     nx_bounds  - unit normals to the given physical boundaries
*/
point_2D **compute_xbounds_dn (int debug, int nx, int ny, grid_der_2D **dgrid,
                               coeffs_1 **coeffs, point_2D **nx_bounds)
{
    /* Return xbounds_dn */
    point_2D                **xbounds_dn;

    /* Local variables */
    point_2D                **dcomp_dx, **dcomp_dy;
    int                     i;


    if (debug == 0) printf ("Debugging from compute_xbounds_dn!!!\n");


    /* Initialization */
    xbounds_dn                          = allocate_2D_point_2D_array ("xbounds_dn", 2, nx);
    dcomp_dx                            = allocate_2D_point_2D_array ("dcomp_dx", 2, nx);
    dcomp_dy                            = allocate_2D_point_2D_array ("dcomp_dy", 2, nx);


    /* Derivatives of (xi, eta) wrt (x, y) */
    for (i = 0; i < nx; i++)
    {
        /* Derivatives for the boundary which maps
           to eta = 0 */
        dcomp_dx[0][i].x                = dgrid[i][0].y_der.y/coeffs[i][0].J;   // dxi/dx
        dcomp_dx[0][i].y                = -dgrid[i][0].x_der.y/coeffs[i][0].J;  // deta/dx
        dcomp_dy[0][i].x                = -dgrid[i][0].y_der.x/coeffs[i][0].J;  // dxi/dy
        dcomp_dy[0][i].y                = dgrid[i][0].x_der.x/coeffs[i][0].J;   // deta/dy

        /* Derivatives wrt for the boundary which maps
           to eta = 1 */
        dcomp_dx[1][i].x                = dgrid[i][ny - 1].y_der.y/coeffs[i][ny - 1].J;   // dxi/dx
        dcomp_dx[1][i].y                = -dgrid[i][ny - 1].x_der.y/coeffs[i][ny - 1].J;  // deta/dx
        dcomp_dy[1][i].x                = -dgrid[i][ny - 1].y_der.x/coeffs[i][ny - 1].J;  // dxi/dy
        dcomp_dy[1][i].y                = dgrid[i][ny - 1].x_der.x/coeffs[i][ny - 1].J;   // deta/dy
    }


    /* Derivatives of (xi, eta) wrt the normals to the
       boundaries which map to eta = 0 and eta = 1 */
    for (i = 0; i < nx; i++)
    {
        /* eta = 0 */
        xbounds_dn[0][i].x              = (nx_bounds[0][i].x * dcomp_dx[0][i].x) +
                                          (nx_bounds[0][i].y * dcomp_dy[0][i].x);   // dxi/dn
        xbounds_dn[0][i].y              = (nx_bounds[0][i].x * dcomp_dx[0][i].y) +
                                          (nx_bounds[0][i].y * dcomp_dy[0][i].y);   // deta/dn

        /* eta = 1 */
        xbounds_dn[1][i].x              = (nx_bounds[1][i].x * dcomp_dx[1][i].x) +
                                          (nx_bounds[1][i].y * dcomp_dy[1][i].x);   // dxi/dn
        xbounds_dn[1][i].y              = (nx_bounds[1][i].x * dcomp_dx[1][i].y) +
                                          (nx_bounds[1][i].y * dcomp_dy[1][i].y);   // deta/dn
        if (debug == 0) printf ("%8.6LF    %8.6LF    %8.6LF    %8.6LF\n",
                                nx_bounds[1][i].x, dcomp_dx[1][i].x, nx_bounds[1][i].y, dcomp_dy[1][i].x);
    }


    /* Free memory */
    free_2D_point_2D_array ("dcomp_dx", 2, dcomp_dx);
    free_2D_point_2D_array ("dcomp_dy", 2, dcomp_dy);


    return xbounds_dn;
}






/*
   Compute the derivatives of the computational
   grid (xi, eta) with respect to the outward
   unit normal to the boundaries of the physical
   domain which map to xi = 0 and xi = 1

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     grid       - (x, y) grid
                     ny_bounds  - unit normals to the given physical boundaries
*/
point_2D **compute_ybounds_dn (int debug, int nx, int ny, grid_der_2D **dgrid,
                               coeffs_1 **coeffs, point_2D **ny_bounds)
{
    /* Return ybounds_dn */
    point_2D                **ybounds_dn;

    /* Local variables */
    point_2D                **dcomp_dx, **dcomp_dy;
    int                     j;


    if (debug == 0) printf ("Debugging from compute_ybounds_dn!!!\n");


    /* Initialization */
    ybounds_dn                          = allocate_2D_point_2D_array ("ybounds_dn", 2, ny);
    dcomp_dx                            = allocate_2D_point_2D_array ("dcomp_dx", 2, ny);
    dcomp_dy                            = allocate_2D_point_2D_array ("dcomp_dy", 2, ny);


    /* Derivatives of (xi, eta) wrt (x, y) */
    for (j = 0; j < ny; j++)
    {
        /* Derivatives for boundary which maps to
           xi = 0 */
        dcomp_dx[0][j].x                = dgrid[0][j].y_der.y/coeffs[0][j].J;   // dxi/dx
        dcomp_dx[0][j].y                = -dgrid[0][j].x_der.y/coeffs[0][j].J;  // deta/dx
        dcomp_dy[0][j].x                = -dgrid[0][j].y_der.x/coeffs[0][j].J;  // dxi/dy
        dcomp_dy[0][j].y                = dgrid[0][j].x_der.x/coeffs[0][j].J;   // deta/dy

        /* Derivatives for boundary which maps to
           xi = 1 */
        dcomp_dx[1][j].x                = dgrid[nx - 1][j].y_der.y/coeffs[nx - 1][j].J;   // dxi/dx
        dcomp_dx[1][j].y                = -dgrid[nx - 1][j].x_der.y/coeffs[nx - 1][j].J;  // deta/dx
        dcomp_dy[1][j].x                = -dgrid[nx - 1][j].y_der.x/coeffs[nx - 1][j].J;  // dxi/dy
        dcomp_dy[1][j].y                = dgrid[nx - 1][j].x_der.x/coeffs[nx - 1][j].J;   // deta/dy
    }


    /* Derivatives of (xi, eta) wrt the normals to the
       boundaries which map to xi = 0 and xi = 1 */
    for (j = 0; j < ny; j++)
    {
        /* xi = 0 */
        ybounds_dn[0][j].x              = (ny_bounds[0][j].x * dcomp_dx[0][j].x) +
                                          (ny_bounds[0][j].y * dcomp_dy[0][j].x);   // dxi/dn
        ybounds_dn[0][j].y              = (ny_bounds[0][j].x * dcomp_dx[0][j].y) +
                                          (ny_bounds[0][j].y * dcomp_dy[0][j].y);   // deta/dn

        /* xi = 1 */
        ybounds_dn[1][j].x              = (ny_bounds[1][j].x * dcomp_dx[1][j].x) +
                                          (ny_bounds[1][j].y * dcomp_dy[1][j].x);   // dxi/dn
        ybounds_dn[1][j].y              = (ny_bounds[1][j].x * dcomp_dx[1][j].y) +
                                          (ny_bounds[1][j].y * dcomp_dy[1][j].y);   // deta/dn
    }


    /* Free memory */
    free_2D_point_2D_array ("dcomp_dx", 2, dcomp_dx);
    free_2D_point_2D_array ("dcomp_dy", 2, dcomp_dy);


    return ybounds_dn;
}






/*
   Compute the values of (p, q) on the boundaries
   of the computational domain. These act as
   boundary conditions for the (p, q) elliptic
   system

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     debug      - switch to print debugging info from function
                     dgrid      - array containing values of first derivatives
                                  of grid (x, y) with respect to (xi, eta)
                     coeffs     - array containing values of first derivative
                                  coefficients
                     d2grid     - array containing values of second derivatives
                                  of grid (x, y) with respect to (xi, eta)
                     pq         - array containing interior (p, q) values
*/
void compute_pq_boundaries (int nx, int ny, int debug, grid_der_2D **dgrid, coeffs_1 **coeffs,
                            grid_dder_2D **d2grid, long double cxi, long double ceta, point_2D ***pq)
{
    /* Local variables */
    int                     i, j, k;
    point_2D                **nx_bounds, **ny_bounds, **xbounds_dn, **ybounds_dn;


    if (debug == 0) printf ("Debugging from compute_pq_boundaries!!!\n");

    /* Compute Laplacian of xi and eta at the
       boundaries which map to eta = 0 and eta = 1 */
    for (j = 0; j < 2; j++)
    {
        k                               = j * (ny - 1);
        for (i = 0; i < nx; i++)
        {
            (*pq)[i][k]                 = compute_laplacians_ij (dgrid[i][k], coeffs[i][k], d2grid[i][k]);
        }
    }

    /* Compute Laplacian of xi and eta at the
       boundaries which map to xi = 0 and xi = 1 */
    for (i = 0; i < 2; i++)
    {
        k                               = i * (nx - 1);
        for (j = 1; j < ny - 1; j++)
        {
            (*pq)[k][j]                 = compute_laplacians_ij (dgrid[k][j], coeffs[k][j], d2grid[k][j]);
        }
    }


    /* Read normals at all the boundaries */
    nx_bounds                           = read_xbound_normals_2D (nx);
    ny_bounds                           = read_ybound_normals_2D (ny);

    /* Compute derivatives of (xi, eta) wrt to the
       normals of the physical domain boundaries */
    xbounds_dn                          = compute_xbounds_dn (debug, nx, ny, dgrid, coeffs, nx_bounds);
    ybounds_dn                          = compute_ybounds_dn (debug, nx, ny, dgrid, coeffs, ny_bounds);


    /* Compute p and q at the boundaries */
    for (i = 0; i < nx; i++)
    {
        /* eta = 0 (j = 0) */
        (*pq)[i][0].x                   = (*pq)[i][0].x - (ceta * xbounds_dn[0][i].x);
        //(*pq)[i][0].y                   = (*pq)[i][0].y + (c * xbounds_dn[0][i].y);

        /* eta = 1 (j = ny - 1) */
        (*pq)[i][ny - 1].x              = (*pq)[i][ny - 1].x - (ceta * xbounds_dn[1][i].x);
        //(*pq)[i][ny - 1].y              = (*pq)[i][ny - 1].y + (c * xbounds_dn[1][i].y);
    }

    for (j = 0; j < ny; j++)

    {
        /* xi = 0 (i = 0) */
        //(*pq)[0][j].x                   = (*pq)[0][j].x + (c * ybounds_dn[0][j].x);
        (*pq)[0][j].y                   = (*pq)[0][j].y - (cxi * ybounds_dn[0][j].y);

        /* xi = 1 (i = nx - 1) */
        //(*pq)[nx - 1][j].x              = (*pq)[nx - 1][j].x + (c * ybounds_dn[1][j].x);
        (*pq)[nx - 1][j].y              = (*pq)[nx - 1][j].y - (cxi * ybounds_dn[1][j].y);
    }


    /* Free memory */
    free_2D_point_2D_array ("nx_bounds", 2, nx_bounds);
    free_2D_point_2D_array ("ny_bounds", 2, ny_bounds);
    free_2D_point_2D_array ("xbounds_dn", 2, xbounds_dn);
    free_2D_point_2D_array ("ybounds_dn", 2, ybounds_dn);
}






/*
   Compute the values of (p, q) at the boundaries
   of the grid generation domain. p at the xi
   constant boundaries and q on the eta constant
   boundaries are computed using TFI instead
   of Laplacians of xi and eta

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dgrid  - array containing derivatives of (x, y)
                              with respect to (xi, eta)
                     coeffs - array containing first derivative based
                              coefficients
                     d2grid - array containing second-order derivatives
                              of (x, y) with respect to (xi, eta)
                     cxi    - normal derivative multiplier on the
                              xi constant boundaries
                     ceta   - normal derivative multiplier on the
                              eta constant boundaries

   TODO: Cannot be used until Jacobian formulation
         is updated
*/
void compute_pq_boundaries_TFI (int nx, int ny, grid_der_2D **dgrid, coeffs_1 **coeffs,
                                grid_dder_2D **d2grid, long double cxi, long double ceta,
                                point_2D ***pq)
{
    /* Local variables */
    int                     i, j;
    point_2D                **nx_bounds, **ny_bounds, **xbounds_dn, **ybounds_dn,
                            temp, **xieta;


    /* Compute Laplacians of eta at the boundaries
       which map to eta = 0, 1 */
    for (j = 0; j < ny; j += ny - 1)
    {
        for (i = 0; i < nx; i++)
        {
            temp                        = compute_laplacians_ij (dgrid[i][j], coeffs[i][j], d2grid[i][j]);
            (*pq)[i][j].x               = temp.x;
        }
    }

    /* Compute Laplacians of xi at the boundaries
       which map to xi = 0, 1 */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 0; j < ny; j++)
        {
            temp                        = compute_laplacians_ij (dgrid[i][j], coeffs[i][j], d2grid[i][j]);
            (*pq)[i][j].y               = temp.y;
        }
    }

    /* Read normals at all boundaries */
    nx_bounds                           = read_xbound_normals_2D (nx);
    ny_bounds                           = read_ybound_normals_2D (ny);

    /* Compute derivatives of (xi, eta) wrt to the
       normals of the physical domain boundaries */
    xbounds_dn                          = compute_xbounds_dn (1, nx, ny, dgrid, coeffs, nx_bounds);
    ybounds_dn                          = compute_ybounds_dn (1, nx, ny, dgrid, coeffs, ny_bounds);


    /* Compute computational grid */
    xieta                               = computational_grid (nx, ny, 0);


    /* Compute p at the boundaries which map
       to eta = 0, 1 */
    for (i = 0; i < nx; i++)
    {
        /* eta = 0 (j = 0) */
        (*pq)[i][0].x                   = (*pq)[i][0].x - (ceta * xbounds_dn[0][i].x);

        /* eta = 1 (j = ny - 1) */
        (*pq)[i][ny - 1].x              = (*pq)[i][ny - 1].x - (ceta * xbounds_dn[1][i].x);
    }

    /* Compute p at the boundaries which map to
       xi = 0, 1 using TFI */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 1; j < ny - 1; j++)
        {
            (*pq)[i][j].x               = ((ONE - xieta[i][j].y) * (*pq)[i][0].x) +
                                          (xieta[i][j].y * (*pq)[i][ny - 1].x);
        }
    }


    /* Compute q at the boundaries which map
       to xi = 0, 1 */
    for (j = 0; j < ny; j++)
    {
        /* xi = 0 (i = 0) */
        (*pq)[0][j].y                   = (*pq)[0][j].y - (cxi * ybounds_dn[0][j].y);

        /* xi = nx - 1 (i = nx - 1) */
        (*pq)[nx - 1][j].y              = (*pq)[nx - 1][j].y - (cxi * ybounds_dn[1][j].y);
    }

    /* Compute q at the boundaries which map to
       eta = 0, 1 using TFI */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j += ny - 1)
        {
            (*pq)[i][j].y               = ((ONE - xieta[i][j].x) * (*pq)[0][j].y) +
                                          (xieta[i][j].x * (*pq)[nx - 1][j].y);
        }
    }


    /* Free memory */
    free_2D_point_2D_array ("nx_bounds", 2, nx_bounds);
    free_2D_point_2D_array ("ny_bounds", 2, ny_bounds);
    free_2D_point_2D_array ("xbounds_dn", 2, xbounds_dn);
    free_2D_point_2D_array ("ybounds_dn", 2, ybounds_dn);
}






/*
   Compute initial distribution of the control
   functions (p, q)

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     cxi    - normal derivative multiplier on the
                              xi constant boundaries
                     ceta   - nomrla derivative multiplier on the
                              eta constant boundaries
                     grid   - initial (x, y) grid
*/
void initial_pq_2D (int bc, int nx, int ny, long double cxi, long double ceta,
                    point_2D **grid, point_2D ***pq)
{
    /* Local variables */
    grid_der_2D             **dgrid;
    grid_dder_2D            **d2grid;
    coeffs_1                **coeffs;
    int                     i, j;
    long double             temp1, temp2;


    /* Initialization */
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);

    /* Compute grid first and second derivatives */
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);
    grid_second_ders_2D (nx, ny, grid, &d2grid);


    /* Compute initial (p, q) in the interior */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp1                       = (-ONE/(coeffs[i][j].J * coeffs[i][j].J)) *
                                          ((coeffs[i][j].alpha * d2grid[i][j].x_der.x) -
                                           (TWO * coeffs[i][j].beta * d2grid[i][j].xy_der.x) +
                                           (coeffs[i][j].gamma * d2grid[i][j].y_der.x));
            temp2                       = (-ONE/(coeffs[i][j].J * coeffs[i][j].J)) *
                                          ((coeffs[i][j].alpha * d2grid[i][j].x_der.y) -
                                           (TWO * coeffs[i][j].beta * d2grid[i][j].xy_der.y) +
                                           (coeffs[i][j].gamma * d2grid[i][j].y_der.y));

            (*pq)[i][j].x               = ((dgrid[i][j].y_der.y * temp1) - (dgrid[i][j].y_der.x * temp2))/(coeffs[i][j].J);
            (*pq)[i][j].y               = ((dgrid[i][j].x_der.x * temp2) - (dgrid[i][j].x_der.y * temp1))/(coeffs[i][j].J);
        }
    }

    /* Compute initial (p, q) at the boundaries */
    if (bc == 1) {
        compute_pq_boundaries (nx, ny, 1, dgrid, coeffs, d2grid, cxi, ceta, pq);
    } else {
        compute_pq_boundaries_TFI (nx, ny, dgrid, coeffs, d2grid, cxi, ceta, pq);
    }

    /* Free memory from derivative arrays */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
}






void initial_pq_2D_TFI (int nx, int ny, long double cxi, long double ceta,
                        point_2D **grid, point_2D ***pq)
{
    /* Local variables */
    grid_der_2D             **dgrid;
    grid_dder_2D            **d2grid;
    coeffs_1                **coeffs;
    int                     i, j;
    long double             temp1, temp2;
    point_2D                **xieta;


    /* Initialization */
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);

    /* Compute initial (p, q) in the interior */
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);
    grid_second_ders_2D (nx, ny, grid, &d2grid);


    /* Compute initial (p, q) at the boundaries */
    compute_pq_boundaries (nx, ny, 1, dgrid, coeffs, d2grid, cxi, ceta, pq);

    /* Computational grid */
    xieta                               = computational_grid (nx, ny, 0);

    /* Interpolate values of p in the interior */
    for (i = 0; i < nx; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp2                       = ONE - xieta[i][j].y;
            (*pq)[i][j].x               = (temp2 * (*pq)[i][0].x) + (xieta[i][j].y * (*pq)[i][ny - 1].x);
        }
    }

    /* Interpolate values of q in the interior */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j++)
        {
            temp1                       = ONE - xieta[i][j].x;
            (*pq)[i][j].y               = (temp1 * (*pq)[0][j].y) + (xieta[i][j].x * (*pq)[nx - 1][j].y);
        }
    }

    /* Free memory from derivative arrays */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
}






/*
   Compute the maximum residual for a
   solver iteration

   Input parameters: nx         - number of points in x
                     ny         - number of points in y
                     prev_grid  - array containing values of (x, y)
                                  at iteration k
                     grid       - array containing values of (x, y)
                                  at iteration k + 1
*/
long double compute_max_residual (int nx, int ny, point_2D **prev_grid,
                                  point_2D **grid)
{
    /* Return */
    long double             max_res;

    /* Local variables */
    int                     i, j;
    long double             prev_xy[nx][ny], current_xy[nx][ny],
                            **res_xy;


    /* Allocate memory for res_xy array */
    res_xy                          = allocate_2D_long_double_array ("res_xy", nx, ny);

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            /* Compute (x, y) residual arrays */
            prev_xy[i][j]           = (long double) sqrt(pow((double) prev_grid[i][j].x, (double) 2.0) +
                                                         pow((double) prev_grid[i][j].y, (double) 2.0));
            current_xy[i][j]        = (long double) sqrt(pow((double) grid[i][j].x, (double) 2.0) +
                                                         pow((double) grid[i][j].y, (double) 2.0));
            res_xy[i][j]            = (long double) fabs ((double) current_xy[i][j] - (double) prev_xy[i][j]);
        }
    }

    max_res                         = array_max_2D_long_double (nx, ny, res_xy);

    /* Free memory from dynamically allocated array */
    free_2D_long_double_array ("res_xy", nx, res_xy);


    return max_res;
}






/*
   Compute initial distribution of the control
   functions (p, q) using a Laplace equation
   for the domain interior

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     cxi    - normal derivative multiplier on the
                              xi constant boundaries
                     ceta   - normal derivative multiplier on the
                              eta constant boundaries
                     grid   - initial (x, y) grid
*/
void initial_pq_2D_laplace (int nx, int ny, long double cxi, long double ceta,
                            point_2D **grid, point_2D ***pq)
{
    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    point_2D                **prev_pq;
    int                     i, j, k;
    const int               niter = 250;
    long double             max_res, g1, g2, g3, g;


    /* Initialization */
    (*pq)                               = allocate_2D_point_2D_array ("pq", nx, ny);
    initial_pq_2D (1, nx, ny, cxi, ceta, grid, pq);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    prev_pq                             = allocate_2D_point_2D_array ("prev_pq", nx, ny);

    /* Compute grid derivatives */
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);

    /* Compute numerical grid */
    long double dxi                     = ONE/((long double) (nx - 1));
    long double deta                    = ONE/((long double) (ny - 1));


    /* Gauss-Seidel iteration for solving for
       initial distribution of (p, q) */
    for (k = 0; k < niter; k++)
    {
        /* Set previous iteration (p, q) */
        equals_2D_point_2D_array (nx, ny, *pq, &prev_pq);

        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                g                       = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
                g1                      = coeffs[i][j].alpha/g;
                g2                      = TWO * (coeffs[i][j].beta/g);
                g3                      = coeffs[i][j].gamma/g;

                (*pq)[i][j].x           = (g1 * (((*pq)[i + 1][j].x + (*pq)[i - 1][j].x)/(dxi * dxi))) -
                                          (g2 * (((*pq)[i + 1][j + 1].x - (*pq)[i + 1][j - 1].x -
                                           (*pq)[i - 1][j + 1].x + (*pq)[i - 1][j - 1].x)/(FOUR * dxi * deta))) +
                                          (g3 * (((*pq)[i][j + 1].x + (*pq)[i][j - 1].x)/(deta * deta)));
                (*pq)[i][j].y           = (g1 * (((*pq)[i + 1][j].y + (*pq)[i - 1][j].y)/(dxi * dxi))) -
                                          (g2 * (((*pq)[i + 1][j + 1].y - (*pq)[i + 1][j - 1].y -
                                           (*pq)[i - 1][j + 1].y + (*pq)[i - 1][j - 1].y)/(FOUR * dxi * deta))) +
                                          (g3 * (((*pq)[i][j + 1].y + (*pq)[i][j - 1].y)/(deta * deta)));
            }
        }

        /* Compute maximum residual for exit condition */
        max_res                     = compute_max_residual (nx, ny, prev_pq, *pq);
        if (max_res < 1E-8) break;

    }


    /* Free memory from derivative arrays */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_point_2D_array ("prev_pq", nx, prev_pq);
}






/*
   Write out the values of the derivative
   deta/dn along the boundaries mapping to
   xi = 0 and xi = 1 for the final grid

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dgrid  - 2D array of type point_2D containing
                              the values of the derivatives of
                              (x, y) wrt (xi, eta)
*/
void print_normals (int nx, int ny, grid_der_2D **dgrid)
{
    point_2D                **nx_bounds, **ny_bounds;
    long double             temp;
    FILE                    *fptr;
    int                     i, j;

    /* Read normals along the boundaries */
    nx_bounds                           = read_xbound_normals_2D (nx);
    ny_bounds                           = read_ybound_normals_2D (ny);

    /* Open file for writing final deta/dn values */
    fptr                                = fopen ("final_detadn.dat", "w");

    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 1; j < ny - 1; j++)
        {
            int i1                      = i/(nx - 1);
            temp                        = -(dgrid[i][j].x_der.y * ny_bounds[i1][j].x) +
                                          (dgrid[i][j].x_der.x * ny_bounds[i1][j].y);
            fprintf (fptr, "%d    %d    %10.6LF    %10.6LF    %10.6LF\n",
                     i, j, ny_bounds[i1][j].x, ny_bounds[i1][j].y, temp);
        }
    }

    /* Close file */
    fclose (fptr);

    /* Open file for writing dxi/dn values */
    fptr                                = fopen ("final_dxidn.dat", "w");

    for (j = 0; j < ny; j += ny - 1)
    {
        for (i = 0; i < nx; i++)
        {
            int j1                      = j/(ny - 1);
            temp                        = (dgrid[i][j].y_der.y * nx_bounds[j1][i].x) -
                                          (dgrid[i][j].y_der.x * nx_bounds[j1][i].y);
            fprintf (fptr, "%d    %d    %10.6LF    %10.6LF    %10.6LF\n",
                     i, j, nx_bounds[j1][i].x, nx_bounds[j1][i].y, temp);
        }
    }

    /* Close file */
    fclose (fptr);

    /* Free memory */
    free_2D_point_2D_array ("nx_bounds", 2, nx_bounds);
    free_2D_point_2D_array ("ny_bounds", 2, ny_bounds);
}






/*
   Compute Thomas-Middlecoff control functions
   for the grid generation domain

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dgrid  - derivatives of (x, y) with respect to
                              (xi, eta)
                     d2grid - second-order derivatives of (x, y) with
                              respect to (xi, eta)
*/
void compute_thomas_middlecoff (int nx, int ny, grid_der_2D **dgrid, grid_dder_2D **d2grid, point_2D ***pq)
{
    /* Local variables */
    int                     i, j;
    //const int               niter = 500;
    point_2D                **xieta, **prev_pq;
    long double             temp1, temp2;//, dxi, deta, g, g1, g2, g3, max_res;


    /* Allocate array memory */
    prev_pq                             = allocate_2D_point_2D_array ("prev_pq", nx, ny);


    /* Computational grid */
    xieta                               = computational_grid (nx, ny, 0);


    /* Compute p on the eta constant boundaries */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j += ny - 1)
        {
            temp1                       = (dgrid[i][j].x_der.x * d2grid[i][j].x_der.x) +
                                          (dgrid[i][j].x_der.y * d2grid[i][j].x_der.y);
            temp2                       = (dgrid[i][j].x_der.x * dgrid[i][j].x_der.x) +
                                          (dgrid[i][j].x_der.y * dgrid[i][j].x_der.y);
            (*pq)[i][j].x               = -temp1/(temp2 * temp2);
        }
    }

    /* Compute p on the xi constant boundaries */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp2                       = ONE - xieta[i][j].y;
            (*pq)[i][j].x               = (temp2 * (*pq)[i][0].x) + (xieta[i][j].y * (*pq)[i][ny - 1].x);
        }
    }

    /* Interpolate values of p in the interior */
    for (i = 0; i < nx; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp1                       = ONE - xieta[i][j].x;
            temp2                       = ONE - xieta[i][j].y;
            (*pq)[i][j].x               = (temp1 * (*pq)[0][j].x) + (xieta[i][j].x * (*pq)[nx - 1][j].x) +
                                          (temp2 * (*pq)[i][0].x) + (xieta[i][j].y * (*pq)[i][ny - 1].x) -
                                          (temp1 * temp2 * (*pq)[0][0].x) - (xieta[i][j].x * temp2 * (*pq)[nx - 1][0].x) -
                                          (temp1 * xieta[i][j].y * (*pq)[0][ny - 1].x) - (xieta[i][j].x * xieta[i][j].y *
                                           (*pq)[nx - 1][ny - 1].x);
        }
    }


    /* Compute q on the xi constant boundaries */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 0; j < ny; j++)
        {
            temp1                       = (dgrid[i][j].y_der.x * d2grid[i][j].y_der.x) +
                                          (dgrid[i][j].y_der.y * d2grid[i][j].y_der.y);
            temp2                       = (dgrid[i][j].y_der.x * dgrid[i][j].y_der.x) +
                                          (dgrid[i][j].y_der.y * dgrid[i][j].y_der.y);
            (*pq)[i][j].y               = -temp1/(temp2 * temp2);
        }
    }

    /* Compute q on the eta constant boundaries */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j += ny - 1)
        {
            temp1                       = ONE - xieta[i][j].x;
            (*pq)[i][j].y               = (temp1 * (*pq)[0][j].y) + (xieta[i][j].x * (*pq)[nx - 1][j].y);
        }
    }

    /* Interpolate values of q in the interior */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j++)
        {
            temp1                       = ONE - xieta[i][j].x;
            temp2                       = ONE - xieta[i][j].y;
            (*pq)[i][j].y               = (temp1 * (*pq)[0][j].y) + (xieta[i][j].x * (*pq)[nx - 1][j].y) +
                                          (temp2 * (*pq)[i][0].y) + (xieta[i][j].y * (*pq)[i][ny - 1].y) -
                                          (temp1 * temp2 * (*pq)[0][0].y) - (xieta[i][j].x * temp2 * (*pq)[nx - 1][0].y) -
                                          (temp1 * xieta[i][j].y * (*pq)[0][ny - 1].y) - (xieta[i][j].x * xieta[i][j].y *
                                           (*pq)[nx - 1][ny - 1].y);
        }
    }


    /* Compute values of (p, q) in the domain
       interior with the Laplace equation solved
       using Gauss-Seidel method */
    //for (k = 0; k < niter; k++)
    //{
    //    /* Set previous iteration (p, q) */
    //    equals_2D_point_2D_array (nx, ny, *pq, &prev_pq);

    //    for (i = 1; i < nx - 1; i++)
    //    {
    //        for (j = 1; j < ny - 1; j++)
    //        {
    //            g                       = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
    //                                             (coeffs[i][j].gamma/(deta * deta)));
    //            g1                      = coeffs[i][j].alpha/g;
    //            g2                      = TWO * (coeffs[i][j].beta/g);
    //            g3                      = coeffs[i][j].gamma/g;

    //            (*pq)[i][j].x           = (g1 * (((*pq)[i + 1][j].x + (*pq)[i - 1][j].x)/(dxi * dxi))) -
    //                                      (g2 * (((*pq)[i + 1][j + 1].x - (*pq)[i + 1][j - 1].x -
    //                                       (*pq)[i - 1][j + 1].x + (*pq)[i - 1][j - 1].x)/(FOUR * dxi * deta))) +
    //                                      (g3 * (((*pq)[i][j + 1].x + (*pq)[i][j - 1].x)/(deta * deta)));
    //            (*pq)[i][j].y           = (g1 * (((*pq)[i + 1][j].y + (*pq)[i - 1][j].y)/(dxi * dxi))) -
    //                                      (g2 * (((*pq)[i + 1][j + 1].y - (*pq)[i + 1][j - 1].y -
    //                                       (*pq)[i - 1][j + 1].y + (*pq)[i - 1][j - 1].y)/(FOUR * dxi * deta))) +
    //                                      (g3 * (((*pq)[i][j + 1].y + (*pq)[i][j - 1].y)/(deta * deta)));
    //        }
    //    }

    //    /* Compute maximum residual for exit condition */
    //    max_res                     = compute_max_residual (nx, ny, prev_pq, *pq);
    //    if (max_res < 1E-8) break;

    //}


    /* Free memory */
    free_2D_point_2D_array ("xieta", nx, xieta);
    free_2D_point_2D_array ("prev_pq", nx, prev_pq);
}






/*
   Compute Thomas-Middlecoff control functions
   for the grid generation domain in (\phi, \psi)
   form

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dgrid  - derivatives of (x, y) with respect to
                              (xi, eta)
                     d2grid - second-order derivatives of (x, y) with
                              respect to (xi, eta)
*/
void compute_thomas_middlecoff_alt (int nx, int ny, grid_der_2D **dgrid, grid_dder_2D **d2grid,
                                    point_2D ***phipsi)
{
    /* Local variables */
    int                     i, j;
    point_2D                **xieta;
    long double             temp1, temp2;


    /* Computational grid */
    xieta                               = computational_grid (nx, ny, 0);


    /* Compute phi on the eta constant boudnaries */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j += ny - 1)
        {
            temp1                       = (dgrid[i][j].x_der.x * d2grid[i][j].x_der.x) +
                                          (dgrid[i][j].x_der.y * d2grid[i][j].x_der.y);
            temp2                       = sqr(dgrid[i][j].x_der.x) + sqr(dgrid[i][j].x_der.y);
            (*phipsi)[i][j].x           = -temp1/temp2;
        }
    }

    /* Compute phi on the xi constant boundaries */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp2                       = ONE - xieta[i][j].y;
            (*phipsi)[i][j].x           = (temp2 * (*phipsi)[i][0].x) + (xieta[i][j].y * (*phipsi)[i][ny - 1].x);
        }
    }

    /* Interpolate values of phi in the interior */
    for (i = 0; i < nx; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp1                       = ONE - xieta[i][j].x;
            temp2                       = ONE - xieta[i][j].y;
            (*phipsi)[i][j].x           = (temp1 * (*phipsi)[0][j].x) + (xieta[i][j].x * (*phipsi)[nx - 1][j].x) +
                                          (temp2 * (*phipsi)[i][0].x) + (xieta[i][j].y * (*phipsi)[i][ny - 1].x) -
                                          (temp1 * temp2 * (*phipsi)[0][0].x) - (xieta[i][j].x * temp2 * (*phipsi)[nx - 1][0].x) -
                                          (temp1 * xieta[i][j].y * (*phipsi)[0][ny - 1].x) - (xieta[i][j].x * xieta[i][j].y *
                                           (*phipsi)[nx - 1][ny - 1].x);
        }
    }


    /* Compute psi on the xi constant boundaries */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 0; j < ny; j++)
        {
            temp1                       = (dgrid[i][j].y_der.x * d2grid[i][j].y_der.x) +
                                          (dgrid[i][j].y_der.y * d2grid[i][j].y_der.y);
            temp2                       = (sqr(dgrid[i][j].y_der.x) + sqr(dgrid[i][j].y_der.y));
            (*phipsi)[i][j].y           = -temp1/temp2;
        }
    }

    /* Compute psi on the eta constant boundaries */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j += ny - 1)
        {
            temp1                       = ONE - xieta[i][j].x;
            (*phipsi)[i][j].y           = (temp1 * (*phipsi)[0][j].y) + (xieta[i][j].x * (*phipsi)[nx - 1][j].y);
        }
    }

    /* Interpolate values of psi in the interior */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j++)
        {
            temp1                       = ONE - xieta[i][j].x;
            temp2                       = ONE - xieta[i][j].y;
            (*phipsi)[i][j].y           = (temp1 * (*phipsi)[0][j].y) + (xieta[i][j].x * (*phipsi)[nx - 1][j].y) +
                                          (temp2 * (*phipsi)[i][0].y) + (xieta[i][j].y * (*phipsi)[i][ny - 1].y) -
                                          (temp1 * temp2 * (*phipsi)[0][0].y) - (xieta[i][j].x * temp2 * (*phipsi)[nx - 1][0].y) -
                                          (temp1 * xieta[i][j].y * (*phipsi)[0][ny - 1].y) - (xieta[i][j].x * xieta[i][j].y *
                                           (*phipsi)[nx - 1][ny - 1].y);
        }
    }


    /* Free memory */
    free_2D_point_2D_array ("xieta", nx, xieta);
}






/*
   Compute Steger-Sorenson control functions
   for the grid generation domain

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     dgrid  - derivatives of (x, y) with respect to
                              (xi, eta)
                     d2grid - second-order derivatives of (x, y) with
                              respect to (xi, eta)
*/
void compute_steger_sorenson (int nx, int ny, grid_der_2D **dgrid, grid_dder_2D **d2grid, point_2D ***pq)
{
    /* Local variables */
    int                     i, j;
    //const int               niter = 500;
    point_2D                **xieta, **prev_pq;
    long double             temp1, temp2, temp3, temp4;//, dxi, deta, g, g1, g2, g3, max_res;


    /* Allocate array memory */
    prev_pq                             = allocate_2D_point_2D_array ("prev_pq", nx, ny);


    /* Computational grid */
    xieta                               = computational_grid (nx, ny, 0);


    /* Compute p on the eta constant boundaries */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j += ny - 1)
        {
            temp1                       = (dgrid[i][j].x_der.x * d2grid[i][j].x_der.x) +
                                          (dgrid[i][j].x_der.y * d2grid[i][j].x_der.y);
            temp2                       = sqr(dgrid[i][j].x_der.x) + sqr(dgrid[i][j].x_der.y);
            temp3                       = (dgrid[i][j].x_der.x * d2grid[i][j].y_der.x) +
                                          (dgrid[i][j].x_der.y * d2grid[i][j].y_der.y);
            temp4                       = sqr(dgrid[i][j].y_der.x) + sqr(dgrid[i][j].y_der.y);
            (*pq)[i][j].x               = -(temp1/temp2) - (temp3/temp4);
        }
    }

    /* Compute p on the xi constant boundaries */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp1                       = (dgrid[i][j].x_der.x * d2grid[i][j].x_der.x) +
                                          (dgrid[i][j].x_der.y * d2grid[i][j].x_der.y);
            temp2                       = sqr(dgrid[i][j].x_der.x) + sqr(dgrid[i][j].x_der.y);
            temp3                       = (dgrid[i][j].x_der.x * d2grid[i][j].y_der.x) +
                                          (dgrid[i][j].x_der.y * d2grid[i][j].y_der.y);
            temp4                       = sqr(dgrid[i][j].y_der.x) + sqr(dgrid[i][j].y_der.y);
            (*pq)[i][j].x               = -(temp1/temp2) - (temp3/temp4);
        }
    }

    /* Interpolate values of p in the interior */
    for (i = 0; i < nx; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            temp1                       = ONE - xieta[i][j].x;
            temp2                       = ONE - xieta[i][j].y;
            (*pq)[i][j].x               = (temp1 * (*pq)[0][j].x) + (xieta[i][j].x * (*pq)[nx - 1][j].x) +
                                          (temp2 * (*pq)[i][0].x) + (xieta[i][j].y * (*pq)[i][ny - 1].x) -
                                          (temp1 * temp2 * (*pq)[0][0].x) - (xieta[i][j].x * temp2 * (*pq)[nx - 1][0].x) -
                                          (temp1 * xieta[i][j].y * (*pq)[0][ny - 1].x) - (xieta[i][j].x * xieta[i][j].y *
                                           (*pq)[nx - 1][ny - 1].x);
        }
    }


    /* Compute q on the xi constant boundaries */
    for (i = 0; i < nx; i += nx - 1)
    {
        for (j = 0; j < ny; j++)
        {
            temp1                       = (dgrid[i][j].y_der.x * d2grid[i][j].y_der.x) +
                                          (dgrid[i][j].y_der.y * d2grid[i][j].y_der.y);
            temp2                       = sqr(dgrid[i][j].y_der.x) + sqr(dgrid[i][j].y_der.y);
            temp3                       = (dgrid[i][j].y_der.x * d2grid[i][j].x_der.x) +
                                          (dgrid[i][j].y_der.y * d2grid[i][j].x_der.y);
            temp4                       = sqr(dgrid[i][j].x_der.x) + sqr(dgrid[i][j].x_der.y);
            (*pq)[i][j].x               = -(temp1/temp2) - (temp3/temp4);
        }
    }

    /* Compute q on the eta constant boundaries */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j += ny - 1)
        {
            temp1                       = (dgrid[i][j].y_der.x * d2grid[i][j].y_der.x) +
                                          (dgrid[i][j].y_der.y * d2grid[i][j].y_der.y);
            temp2                       = sqr(dgrid[i][j].y_der.x) + sqr(dgrid[i][j].y_der.y);
            temp3                       = (dgrid[i][j].y_der.x * d2grid[i][j].x_der.x) +
                                          (dgrid[i][j].y_der.y * d2grid[i][j].x_der.y);
            temp4                       = sqr(dgrid[i][j].x_der.x) + sqr(dgrid[i][j].x_der.y);
            (*pq)[i][j].x               = -(temp1/temp2) - (temp3/temp4);
        }
    }

    /* Interpolate values of q in the interior */
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 0; j < ny; j++)
        {
            temp1                       = ONE - xieta[i][j].x;
            temp2                       = ONE - xieta[i][j].y;
            (*pq)[i][j].y               = (temp1 * (*pq)[0][j].y) + (xieta[i][j].x * (*pq)[nx - 1][j].y) +
                                          (temp2 * (*pq)[i][0].y) + (xieta[i][j].y * (*pq)[i][ny - 1].y) -
                                          (temp1 * temp2 * (*pq)[0][0].y) - (xieta[i][j].x * temp2 * (*pq)[nx - 1][0].y) -
                                          (temp1 * xieta[i][j].y * (*pq)[0][ny - 1].y) - (xieta[i][j].x * xieta[i][j].y *
                                           (*pq)[nx - 1][ny - 1].y);
        }
    }


    /* Free memory */
    free_2D_point_2D_array ("xieta", nx, xieta);
    free_2D_point_2D_array ("prev_pq", nx, prev_pq);
}






/*
   Compute the Euclidean norm of the residual
   for a Newton iteration. The residual is
   computed as the distance between corresponding
   points between the previous grid and the
   current grid

   Input parameters: nx     - number of x points
                     ny     - number of y points
                     del_u  - updates of x and y obtained from
                              Newton's method
*/
long double compute_residual_norm (int nx, int ny, point_2D *del_u)
{
    /* Return res_norm */
    long double             res_norm;

    /* Local variables */
    long double             *dist;


    /* Initialize residual norm as 0 */
    res_norm                            = ZERO;

    /* Allocate memory */
    dist                                = allocate_1D_long_double_array ("dist", nx * ny);

    for (int i = 0; i < nx * ny; i++)
    {
        dist[i]                         = (long double) sqrt((double) ((del_u[i].x * del_u[i].x) +
                                                                       (del_u[i].y * del_u[i].y)));
        res_norm                        += dist[i] * dist[i];
    }
    res_norm                           = (long double) sqrt((double) res_norm);

    /* Free memory */
    free_1D_long_double_array ("dist", dist);


    return res_norm;
}






/*
   Elliptic grid generation without any control
   functions using Newton's method

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial grid
                     niter      - number of Newton iterations
*/
point_2D **elliptic_grid_2D (int nx, int ny, point_2D **init_grid, int niter)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, **lhs_array, *rhs_array, *update, res_norm;
    point_2D                *fu, *del_u;
    ell_jacobian_2D         **dFdu;
    char                    file1[256];
    FILE                    *fptr;
    int                     iiter;


    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("elliptic_grid", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);

    /* Allocate memory for Newton's method arrays */
    dFdu                                = allocate_2D_ell_jacobian_2D_array ("dFdu", nx * ny, nx * ny);
    fu                                  = allocate_1D_point_2D_array ("fu", nx * ny);
    lhs_array                           = allocate_2D_long_double_array ("lhs_array", 2 * nx * ny, 2 * nx * ny);
    rhs_array                           = allocate_1D_long_double_array ("rhs_array", 2 * nx * ny);
    update                              = allocate_1D_long_double_array ("update", 2 * nx * ny);
    del_u                               = allocate_1D_point_2D_array ("del_u", nx * ny);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);

    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Solve for elliptic grids with no control functions */
    fptr                                = fopen ("residuals.dat", "w");
    write_2D_singleblock_vts ("initial_grid.vts", nx, ny, grid, 1); // Write initial grid to Paraview file
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid); // Write initial grid to .dat files
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Compute initial grid derivatives and
        first derivative coefficients */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        /* Newton's method */
        ell_construct_newton_matrix (nx, ny, dxi, deta, dgrid, coeffs, d2grid, &dFdu);
        ell_construct_newton_rhs_vector (nx, ny, grid, coeffs, d2grid, &fu, 0);
        gmres_linear_system_solve (0, nx, ny, dFdu, fu, &del_u);
        res_norm                        = compute_residual_norm (nx, ny, del_u);
        fprintf (fptr, "%d  %22.16LF\n", iiter + 1, res_norm);

        /* Update 2D (x, y) arrays */
        ell_update_solution (nx, ny, del_u, &grid);

        /* Compute new derivatives */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);
    }
    snprintf (file1, 256, "final_grid_%d.vts", iiter + 1);
    write_2D_singleblock_vts (file1, nx, ny, grid, 1); // Write final grid to Paraview file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files
    fclose (fptr);


    /* Free memory from grid arrays */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);

    /* Free memory from Newton's method arrays */
    free_2D_ell_jacobian_2D_array ("dFdu", nx * ny, dFdu);
    free_1D_point_2D_array ("fu", fu);
    free_2D_long_double_array ("lhs_array", 2 * nx * ny, lhs_array);
    free_1D_long_double_array ("rhs_array", rhs_array);
    free_1D_long_double_array ("update", update);
    free_1D_point_2D_array ("del_u", del_u);


    return grid;
}






/*
   Elliptic grid generation with zero control
   functions using pointwise iteration

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial (x, y) grid
                     niter      - number of Newton iterations
*/
point_2D **elliptic_grid_2D_point (int nx, int ny, point_2D **init_grid, int niter)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, g, g1, g2, g3;
    const long double       omega = 1.0;
    point_2D                *del_u, **prev_grid;
    char                    file1[256];
    //FILE                    *fptr;
    int                     iiter, i, j, k;


    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("elliptic_grid", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    del_u                               = allocate_1D_point_2D_array ("del_u", nx * ny);
    prev_grid                           = allocate_2D_point_2D_array ("prev_grid", nx, ny);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);


    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Solve for elliptic grid with control functions */
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid);   // Write initial grid to .dat files
    write_2D_singleblock_vts ("initial_grid.vts", nx, ny, grid, 1);
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Set previous iteration (x, y) */
        equals_2D_point_2D_array (nx, ny, grid, &prev_grid);

        /* Grid derivatives and coefficients */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        /* Compute new interior grid values */
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                g                       = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
                g1                      = coeffs[i][j].alpha/g;
                g2                      = (TWO * coeffs[i][j].beta)/g;
                g3                      = coeffs[i][j].gamma/g;

                grid[i][j].x            = (g1 * ((grid[i + 1][j].x + grid[i - 1][j].x)/(dxi * dxi))) -
                                          (g2 * d2grid[i][j].xy_der.x) +
                                          (g3 * ((grid[i][j + 1].x + grid[i][j - 1].x)/(deta * deta)));
                grid[i][j].y            = (g1 * ((grid[i + 1][j].y + grid[i - 1][j].y)/(dxi * dxi))) -
                                          (g2 * d2grid[i][j].xy_der.y) +
                                          (g3 * ((grid[i][j + 1].y + grid[i][j - 1].y)/(deta * deta)));
            }
        }

        /* Under-relaxation */
        for (k = 0; k < nx * ny; k++)
        {
            i                           = k%nx;
            j                           = k/nx;
            grid[i][j].x                = (omega * grid[i][j].x) + ((ONE - omega) * prev_grid[i][j].x);
            grid[i][j].y                = (omega * grid[i][j].y) + ((ONE - omega) * prev_grid[i][j].y);

            /* Final iteration update vector */
            if (iiter == niter - 1)
            {
                del_u[k].x                  = grid[i][j].x - prev_grid[i][j].x;
                del_u[k].y                  = grid[i][j].y - prev_grid[i][j].y;
            }
        }
    }
    res_norm                                = compute_residual_norm (nx, ny, del_u); // Final iteration residual
    printf ("res_norm = %22.16LF\n", res_norm);
    snprintf (file1, 256, "final_grid_%d.vts", iiter);
    write_2D_singleblock_vts (file1, nx, ny, grid, 1);
    write_2D_singleblock_plot3D ("final_grid.x", nx, ny, grid);
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files


    /* Free memory */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_1D_point_2D_array ("del_u", del_u);
    free_2D_point_2D_array ("prev_grid", nx, prev_grid);


    return grid;
}






/*
   Elliptic grid generation with non-zero control
   functions using Newton's method

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial (x, y) grid
                     niter      - number of Newton iterations
*/
point_2D **poisson_grid_2D (int nx, int ny, point_2D **init_grid, int niter)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, **lhs_array, *rhs_array, *update, res_norm;
    point_2D                *fu, *del_u, **pq;
    ell_jacobian_2D         **dFdu;
    char                    file1[256];
    FILE                    *fptr;
    int                     iiter;


    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("poisson_grid", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    pq                                  = allocate_2D_point_2D_array ("pq", nx, ny);

    /* Allocate memory for Newton's method arrays */
    dFdu                                = allocate_2D_ell_jacobian_2D_array ("dFdu", nx * ny, nx * ny);
    fu                                  = allocate_1D_point_2D_array ("fu", nx * ny);
    lhs_array                           = allocate_2D_long_double_array ("lhs_array", 2 * nx * ny, 2 * nx * ny);
    rhs_array                           = allocate_1D_long_double_array ("rhs_array", 2 * nx * ny);
    update                              = allocate_1D_long_double_array ("lhs_array", 2 * nx * ny);
    del_u                               = allocate_1D_point_2D_array ("del_u", nx * ny);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);


    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Solve for elliptic grid with control functions */
    fptr                                = fopen ("residuals.dat", "w"); // File for writing residuals
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid); // Write initial grid to .dat files
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Grid derivatives and coefficients */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        /* Compute control functions (p, q) */
        compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &pq);
        if (iiter == 0)
            write_point_2D_to_vts ("initial_pq.vts", nx, ny, grid, "pq", pq); // Write initial (p, q) to a Paraview file

        /* Newton's method */
        psn_construct_newton_matrix (nx, ny, dxi, deta, dgrid, coeffs, d2grid, pq, &dFdu);
        psn_construct_newton_rhs_vector (nx, ny, grid, dgrid, coeffs, d2grid, pq, &fu);
        convert_ell_jacobian_2D (nx * ny, dFdu, &lhs_array);
        convert_point_2D (nx * ny, fu, &rhs_array);
        LU_linear_system_solve (2 * nx * ny, lhs_array, rhs_array, &update);
        convert_to_point_2D (nx * ny, update, &del_u);
        //gmres_linear_system_solve (0, nx, ny, dFdu, fu, &del_u);
        res_norm                        = compute_residual_norm (nx, ny, del_u);
        fprintf (fptr, "%d    %22.16LF\n", iiter + 1, res_norm);

        /* Update 2D (x, y) arrays */
        ell_update_solution (nx, ny, del_u, &grid);
    }
    snprintf (file1, 256, "final_pq_%d.vts", iiter);
    write_point_2D_to_vts (file1, nx, ny, grid, "pq", pq);  // Write final (p, q) to a Paraview file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files
    fclose (fptr);  // Close residual output file


    /* Free memory from grid arrays */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_point_2D_array ("pq", nx, pq);

    /* Free memory from Newton's method arrays */
    free_2D_ell_jacobian_2D_array ("dFdu", nx * ny, dFdu);
    free_1D_point_2D_array ("fu", fu);
    free_2D_long_double_array ("lhs_array", 2 * nx * ny, lhs_array);
    free_1D_long_double_array ("rhs_array", rhs_array);
    free_1D_long_double_array ("update", update);
    free_1D_point_2D_array ("del_u", del_u);


    return grid;
}






/*
   Elliptic grid generation with non-zero control
   functions using pointwise iteration

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial (x, y) grid
                     niter      - number of Newton iterations
*/
point_2D **poisson_grid_2D_point (int nx, int ny, point_2D **init_grid, int niter)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, g, g1, g2, g3, g4;
    const long double       omega = 1;
    point_2D                *del_u, **prev_grid, **pq;
    char                    file1[256];
    //FILE                    *fptr;
    int                     iiter, i, j, k;


    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("poisson_grid", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    del_u                               = allocate_1D_point_2D_array ("del_u", nx * ny);
    prev_grid                           = allocate_2D_point_2D_array ("prev_grid", nx, ny);
    pq                                  = allocate_2D_point_2D_array ("pq", nx, ny);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);


    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Solve for elliptic grid with control functions */
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid);   // Write initial grid to .dat files
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Set previous iteration (x, y) */
        equals_2D_point_2D_array (nx, ny, grid, &prev_grid);

        /* Grid derivatives and coefficients */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        /* Compute control functions (p, q) */
        compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &pq);
        //compute_steger_sorenson (nx, ny, dgrid, d2grid, &pq);
        if (iiter == 0)
            write_point_2D_to_vts ("initial_pq.vts", nx, ny, grid, "pq", pq);   // Write initial (p, q) to a Paraview file

        /* Compute new interior grid values */
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                g                       = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
                g1                      = coeffs[i][j].alpha/g;
                g2                      = (TWO * coeffs[i][j].beta)/g;
                g2                      = 0.0;
                g3                      = coeffs[i][j].gamma/g;
                g4                      = (coeffs[i][j].J * coeffs[i][j].J)/g;

                grid[i][j].x            = (g1 * ((grid[i + 1][j].x + grid[i - 1][j].x)/(dxi * dxi))) -
                                          (g2 * d2grid[i][j].xy_der.x) +
                                          (g3 * ((grid[i][j + 1].x + grid[i][j - 1].x)/(deta * deta))) +
                                          (g4 * ((pq[i][j].x * dgrid[i][j].x_der.x) + (pq[i][j].y * dgrid[i][j].y_der.x)));
                grid[i][j].y            = (g1 * ((grid[i + 1][j].y + grid[i - 1][j].y)/(dxi * dxi))) -
                                          (g2 * d2grid[i][j].xy_der.y) +
                                          (g3 * ((grid[i][j + 1].y + grid[i][j - 1].y)/(deta * deta))) +
                                          (g4 * ((pq[i][j].x * dgrid[i][j].x_der.y) + (pq[i][j].y * dgrid[i][j].y_der.y)));
            }
        }

        /* Under-relaxation */
        for (k = 0; k < nx * ny; k++)
        {
            i                           = k%nx;
            j                           = k/nx;
            grid[i][j].x                = (omega * grid[i][j].x) + ((ONE - omega) * prev_grid[i][j].x);
            grid[i][j].y                = (omega * grid[i][j].y) + ((ONE - omega) * prev_grid[i][j].y);

            /* Final iteration update vector */
            if (iiter == niter - 1)
            {
                del_u[k].x                  = grid[i][j].x - prev_grid[i][j].x;
                del_u[k].y                  = grid[i][j].y - prev_grid[i][j].y;
            }
        }
    }
    res_norm                                = compute_residual_norm (nx, ny, del_u); // Final iteration residual
    printf ("res_norm = %22.16LF\n", res_norm);
    snprintf (file1, 256, "final_pq_%d.vts", iiter);
    write_point_2D_to_vts (file1, nx, ny, grid, "pq", pq);  // Write final (p, q) to a Paraview file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files


    /* Free memory */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_1D_point_2D_array ("del_u", del_u);
    free_2D_point_2D_array ("prev_grid", nx, prev_grid);
    free_2D_point_2D_array ("pq", nx, pq);


    return grid;
}






/*
   Elliptic grid generation with non-zero control
   functions using pointwise iteration. Uses the
   (\phi, \psi) formulation of the control functions
   rather than (p, q)

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial (x, y) grid
                     niter      - number of iterations
*/
point_2D **poisson_grid_2D_point_alt (int nx, int ny, point_2D **init_grid, int niter)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, g, g1, g2, g3;
    const long double       omega = 0.001;
    point_2D                *del_u, **prev_grid, **phipsi;
    char                    file1[256];
    //FILE                    *fptr;
    int                     iiter, i, j, k;


    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("poisson_grid", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    del_u                               = allocate_1D_point_2D_array ("del_u", nx * ny);
    prev_grid                           = allocate_2D_point_2D_array ("prev_grid", nx, ny);
    phipsi                              = allocate_2D_point_2D_array ("phipsi", nx, ny);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);


    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Solve for elliptic grid with control functions */
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid);   // Write initial grid to .dat files
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Set previous iteration (x, y) */
        equals_2D_point_2D_array (nx, ny, grid, &prev_grid);

        /* Grid derivatives and coefficients */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        /* Compute control functions (phi, psi) */
        compute_thomas_middlecoff_alt (nx, ny, dgrid, d2grid, &phipsi);
        //compute_steger_sorenson (nx, ny, dgrid, d2grid, &phipsi);
        if (iiter == 0)
            write_point_2D_to_vts ("initial_pq.vts", nx, ny, grid, "phipsi", phipsi);   // Write initial CF to a Paraview file

        /* Compute new interior grid values */
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                g                       = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
                g1                      = coeffs[i][j].alpha/g;
                g2                      = (TWO * coeffs[i][j].beta)/g;
                g3                      = coeffs[i][j].gamma/g;

                grid[i][j].x            = (g1 * ((grid[i + 1][j].x + grid[i - 1][j].x)/(sqr(dxi)))) -
                                          (g2 * d2grid[i][j].xy_der.x) +
                                          (g3 * ((grid[i][j + 1].x + grid[i][j - 1].x)/(sqr(deta)))) +
                                          (coeffs[i][j].alpha * phipsi[i][j].x * dgrid[i][j].x_der.x) +
                                          (coeffs[i][j].gamma * phipsi[i][j].y * dgrid[i][j].y_der.x);
                grid[i][j].y            = (g1 * ((grid[i + 1][j].y + grid[i - 1][j].y)/(sqr(dxi)))) -
                                          (g2 * d2grid[i][j].xy_der.y) +
                                          (g3 * ((grid[i][j + 1].y + grid[i][j - 1].y)/(sqr(deta)))) +
                                          (coeffs[i][j].alpha * phipsi[i][j].x * dgrid[i][j].x_der.y) +
                                          (coeffs[i][j].gamma * phipsi[i][j].y * dgrid[i][j].y_der.y);
            }
        }

        /* Under-relaxation */
        for (k = 0; k < nx * ny; k++)
        {
            i                           = k%nx;
            j                           = k/nx;
            grid[i][j].x                = (omega * grid[i][j].x) + ((ONE - omega) * prev_grid[i][j].x);
            grid[i][j].y                = (omega * grid[i][j].y) + ((ONE - omega) * prev_grid[i][j].y);

            /* Update vector */
            del_u[k].x                  = grid[i][j].x - prev_grid[i][j].x;
            del_u[k].y                  = grid[i][j].y - prev_grid[i][j].y;
        }

        res_norm                                = compute_residual_norm (nx, ny, del_u); // Iteration residual
        printf ("Iteration %d res_norm = %22.16LF\n", iiter, res_norm);
    }
    snprintf (file1, 256, "final_pq_%d.vts", iiter);
    write_point_2D_to_vts (file1, nx, ny, grid, "phipsi", phipsi);  // Write final (phi, psi) to a Paraview file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files


    /* Free memory */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_1D_point_2D_array ("del_u", del_u);
    free_2D_point_2D_array ("prev_grid", nx, prev_grid);
    free_2D_point_2D_array ("phipsi", nx, phipsi);


    return grid;
}






point_2D **biharmonic_grid_2D (int nx, int ny, point_2D **init_grid, int niter)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    point_2D                **pq;
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid, **d2pq;
    long double             dxi, deta, **lhs_array, *rhs_array, *update, cxi, ceta, res_norm;
    const long double       temp1 = 1.0, temp2 = 0.0, temp3 = 1.0, temp4 = 0.0, damping = 1.0;
    bh_sol_2D               *fu, *del_u;
    bh_jacobian_2D          **dFdu;
    char                    file1[256];
    int                     iiter;


    cxi                                 = temp1 * ((long double) pow(10.0, (double) temp2));
    ceta                                = temp3 * ((long double) pow(10.0, (double) temp4));

    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("biharmonic_grid", nx, ny);
    pq                                  = allocate_2D_point_2D_array ("pq", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    d2pq                                = allocate_2D_grid_dder_2D_array ("d2pq", nx, ny);

    /* Allocate memory for Newton's method arrays */
    dFdu                                = allocate_2D_bh_jacobian_2D_array ("dFdu", nx * ny, nx * ny);
    fu                                  = allocate_1D_bh_sol_2D_array ("fu", nx * ny);
    lhs_array                           = allocate_2D_long_double_array ("rhs_array", 4 * nx * ny, 4 * nx * ny);
    rhs_array                           = allocate_1D_long_double_array ("rhs_array", 4 * nx * ny);
    update                              = allocate_1D_long_double_array ("update", 4 * nx * ny);
    del_u                               = allocate_1D_bh_sol_2D_array ("del_u", nx * ny);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);

    /* Compute initial (p, q), grid derivatives
       and 1st derivative related coefficients */
    initial_pq_2D (1, nx, ny, cxi, ceta, grid, &pq);
    initial_pq_2D_laplace (nx, ny, cxi, ceta, grid, &pq);
    write_point_2D_to_vts ("initial_pq.vts", nx, ny, grid, "pq", pq);
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);
    grid_second_ders_2D (nx, ny, grid, &d2grid);
    pq_second_ders_2D (nx, ny, pq, &d2pq);


    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));

    /* Solve for biharmonic grid */
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid); // Write initial grid to .dat files
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Newton's method */
        bh_construct_newton_matrix (nx, ny, dxi, deta, pq, dgrid, coeffs, d2grid, d2pq, cxi, ceta, &dFdu);
        bh_construct_newton_rhs_vector (nx, ny, grid, pq, dgrid, coeffs, d2grid, d2pq, cxi, ceta, &fu);
        convert_bh_jacobian_2D (nx * ny, dFdu, &lhs_array);
        convert_bh_sol_2D (nx * ny, fu, &rhs_array);
        LU_linear_system_solve (4 * nx * ny, lhs_array, rhs_array, &update);
        convert_to_bh_sol_2D (nx * ny, update, &del_u);

        /* Update (x, y) and (p, q) 2D arrays */
        bh_update_solution (nx, ny, damping, del_u, &grid, &pq);

        /* Compute new derivatives for next iteration */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);
        pq_second_ders_2D (nx, ny, pq, &d2pq);
    }
    snprintf (file1, 256, "final_pq_%d.vts", iiter);
    write_point_2D_to_vts (file1, nx, ny, grid, "pq", pq);  // Write final (p, q) to a Paraview file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files


    /* Free memory from grid arrays */
    free_2D_point_2D_array ("pq", nx, pq);
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_grid_dder_2D_array ("d2pq", nx, d2pq);

    /* Free memory from Newton's method arrays */
    free_2D_bh_jacobian_2D_array ("dFdu", nx * ny, dFdu);
    free_1D_bh_sol_2D_array ("F(u)", fu);
    free_2D_long_double_array ("lhs_array", nx * ny, lhs_array);
    free_1D_long_double_array ("rhs_array", rhs_array);
    free_1D_long_double_array ("update", update);
    free_1D_bh_sol_2D_array ("del_u", del_u);


    return grid;
}
