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

    TODO: Add checks for nx >= 3 and ny >= 3
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

    TODO: Add checks for nx >= 4 and ny >= 4
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
   Compute the Steger-Sorenson control functions on the (xi, eta)
   computational-grid boundaries, then TFI to the interior:

       phi = -(r_xi  . r_xixi)  / gamma  -  (r_xi  . r_etaeta) / alpha
       psi = -(r_eta . r_etaeta) / alpha  -  (r_eta . r_xixi)  / gamma

   These come from dotting the elliptic generation system

       alpha r_xixi - 2 beta r_xieta + gamma r_etaeta
                                    + alpha gamma (P r_xi + Q r_eta) = 0

   with r_xi and r_eta and assuming orthogonality (beta = 0). The
   beta = 0 assumption is a recipe-level approximation applied
   uniformly at every boundary point -- it is NOT specific to corners.
   Where the user-supplied boundaries are non-orthogonal, the formula
   is approximate everywhere on the boundary, not just at the corner.

   FD-error structure (verified empirically on a quartic test grid;
   see ~/Desktop/temp/ss_corner_pre_derivation.txt and the verification
   driver in the same directory):

     - phi gradients with respect to r_xixi and r_etaeta are both
       proportional to r_xi. The one-sided FD^2 error in r_xixi
       (occurring on i = 0 and i = nx-1 boundaries) is amplified,
       while the one-sided FD^2 error in r_etaeta (j = 0, ny-1) is
       largely suppressed by the dot-product structure. Result:
       i-edges dominate phi error; j-edges are quiet for phi.
     - By symmetry (psi gradients proportional to r_eta), j-edges
       dominate psi error and i-edges are quiet for psi.
     - Corners get one-sided FD^2 in both directions and are joint
       worst, but only marginally so vs. the dominant edge for each
       function -- they are NOT uniquely pathological.
     - Interior values are TFI of boundary values. The TFI
       non-separability error (TFI of analytical-boundary-phi minus
       analytical interior phi) is h-independent; for non-trivial
       geometries it can dominate over the O(h^2) FD-propagation
       term at modest resolutions.

   Input parameters: nx     - number of xi points
                     ny     - number of eta points
                     mode   - 0: divide by gamma, alpha to produce
                                 (P, Q) for the J^2-P solver convention
                              nonzero: produce (phi, psi) for the
                                 alpha-phi solver convention
                     coeffs - per-point alpha, beta, gamma, J;
                              accessed only when mode = 0
                     dgrid  - per-point first derivatives of (x, y)
                              with respect to (xi, eta)
                     d2grid - per-point second derivatives of (x, y)
                              with respect to (xi, eta)
   Output:           pq     - per-point (phi, psi) (mode != 0) or
                              (P, Q) (mode = 0) values
*/
void compute_steger_sorenson (int nx, int ny, int mode, coeffs_1 **coeffs, grid_der_2D **dgrid,
                              grid_dder_2D **d2grid, point_2D ***pq)
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
            (*pq)[i][j].y               = -(temp1/temp2) - (temp3/temp4);
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
            (*pq)[i][j].y               = -(temp1/temp2) - (temp3/temp4);
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


    /* Compute (P, Q) representation if needed */
    if (mode == 0)
    {
        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                (*pq)[i][j].x           = (*pq)[i][j].x/coeffs[i][j].gamma;
                (*pq)[i][j].y           = (*pq)[i][j].y/coeffs[i][j].alpha;
            }
        }
    }


    /* Free memory */
    free_2D_point_2D_array ("xieta", nx, xieta);
    free_2D_point_2D_array ("prev_pq", nx, prev_pq);
}






/*
   Compute the maximum per-point Euclidean magnitude of a
   point_2D array of length nx*ny:

       max_norm = max_i sqrt(v[i].x^2 + v[i].y^2)

   Two intended uses inside the elliptic and Poisson solvers:

     1. Newton step norm. Pass the Newton update del_u; the return
        value is the largest single grid-point movement this
        iteration. The caller divides by a reference length
        (bounding-box diagonal, mean spacing, etc) to form a
        dimensionless, grid-size-independent step tolerance --
        the primary stopping criterion for the Newton loop.

     2. Equation residual norm. Pass the Newton RHS fu; the return
        value is the largest per-point residual of the discretized
        equation. This is checked once after the Newton loop exits
        as an honest "did we actually solve it" diagnostic.

   The function itself is just a max-norm over a point_2D array;
   the meaning of the result depends on what the caller passes in

   Input parameters: nx  - number of x points
                     ny  - number of y points
                     v   - flat (length nx*ny) array of point_2D
                           values (typically Newton update del_u
                           or RHS fu)
*/
long double compute_residual_norm (int nx, int ny, point_2D *v)
{
    /* Return max_norm */
    long double             max_norm;

    /* Local variables */
    int                     i;
    long double             norm_i;


    /* Initialize residual norm as 0 */
    max_norm                            = ZERO;
    for (i = 0; i < nx * ny; i++)
    {
        norm_i = sqrtl ((v[i].x * v[i].x) + (v[i].y * v[i].y));
        if (norm_i > max_norm)
        {
            max_norm = norm_i;
        }
    }


    return max_norm;
}




/*
   Compute the bounding-box diagonal of a 2D grid:
       L_ref = sqrt((max_x - min_x)^2 + (max_y - min_y)^2)

   Used as a reference length to non-dimensionalize the Newton
   step norm in elliptic / Poisson solvers, so the convergence
   tolerance becomes "max grid-point movement as a fraction of
   domain size"

   Input parameters: nx   - number of x points
                     ny   - number of y points
                     grid - 2D array of point_2D values
*/
static long double compute_bbox_diagonal_2D (int nx, int ny, point_2D **grid)
{
    long double             min_x, max_x, min_y, max_y, dx, dy;
    int                     i, j;


    min_x = max_x = grid[0][0].x;
    min_y = max_y = grid[0][0].y;

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            if (grid[i][j].x < min_x) min_x = grid[i][j].x;
            if (grid[i][j].x > max_x) max_x = grid[i][j].x;
            if (grid[i][j].y < min_y) min_y = grid[i][j].y;
            if (grid[i][j].y > max_y) max_y = grid[i][j].y;
        }
    }

    dx = max_x - min_x;
    dy = max_y - min_y;
    return sqrtl ((dx * dx) + (dy * dy));
}




/*
   Tridiagonal solver via Thomas algorithm. Solves
       sub[i] * x[i-1] + diag[i] * x[i] + sup[i] * x[i+1] = rhs[i]
   for i = 0..n-1, with sub[0] and sup[n-1] required to be zero
   (boundary contributions should be folded into rhs by the caller).

   Two right-hand-sides supported simultaneously (rhs_x, rhs_y) — they
   share the same matrix factorization, so solving x and y components
   along a single grid line costs only one forward sweep.

   The arrays diag, rhs_x, rhs_y are MODIFIED in place by the forward
   sweep. Outputs are written to sol_x and sol_y. sub and sup are
   unchanged.

   Input/output parameters:
        n      - number of unknowns
        sub    - sub-diagonal of length n (sub[0] must be 0)
        diag   - diagonal of length n (overwritten)
        sup    - super-diagonal of length n (sup[n-1] must be 0)
        rhs_x  - x-component RHS, length n (overwritten)
        rhs_y  - y-component RHS, length n (overwritten)
        sol_x  - x-component solution, length n
        sol_y  - y-component solution, length n
*/
static void tridiag_solve_2rhs (int n, long double *sub, long double *diag, long double *sup,
                                long double *rhs_x, long double *rhs_y,
                                long double *sol_x, long double *sol_y)
{
    int                     i;
    long double             m;


    /* Forward sweep: eliminate sub-diagonal */
    for (i = 1; i < n; i++)
    {
        m                               = sub[i]/diag[i - 1];
        diag[i]                        -= m * sup[i - 1];
        rhs_x[i]                       -= m * rhs_x[i - 1];
        rhs_y[i]                       -= m * rhs_y[i - 1];
    }


    /* Back substitution */
    sol_x[n - 1]                        = rhs_x[n - 1]/diag[n - 1];
    sol_y[n - 1]                        = rhs_y[n - 1]/diag[n - 1];
    for (i = n - 2; i >= 0; i--)
    {
        sol_x[i]                        = (rhs_x[i] - (sup[i] * sol_x[i + 1]))/diag[i];
        sol_y[i]                        = (rhs_y[i] - (sup[i] * sol_y[i + 1]))/diag[i];
    }
}






/*
   Elliptic grid generation without any control
   functions using Newton's method

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial grid
                     niter      - number of Newton iterations
                     params     - solver tunables (uses tol_step, tol_resid)
*/
point_2D **elliptic_grid_2D (int nx, int ny, point_2D **init_grid, int niter, solver_params params)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, L_ref, fu_norm;
    point_2D                *fu, *del_u;
    ell_jacobian_2D         **dFdu;
    char                    file1[256];
    FILE                    *fptr;
    int                     iiter, i, j;


    /* Allocate memory for grid arrays */
    grid    = allocate_2D_point_2D_array ("elliptic_grid", nx, ny);
    dgrid   = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs  = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid  = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);

    /* Allocate Newton's method arrays — interior unknowns only.
       Boundaries are Dirichlet constants and do not enter the linear system. */
    const int n_int = (nx - 2) * (ny - 2);
    dFdu    = allocate_2D_ell_jacobian_2D_array ("dFdu", n_int, n_int);
    fu      = allocate_1D_point_2D_array ("fu", n_int);
    del_u   = allocate_1D_point_2D_array ("del_u", n_int);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);

    /* Computational grid spacing */
    dxi  = ONE/((long double) (nx - 1));
    deta = ONE/((long double) (ny - 1));

    /* Reference length for non-dimensional convergence test */
    L_ref = compute_bbox_diagonal_2D (nx, ny, grid);


    /* Solve for elliptic grids with no control functions */
    fptr = fopen ("residuals.dat", "w");
    write_2D_singleblock_vts ("initial_grid.vts", nx, ny, grid, 1); // Write initial grid to Paraview file
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid); // Write initial grid to .dat files
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Compute grid derivatives and first derivative coefficients */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        /* Newton's method */
        ell_construct_newton_matrix (nx, ny, dxi, deta, dgrid, coeffs, d2grid, &dFdu);
        ell_construct_newton_rhs_vector (nx, ny, grid, coeffs, d2grid, &fu, 0);
        gmres_linear_system_solve (0, nx, ny, dFdu, fu, &del_u);

        /* Step-norm convergence check over interior unknowns (relative to bbox diagonal) */
        res_norm = compute_residual_norm (nx - 2, ny - 2, del_u);
        fprintf (fptr, "%d  %22.16LE\n", iiter + 1, res_norm/L_ref);

        /* Update 2D (x, y) arrays */
        ell_update_solution (nx, ny, del_u, &grid);

        /* Break if Newton step is small relative to domain size */
        if (res_norm/L_ref < params.tol_step)
            break;
    }

    /* Post-loop equation-residual diagnostic */
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);
    grid_second_ders_2D (nx, ny, grid, &d2grid);
    ell_construct_newton_rhs_vector (nx, ny, grid, coeffs, d2grid, &fu, 0);
    /* Normalize per-point |F| by the local diagonal coefficient g so the
       residual is on the same scale as the step-norm (max-norm over
       interior). Without normalization F has alpha/dxi^2 magnitude and
       fires false-positive warnings near step-norm convergence.
       fu is sized (nx-2)*(ny-2); index via k = (i-1) + (j-1)*(nx-2). */
    fu_norm = ZERO;
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            long double         g_ij, norm_ij;
            int                 k = (i - 1) + ((j - 1) * (nx - 2));
            g_ij                        = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
            norm_ij                     = sqrtl ((fu[k].x * fu[k].x) + (fu[k].y * fu[k].y))/g_ij;
            if (norm_ij > fu_norm) fu_norm = norm_ij;
        }
    }
    if (fu_norm > params.tol_resid)
        fprintf (stderr, "Warning: elliptic_grid_2D equation residual %.6Le above tol %.6Le after %d iterations\n",
                 fu_norm, params.tol_resid, iiter + 1);

    snprintf (file1, 256, "final_grid_%d.vts", iiter + 1);
    write_2D_singleblock_vts (file1, nx, ny, grid, 1); // Write final grid to Paraview file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files
    fclose (fptr);


    /* Free memory from grid arrays */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);

    /* Free memory from Newton's method arrays (interior-only sized) */
    free_2D_ell_jacobian_2D_array ("dFdu", n_int, dFdu);
    free_1D_point_2D_array ("fu", fu);
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
                     params     - solver tunables (uses omega)
*/
point_2D **elliptic_grid_2D_point (int nx, int ny, point_2D **init_grid, int niter, solver_params params)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, g, g1, g2, g3;
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
            grid[i][j].x                = (params.omega * grid[i][j].x) + ((ONE - params.omega) * prev_grid[i][j].x);
            grid[i][j].y                = (params.omega * grid[i][j].y) + ((ONE - params.omega) * prev_grid[i][j].y);

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
                     params     - solver tunables (uses tol_step, tol_resid,
                                  ramp_iters, lambda_max)
*/
point_2D **poisson_grid_2D (int nx, int ny, point_2D **init_grid, int niter, solver_params params)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, L_ref, fu_norm;
    point_2D                *fu, *del_u, **pq, **prev_pq;
    ell_jacobian_2D         **dFdu;
    char                    file1[256];
    FILE                    *fptr;
    int                     iiter, i, j;
    long double             lambda;   /* TM->SS ramp parameter, 0..1 */


    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("poisson_grid", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    pq                                  = allocate_2D_point_2D_array ("pq", nx, ny);
    prev_pq                             = allocate_2D_point_2D_array ("prev_pq", nx, ny);

    /* Allocate Newton's method arrays — interior unknowns only.
       Boundaries are Dirichlet constants and do not enter the linear system. */
    const int n_int = (nx - 2) * (ny - 2);
    dFdu                                = allocate_2D_ell_jacobian_2D_array ("dFdu", n_int, n_int);
    fu                                  = allocate_1D_point_2D_array ("fu", n_int);
    del_u                               = allocate_1D_point_2D_array ("del_u", n_int);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);


    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Reference length for non-dimensional convergence test */
    L_ref                               = compute_bbox_diagonal_2D (nx, ny, grid);


    /* Solve for Poisson grid with control functions */
    fptr                                = fopen ("residuals.dat", "w"); // File for writing residuals
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid); // Write initial grid to .dat files
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Grid derivatives and coefficients */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        /* Compute control functions (p, q) by ramping from TM to SS:
              pq = (1 - lambda) * pq_TM + lambda * pq_SS
           lambda ramps linearly from 0 (pure TM) to 1 (pure SS) over the first
           ramp_iters outer iterations, then is capped at lambda_max. The cap
           keeps the iteration short of the full SS target, which is unreachable
           on some geometries (e.g. re-entrant corners). prev_pq is reused as
           scratch for pq_TM. */
        lambda                          = (iiter < params.ramp_iters)
                                          ? ((long double) iiter)/((long double) params.ramp_iters)
                                          : ONE;
        if (lambda > params.lambda_max) lambda = params.lambda_max;
        compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &prev_pq);                /* pq_TM */
        compute_steger_sorenson  (nx, ny, 0, coeffs, dgrid, d2grid, &pq);           /* pq_SS */
        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                pq[i][j].x              = ((ONE - lambda) * prev_pq[i][j].x) + (lambda * pq[i][j].x);
                pq[i][j].y              = ((ONE - lambda) * prev_pq[i][j].y) + (lambda * pq[i][j].y);
            }
        }
        if (iiter == 0)
        {
            write_2D_singleblock_vts ("initial_pq.vts", nx, ny, grid, 0);
            write_point_2D_to_vts ("initial_pq.vts", nx, ny, "pq", pq, 1); // Write initial (p, q) to a Paraview file
        }

        /* Newton's method */
        psn_construct_newton_matrix (nx, ny, dxi, deta, dgrid, coeffs, d2grid, pq, &dFdu);
        psn_construct_newton_rhs_vector (nx, ny, grid, dgrid, coeffs, d2grid, pq, &fu);
        gmres_linear_system_solve (0, nx, ny, dFdu, fu, &del_u);

        /* Step-norm convergence check over interior unknowns (relative to bbox diagonal) */
        res_norm                        = compute_residual_norm (nx - 2, ny - 2, del_u);
        fprintf (fptr, "%d  %22.16LE\n", iiter + 1, res_norm/L_ref);

        /* Update 2D (x, y) arrays */
        ell_update_solution (nx, ny, del_u, &grid);

        /* Break if Newton step is small relative to domain size */
        if (res_norm/L_ref < params.tol_step)
            break;
    }

    /* Post-loop equation-residual diagnostic. Use the SAME blended (P, Q)
       the iteration converged with (at the final, capped lambda). Computing
       residual against pure SS would give a misleading warning since the
       iteration's fixed point is the blended-target one, not pure SS. */
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);
    grid_second_ders_2D (nx, ny, grid, &d2grid);
    {
        long double lambda_final        = (iiter < params.ramp_iters)
                                          ? ((long double) iiter)/((long double) params.ramp_iters)
                                          : ONE;
        if (lambda_final > params.lambda_max) lambda_final = params.lambda_max;
        compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &prev_pq);
        compute_steger_sorenson  (nx, ny, 0, coeffs, dgrid, d2grid, &pq);
        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                pq[i][j].x              = ((ONE - lambda_final) * prev_pq[i][j].x) + (lambda_final * pq[i][j].x);
                pq[i][j].y              = ((ONE - lambda_final) * prev_pq[i][j].y) + (lambda_final * pq[i][j].y);
            }
        }
    }
    psn_construct_newton_rhs_vector (nx, ny, grid, dgrid, coeffs, d2grid, pq, &fu);
    /* Normalize per-point |F| by the local diagonal coefficient g so the
       residual is on the same scale as the step-norm (max-norm over
       interior). Without normalization F has alpha/dxi^2 magnitude and
       fires false-positive warnings near step-norm convergence.
       fu is sized (nx-2)*(ny-2); index via k = (i-1) + (j-1)*(nx-2). */
    fu_norm                             = ZERO;
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            long double         g_ij, norm_ij;
            int                 k = (i - 1) + ((j - 1) * (nx - 2));
            g_ij                        = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
            norm_ij                     = sqrtl ((fu[k].x * fu[k].x) + (fu[k].y * fu[k].y))/g_ij;
            if (norm_ij > fu_norm) fu_norm = norm_ij;
        }
    }
    if (fu_norm > params.tol_resid)
        fprintf (stderr, "Warning: poisson_grid_2D equation residual %.6Le above tol %.6Le after %d iterations\n",
                 fu_norm, params.tol_resid, iiter + 1);

    snprintf (file1, 256, "final_pq_%d.vts", iiter + 1);
    write_2D_singleblock_vts (file1, nx, ny, grid, 0);
    write_point_2D_to_vts (file1, nx, ny, "pq", pq, 1);  // Write final (p, q) to a Paraview file
    write_2D_singleblock_plot3D ("final_grid.x", nx, ny, grid); // Write final grid to PLOT3D file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files
    fclose (fptr);  // Close residual output file


    /* Free memory from grid arrays */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_2D_point_2D_array ("pq", nx, pq);
    free_2D_point_2D_array ("prev_pq", nx, prev_pq);

    /* Free memory from Newton's method arrays (interior-only sized) */
    free_2D_ell_jacobian_2D_array ("dFdu", n_int, dFdu);
    free_1D_point_2D_array ("fu", fu);
    free_1D_point_2D_array ("del_u", del_u);


    return grid;
}






/*
   Elliptic grid generation with Thomas-Middlecoff
   control functions using Gauss-Seidel iteration.
   Steger-Sorenson control functions should not
   be used here (see poisson_grid_2D_point_lim)

   Convergence: in-loop break when the max step-norm divided by
   the bounding-box diagonal falls below params.tol_step (default
   1e-8L). Post-loop equation residual is checked against
   params.tol_resid (default 1e-6L); a warning is emitted to
   stderr if the iteration did not reduce the residual below it.

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial (x, y) grid
                     niter      - maximum number of iterations
                     params     - solver tunables (uses omega, tol_step, tol_resid)
*/
point_2D **poisson_grid_2D_point (int nx, int ny, point_2D **init_grid, int niter, solver_params params)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, L_ref, fu_norm, g, g1, g2, g3, g4;
    point_2D                *del_u, **prev_grid, **pq;
    char                    file1[256];
    FILE                    *fptr;
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


    /* Reference length for non-dimensional convergence test */
    L_ref                               = compute_bbox_diagonal_2D (nx, ny, grid);


    /* Solve for Poisson grid via pointwise iteration */
    fptr                                = fopen ("residuals.dat", "w");
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
        if (iiter == 0)
        {
            write_2D_singleblock_vts ("initial_pq.vts", nx, ny, grid, 0);
            write_point_2D_to_vts ("initial_pq.vts", nx, ny, "pq", pq, 1);   // Write initial (p, q) to a Paraview file
        }

        /* Compute new interior grid values */
        for (i = 1; i < nx - 1; i++)
        {
            for (j = 1; j < ny - 1; j++)
            {
                g                       = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
                g1                      = coeffs[i][j].alpha/g;
                //g2                      = (TWO * coeffs[i][j].beta)/g;
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

        /* Under-relaxation and step-vector population */
        for (k = 0; k < nx * ny; k++)
        {
            i                           = k%nx;
            j                           = k/nx;
            grid[i][j].x                = (params.omega * grid[i][j].x) + ((ONE - params.omega) * prev_grid[i][j].x);
            grid[i][j].y                = (params.omega * grid[i][j].y) + ((ONE - params.omega) * prev_grid[i][j].y);
            del_u[k].x                  = grid[i][j].x - prev_grid[i][j].x;
            del_u[k].y                  = grid[i][j].y - prev_grid[i][j].y;
        }

        /* Step-norm convergence check (relative to bbox diagonal) */
        res_norm                        = compute_residual_norm (nx, ny, del_u);
        fprintf (fptr, "%d  %22.16LE\n", iiter + 1, res_norm/L_ref);

        /* Break if step is small relative to domain size */
        if (res_norm/L_ref < params.tol_step)
            break;
    }

    /* Post-loop equation-residual diagnostic. Evaluates the
       beta = 0 simplified Poisson system that the iteration
       actually solves:
           F = alpha r_xixi + gamma r_etaeta + J^2 (P r_xi + Q r_eta)
       Per-point F is divided by the iteration's diagonal coefficient
       g = 2 (alpha/dxi^2 + gamma/deta^2), yielding a residual whose
       scale matches the step-norm metric (since at each pointwise
       step |F| ~ g * |step|). Max-norm taken over interior points
       only; boundaries are Dirichlet and never moved by the iteration */
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);
    grid_second_ders_2D (nx, ny, grid, &d2grid);
    compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &pq);
    fu_norm                             = ZERO;
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            long double         F_x, F_y, src_x, src_y, norm_ij, g_ij;
            g_ij                        = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
            src_x                       = (pq[i][j].x * dgrid[i][j].x_der.x) +
                                          (pq[i][j].y * dgrid[i][j].y_der.x);
            src_y                       = (pq[i][j].x * dgrid[i][j].x_der.y) +
                                          (pq[i][j].y * dgrid[i][j].y_der.y);
            F_x                         = (coeffs[i][j].alpha * d2grid[i][j].x_der.x) +
                                          (coeffs[i][j].gamma * d2grid[i][j].y_der.x) +
                                          (coeffs[i][j].J * coeffs[i][j].J * src_x);
            F_y                         = (coeffs[i][j].alpha * d2grid[i][j].x_der.y) +
                                          (coeffs[i][j].gamma * d2grid[i][j].y_der.y) +
                                          (coeffs[i][j].J * coeffs[i][j].J * src_y);
            norm_ij                     = sqrtl ((F_x * F_x) + (F_y * F_y))/g_ij;
            if (norm_ij > fu_norm) fu_norm = norm_ij;
        }
    }
    if (fu_norm > params.tol_resid)
        fprintf (stderr, "Warning: poisson_grid_2D_point equation residual %.6Le above tol %.6Le after %d iterations\n",
                 fu_norm, params.tol_resid, iiter + 1);

    snprintf (file1, 256, "final_pq_%d.vts", iiter + 1);
    write_2D_singleblock_vts (file1, nx, ny, grid, 0);
    write_point_2D_to_vts (file1, nx, ny, "pq", pq, 1);  // Write final (p, q) to a Paraview file
    write_2D_singleblock_plot3D ("final_grid.x", nx, ny, grid); // Write final grid to PLOT3D file
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid); // Write final grid to .dat files
    fclose (fptr);


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
   Elliptic grid generation with non-zero control functions using
   line-implicit (LIM) iteration on xi-lines. Uses the (P, Q)
   formulation of the control functions with the J^2 P/Q convention.

   Each iteration, for each j (interior eta-line), build a tridiagonal
   system over i = 1..nx-2 with the xi-direction discretization
   IMPLICIT and the eta-direction terms EXPLICIT (using current
   grid values, with Gauss-Seidel ordering j = 1, 2, ..., ny-2).
   Solving the tridiagonal couples all xi-line unknowns simultaneously
   instead of one-at-a-time, which dramatically reduces the spectral
   radius of the iteration matrix when the source-term coupling is
   strong (as it is for SS control functions).

   beta = 0 is forced (matches the TM and SS recipe assumption + the
   convention in poisson_grid_2D_point). Cross-derivative term is
   dropped from the discretization.

   Convergence: in-loop break when max step-norm divided by the
   bounding-box diagonal falls below params.tol_step (default
   1e-8L). Post-loop equation residual checked against
   params.tol_resid (default 1e-6L); a warning is emitted to
   stderr if residual is above.

   Input parameters: nx         - number of x points
                     ny         - number of y points
                     init_grid  - initial (x, y) grid
                     niter      - maximum number of iterations
                     params     - solver tunables (uses tol_step, tol_resid,
                                  inner_niter, ramp_iters, lambda_max)
*/
point_2D **poisson_grid_2D_point_lim (int nx, int ny, point_2D **init_grid, int niter, solver_params params)
{
    /* Return grid */
    point_2D                **grid;

    /* Local variables */
    grid_der_2D             **dgrid;
    coeffs_1                **coeffs;
    grid_dder_2D            **d2grid;
    long double             dxi, deta, res_norm, L_ref, fu_norm;
    point_2D                *del_u, **prev_grid, **pq, **prev_pq;
    char                    file1[256];
    FILE                    *fptr;
    int                     iiter, i, j, k;

    /* Tridiagonal work arrays for xi-line solves (length nx-2) */
    long double             *tri_sub, *tri_diag, *tri_sup, *tri_rhs_x, *tri_rhs_y;
    long double             *tri_sol_x, *tri_sol_y;
    int                     nxi;

    int                     iiter_inner;
    long double             lambda;   /* TM->SS ramp parameter, 0..1 */


    nxi                                 = nx - 2; // Number of unknowns along a xi-line


    /* Allocate memory for grid arrays */
    grid                                = allocate_2D_point_2D_array ("poisson_grid", nx, ny);
    dgrid                               = allocate_2D_grid_der_2D_array ("dgrid", nx, ny);
    coeffs                              = allocate_2D_coeffs_1_array ("coeffs", nx, ny);
    d2grid                              = allocate_2D_grid_dder_2D_array ("d2grid", nx, ny);
    del_u                               = allocate_1D_point_2D_array ("del_u", nx * ny);
    prev_grid                           = allocate_2D_point_2D_array ("prev_grid", nx, ny);
    pq                                  = allocate_2D_point_2D_array ("pq", nx, ny);
    prev_pq                             = allocate_2D_point_2D_array ("prev_pq", nx, ny);

    /* Allocate tridiagonal work arrays */
    tri_sub                             = allocate_1D_long_double_array ("tri_sub",   nxi);
    tri_diag                            = allocate_1D_long_double_array ("tri_diag",  nxi);
    tri_sup                             = allocate_1D_long_double_array ("tri_sup",   nxi);
    tri_rhs_x                           = allocate_1D_long_double_array ("tri_rhs_x", nxi);
    tri_rhs_y                           = allocate_1D_long_double_array ("tri_rhs_y", nxi);
    tri_sol_x                           = allocate_1D_long_double_array ("tri_sol_x", nxi);
    tri_sol_y                           = allocate_1D_long_double_array ("tri_sol_y", nxi);


    /* Set grid array equal to initial grid */
    equals_2D_point_2D_array (nx, ny, init_grid, &grid);


    /* Computational grid spacing */
    dxi                                 = ONE/((long double) (nx - 1));
    deta                                = ONE/((long double) (ny - 1));


    /* Reference length for non-dimensional convergence test */
    L_ref                               = compute_bbox_diagonal_2D (nx, ny, grid);


    /* Solve for Poisson grid via xi-line LIM iteration */
    fptr                                = fopen ("residuals.dat", "w");
    write_2D_point_2D ("init_x.dat", "init_y.dat", nx, ny, grid);
    for (iiter = 0; iiter < niter; iiter++)
    {
        /* Save outer-start grid for outer step-norm */
        equals_2D_point_2D_array (nx, ny, grid, &prev_grid);

        /* Compute control functions (p, q) ONCE per outer iter, ramping from TM to SS:
              pq = (1 - lambda) * pq_TM + lambda * pq_SS
           where lambda ramps linearly from 0 (pure TM) to 1 (pure SS) over the first
           ramp_iters outer iterations. The ramp gives the iteration time to evolve
           the grid toward TM's gentler target before the more aggressive SS clustering
           kicks in. Once lambda reaches 1, the iteration runs pure SS for the
           remainder. Coefficients dgrid / coeffs / d2grid are recomputed inside the
           inner loop (they depend on the moving grid); the (P, Q) computed here
           stays FROZEN throughout the inner_niter sub-iterations.

           prev_pq is reused as scratch space for pq_TM here. */
        grid_first_ders_2D (nx, ny, grid, &dgrid);
        first_der_coefficients (nx, ny, dgrid, &coeffs);
        grid_second_ders_2D (nx, ny, grid, &d2grid);

        lambda                          = (iiter < params.ramp_iters)
                                          ? ((long double) iiter)/((long double) params.ramp_iters)
                                          : ONE;
        if (lambda > params.lambda_max) lambda = params.lambda_max;
        compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &prev_pq); // pq_TM
        compute_steger_sorenson  (nx, ny, 0, coeffs, dgrid, d2grid, &pq); // pq_SS
        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                pq[i][j].x              = ((ONE - lambda) * prev_pq[i][j].x) + (lambda * pq[i][j].x);
                pq[i][j].y              = ((ONE - lambda) * prev_pq[i][j].y) + (lambda * pq[i][j].y);
            }
        }
        if (iiter == 0)
        {
            write_2D_singleblock_vts ("initial_pq.vts", nx, ny, grid, 0);
            write_point_2D_to_vts ("initial_pq.vts", nx, ny, "pq", pq, 1);
        }

        /* Inner sub-iteration: inner_niter LIM steps with FROZEN (P, Q).
           Geometric coefficients (alpha, beta, gamma, J) ARE recomputed each
           inner step since they depend on the (moving) grid; only (P, Q) is
           held fixed. This is the iteration that the user-facing niter
           controls effectively up to (niter * inner_niter) total LIM steps. */
        for (iiter_inner = 0; iiter_inner < params.inner_niter; iiter_inner++)
        {
            /* Recompute geometric coefficients from current grid */
            grid_first_ders_2D (nx, ny, grid, &dgrid);
            first_der_coefficients (nx, ny, dgrid, &coeffs);
            grid_second_ders_2D (nx, ny, grid, &d2grid);

            /* xi-line LIM sweep: for each j, build and solve a tridiagonal
               in i = 1..nx-2. Ordering j = 1, 2, ..., ny-2 gives Gauss-Seidel
               coupling between adjacent eta-lines */
            for (j = 1; j < ny - 1; j++)
            {
                /* Build tridiagonal coefficients along the xi-line at this j.
                   Discretization (beta = 0): solving
                       alpha r_xixi + gamma r_etaeta + J^2 (P r_xi + Q r_eta) = 0
                   Collect r_{i,j} terms on the diagonal. For each i = 1..nx-2:
                       sub[idx] = -alpha/dxi^2 + J^2 P/(2 dxi)
                       diag[idx] = g = 2(alpha/dxi^2 + gamma/deta^2)
                       sup[idx] = -alpha/dxi^2 - J^2 P/(2 dxi)
                       rhs[idx] = (gamma/deta^2 - J^2 Q/(2 deta)) r_{i,j-1}
                                + (gamma/deta^2 + J^2 Q/(2 deta)) r_{i,j+1}     */
                for (i = 1; i < nx - 1; i++)
                {
                    int idx                 = i - 1;
                    long double Jsq         = coeffs[i][j].J * coeffs[i][j].J;
                    long double a_dxi2      = coeffs[i][j].alpha/(dxi * dxi);
                    long double g_deta2     = coeffs[i][j].gamma/(deta * deta);
                    long double JsqP_2dxi   = (Jsq * pq[i][j].x)/(TWO * dxi);
                    long double JsqQ_2deta  = (Jsq * pq[i][j].y)/(TWO * deta);
                    long double rhs_coef_jm = g_deta2 - JsqQ_2deta;
                    long double rhs_coef_jp = g_deta2 + JsqQ_2deta;

                    tri_sub[idx]            = -a_dxi2 + JsqP_2dxi;
                    tri_diag[idx]           = TWO * (a_dxi2 + g_deta2);
                    tri_sup[idx]            = -a_dxi2 - JsqP_2dxi;
                    tri_rhs_x[idx]          = (rhs_coef_jm * grid[i][j - 1].x) +
                                              (rhs_coef_jp * grid[i][j + 1].x);
                    tri_rhs_y[idx]          = (rhs_coef_jm * grid[i][j - 1].y) +
                                              (rhs_coef_jp * grid[i][j + 1].y);
                }

                /* Apply Dirichlet BCs at i=0 (folded into rhs[0]) and i=nx-1
                   (folded into rhs[nxi-1]). After folding, sub[0] = sup[nxi-1] = 0 */
                tri_rhs_x[0]               -= tri_sub[0]      * grid[0][j].x;
                tri_rhs_y[0]               -= tri_sub[0]      * grid[0][j].y;
                tri_sub[0]                  = 0.0L;
                tri_rhs_x[nxi - 1]         -= tri_sup[nxi - 1] * grid[nx - 1][j].x;
                tri_rhs_y[nxi - 1]         -= tri_sup[nxi - 1] * grid[nx - 1][j].y;
                tri_sup[nxi - 1]            = 0.0L;

                /* Solve tridiagonal for x and y components simultaneously */
                tridiag_solve_2rhs (nxi, tri_sub, tri_diag, tri_sup,
                                    tri_rhs_x, tri_rhs_y, tri_sol_x, tri_sol_y);

                /* Write solution back to grid (Gauss-Seidel: in-place update) */
                for (i = 1; i < nx - 1; i++)
                {
                    grid[i][j].x            = tri_sol_x[i - 1];
                    grid[i][j].y            = tri_sol_y[i - 1];
                }
            }
        }

        /* Outer step-vector population (grid - outer-start grid) */
        for (k = 0; k < nx * ny; k++)
        {
            i                           = k%nx;
            j                           = k/nx;
            del_u[k].x                  = grid[i][j].x - prev_grid[i][j].x;
            del_u[k].y                  = grid[i][j].y - prev_grid[i][j].y;
        }

        /* Step-norm convergence check (relative to bbox diagonal) */
        res_norm                        = compute_residual_norm (nx, ny, del_u);
        fprintf (fptr, "%d  %22.16LE\n", iiter + 1, res_norm/L_ref);

        /* Break if step is small relative to domain size */
        if (res_norm/L_ref < params.tol_step)
            break;
    }

    /* Post-loop equation-residual diagnostic (beta=0 simplified Poisson).
       Use the SAME blended (P, Q) the iteration solved with — at the final
       lambda (capped at params.lambda_max). Computing residual against
       pure SS would give a misleading warning since the iteration's
       fixed point is the blended-target one, not the pure-SS target. */
    grid_first_ders_2D (nx, ny, grid, &dgrid);
    first_der_coefficients (nx, ny, dgrid, &coeffs);
    grid_second_ders_2D (nx, ny, grid, &d2grid);
    {
        long double lambda_final        = (iiter < params.ramp_iters)
                                          ? ((long double) iiter)/((long double) params.ramp_iters)
                                          : ONE;
        if (lambda_final > params.lambda_max) lambda_final = params.lambda_max;
        compute_thomas_middlecoff (nx, ny, dgrid, d2grid, &prev_pq);
        compute_steger_sorenson  (nx, ny, 0, coeffs, dgrid, d2grid, &pq);
        for (i = 0; i < nx; i++)
        {
            for (j = 0; j < ny; j++)
            {
                pq[i][j].x              = ((ONE - lambda_final) * prev_pq[i][j].x) + (lambda_final * pq[i][j].x);
                pq[i][j].y              = ((ONE - lambda_final) * prev_pq[i][j].y) + (lambda_final * pq[i][j].y);
            }
        }
    }
    fu_norm                             = ZERO;
    for (i = 1; i < nx - 1; i++)
    {
        for (j = 1; j < ny - 1; j++)
        {
            long double         F_x, F_y, src_x, src_y, norm_ij, g_ij;
            g_ij                        = TWO * ((coeffs[i][j].alpha/(dxi * dxi)) +
                                                 (coeffs[i][j].gamma/(deta * deta)));
            src_x                       = (pq[i][j].x * dgrid[i][j].x_der.x) +
                                          (pq[i][j].y * dgrid[i][j].y_der.x);
            src_y                       = (pq[i][j].x * dgrid[i][j].x_der.y) +
                                          (pq[i][j].y * dgrid[i][j].y_der.y);
            F_x                         = (coeffs[i][j].alpha * d2grid[i][j].x_der.x) +
                                          (coeffs[i][j].gamma * d2grid[i][j].y_der.x) +
                                          (coeffs[i][j].J * coeffs[i][j].J * src_x);
            F_y                         = (coeffs[i][j].alpha * d2grid[i][j].x_der.y) +
                                          (coeffs[i][j].gamma * d2grid[i][j].y_der.y) +
                                          (coeffs[i][j].J * coeffs[i][j].J * src_y);
            norm_ij                     = sqrtl ((F_x * F_x) + (F_y * F_y))/g_ij;
            if (norm_ij > fu_norm) fu_norm = norm_ij;
        }
    }
    if (fu_norm > params.tol_resid)
        fprintf (stderr, "Warning: poisson_grid_2D_point_lim equation residual %.6Le above tol %.6Le after %d iterations\n",
                 fu_norm, params.tol_resid, iiter + 1);

    snprintf (file1, 256, "final_pq_%d.vts", iiter + 1);
    write_2D_singleblock_vts (file1, nx, ny, grid, 0);
    write_point_2D_to_vts (file1, nx, ny, "pq", pq, 1);
    write_2D_singleblock_plot3D ("final_grid.x", nx, ny, grid);
    write_2D_point_2D ("final_x.dat", "final_y.dat", nx, ny, grid);
    fclose (fptr);


    /* Free memory */
    free_2D_grid_der_2D_array ("dgrid", nx, dgrid);
    free_2D_coeffs_1_array ("coeffs", nx, coeffs);
    free_2D_grid_dder_2D_array ("d2grid", nx, d2grid);
    free_1D_point_2D_array ("del_u", del_u);
    free_2D_point_2D_array ("prev_grid", nx, prev_grid);
    free_2D_point_2D_array ("pq", nx, pq);
    free_2D_point_2D_array ("prev_pq", nx, prev_pq);
    free_1D_long_double_array ("tri_sub",   tri_sub);
    free_1D_long_double_array ("tri_diag",  tri_diag);
    free_1D_long_double_array ("tri_sup",   tri_sup);
    free_1D_long_double_array ("tri_rhs_x", tri_rhs_x);
    free_1D_long_double_array ("tri_rhs_y", tri_rhs_y);
    free_1D_long_double_array ("tri_sol_x", tri_sol_x);
    free_1D_long_double_array ("tri_sol_y", tri_sol_y);


    return grid;
}





/*
    Compute three grid quality metrics:

    - Orthogonality metric for each point as 
      |\beta|/sqrt(\alpha * \gamma) and statistics to 
      summarize grid quality post-run. The metric is 
      the cosine of the angle between r_\xi and 
      r_\eta: 0 implies full orthogonality and 1 
      implies a folded or coincident grid 
    - Aspect ratio: max(sqrt(\gamma/\alpha), 
      sqrt(\alpha/\gamma)) for each point. This should 
      always be >= 1, with 1 being isotropic. The value 
      for a degenerate point is set as 0 
    - Jacobian: minimum, maximum and number of folded 
      cells

    Input parameters: nx     - number of xi points 
                      ny     - number of eta points 
                      coeffs - array containing first derivative 
                               based coefficients 
*/
grid_quality_2D compute_grid_quality_2D (int nx, int ny, coeffs_1 **coeffs) 
{
    /* Return q */  
    grid_quality_2D         q;

    // Local variables 
    int                     i, j, is_interior;
    long double             ag, ortho_val, aspect_val, ratio, weight, 
                            sum_weighted, weight_sum, J; 


    // Allocate q component array 
    q.ortho = allocate_2D_long_double_array ("q.ortho", nx, ny);
    q.aspect_ratio = allocate_2D_long_double_array ("q.aspect_ratio", nx, ny);

    // Initialize diagnostic summary quantities
    q.ortho_max = ZERO;
    q.ortho_mean = ZERO; 
    q.n_degenerate = 0;

    q.aspect_ratio_max = ZERO;
    q.aspect_ratio_min = HUGE_VALL;

    q.J_min = HUGE_VALL;
    q.J_max = -HUGE_VALL;
    q.n_folded = 0;

    // Initialize local quantities 
    sum_weighted = ZERO;
    weight_sum = ZERO;


    // Loop over all points 
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            is_interior = (i > 0) && (i < nx - 1) && (j > 0) && (j < ny - 1);
            ag = coeffs[i][j].alpha * coeffs[i][j].gamma;
            J = coeffs[i][j].J;

            // Jacobian stats 
            if (is_interior)
            {
                if (J < q.J_min) q.J_min = J;
                if (J > q.J_max) q.J_max = J;
                if (J <= ZERO) q.n_folded += 1;
            }

            if (ag <= ZERO)     // Degenerate point
            {
                q.ortho[i][j] = ONE; 
                q.aspect_ratio[i][j] = ZERO;
                if (is_interior) q.n_degenerate += 1;
            }
            else
            {
                ortho_val = fabsl (coeffs[i][j].beta)/sqrtl (ag);

                // Symmetric aspect ratio 
                ratio = coeffs[i][j].gamma/coeffs[i][j].alpha;
                if (ratio < ONE) ratio = ONE/ratio;
                aspect_val = sqrtl (ratio);

                q.ortho[i][j] = ortho_val;
                q.aspect_ratio[i][j] = aspect_val;

                if (is_interior)
                {
                    weight = fabsl (J);
                    sum_weighted += ortho_val * weight;
                    weight_sum += weight;
                    if (ortho_val > q.ortho_max) q.ortho_max = ortho_val;
                    if (aspect_val > q.aspect_ratio_max) q.aspect_ratio_max = aspect_val;
                    if (aspect_val < q.aspect_ratio_min) q.aspect_ratio_min = aspect_val;
                }
            }
        }
    }

    if (weight_sum > ZERO) q.ortho_mean = sum_weighted/weight_sum;

    return q;
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
    write_2D_singleblock_vts ("initial_pq.vts", nx, ny, grid, 0);
    write_point_2D_to_vts ("initial_pq.vts", nx, ny, "pq", pq, 1);
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
    write_2D_singleblock_vts (file1, nx, ny, grid, 0);
    write_point_2D_to_vts (file1, nx, ny, "pq", pq, 1);  // Write final (p, q) to a Paraview file
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
