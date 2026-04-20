#ifndef TYPES_H_
#define TYPES_H_



/*
   Data type for 2D coordinates

   point_2D.x   - x coordinate
   point_2D.y   - y coordinate
*/
typedef struct
{
    long double             x;
    long double             y;
} point_2D;



/*
   Data type for 2D structured grids

   grid_2D.nblks        - number of blocks
   grid_2D.nx           - number of x points in each block
   grid_2D.ny           - number of y points in each block
   grid_2D.xi_bounds    - xi boundaries coordinates for each block
   grid_2D.eta_bounds   - eta boundaries coordinates for each block
   grid_2D.grid         - all grid coordinates
*/
typedef struct
{
    long                    nx;
    long                    ny;
    point_2D                **x_bounds;
    point_2D                **y_bounds;
    point_2D                **grid;
} grid_2D;



/*
   Data type to store the derivative of
   a grid point (x, y) with respect to the
   computational grid coordinates (xi, eta).
   Alternative to allocating a 3D point_2D
   array

   grid_der_2D.x_der    - derivative wrt "x" coordinate (xi)
   grid_der_2D.y_der    - derivative wrt "y" coordinate (eta)
*/
typedef struct
{
    point_2D                x_der;
    point_2D                y_der;
} grid_der_2D;



/*
   Data type to store various quantities
   computed using only first derivatives
   needed for the biharmonic equation
*/
typedef struct
{
    long double             alpha;  // dx_deta**2 + dy_deta**2
    long double             beta;   // dx_dxi * dx_deta + dy_dxi * dy_deta
    long double             gamma;  // dx_dxi**2 + dy_dxi**2
    long double             J;      // Jacobian = dx_dxi * dy_deta - dx_deta * dy_dxi
} coeffs_1;



/*
   Data type to store the second derivatives
   of a grid point (x, y) with respect to the
   computational grid coordinates (xi, eta)

   grid_dder_2D.x_der   - derivative wrt "x" coordinate (xi)
   grid_dder_2D.xy_der  - mixed derivative
   grid_dder_2D.y_der   - derivative wrt "y" coordinate (eta)
*/
typedef struct
{
    point_2D                x_der;
    point_2D                xy_der;
    point_2D                y_der;
} grid_dder_2D;



/*
   Data type to store all components of
   the Jacobian of a elliptic grid
   generation system. This is used to frame
   the Newton problem for elliptic grid
   generation
*/
typedef struct
{
    point_2D                ddx;
    point_2D                ddy;
} ell_jacobian_2D;



/*
   Data type to store all components
   of a biharmonic grid generation solution.
   This is used to frame the Newton problem
*/
typedef struct
{
    long double             x;
    long double             y;
    long double             p;
    long double             q;
} bh_sol_2D;



/*
   Data type to store all components of
   the Jacobian of a biharmonic grid
   generation system. This is used to frame
   the Newton problem for biharmonic grid
   generation
*/
typedef struct
{
    bh_sol_2D               ddx;
    bh_sol_2D               ddy;
    bh_sol_2D               ddp;
    bh_sol_2D               ddq;
} bh_jacobian_2D;



#endif  // TYPES_H_
