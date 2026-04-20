#ifndef GRID_SOLVERS_H_
#define GRID_SOLVERS_H_
#include "types.h"



/*
   Compute an interpolated 2D algebraic grid given
   boundary points
*/
point_2D **algebraic_grid_2D (int, int, point_2D **, point_2D **);

/* Compute 2D elliptic grid */
point_2D **elliptic_grid_2D (int, int, point_2D **, int);
point_2D **elliptic_grid_2D_point (int, int, point_2D **, int);

/* Compute 2D Poisson grid */
point_2D **poisson_grid_2D (int, int, point_2D **, int);
point_2D **poisson_grid_2D_point (int, int, point_2D **, int);
point_2D **poisson_grid_2D_point_alt (int, int, point_2D **, int);

/* Compute 2D biharmonic grid */
point_2D **biharmonic_grid_2D (int, int, point_2D **, int);



#endif // GRID_SOLVERS_H_
