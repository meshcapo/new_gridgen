#ifndef GEOMETRY_H_
#define GEOMETRY_H_
#include "types.h"



/*
   Get the points along the x boundaries (which map to
   xi = 0 and xi = 1) and y boundaries (which map to
   eta = 0 and eta = 1) from a given grid
*/
point_2D **get_x_bounds_2D (int, int, point_2D **);
point_2D **get_y_bounds_2D (int, int, point_2D **);

/*
   Functions to cluster points using either exponential
   or hyperbolic tangent clustering
*/
long double *exp_clustering (int, int, long double, long double, long double);
long double *tanh_clustering (int, int, long double, long double, long double);

/*
   Create boundaries for a geometry
*/
point_2D **create_x_bounds_2D (int);
point_2D **create_y_bounds_2D (int);



#endif // GEOMETRY_H_
