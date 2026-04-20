#ifndef GLOBVAR_H_
#define GLOBVAR_H_



extern char                 *mode;                                  // Either 2D or 3D
extern long                 nblocks;                                // Number of blocks in structured grid
extern long                 *npts_x, *npts_y, *npts_z;              // Number of points in x, y and z for all blocks
extern char                 *init_gridname;                         // Name of the file containing initial grid

extern long double          ZERO, HALF, ONE, TWO, THREE, FOUR, FIVE,
                            SIX, SEVEN, EIGHT, NINE, TEN;



#endif // GLOBVAR_H_
