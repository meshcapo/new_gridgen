#ifndef GLOBVAR_H_
#define GLOBVAR_H_



extern char                 *mode;                                  // Either 2D or 3D
extern long                 nblocks;                                // Number of blocks in structured grid
extern long                 *npts_x, *npts_y, *npts_z;              // Number of points in x, y and z for all blocks
extern char                 *init_gridname;                         // Name of the file containing initial grid

static const long double    ZERO  = 0.0L;
static const long double    HALF  = 0.5L;
static const long double    ONE   = 1.0L;
static const long double    TWO   = 2.0L;
static const long double    THREE = 3.0L;
static const long double    FOUR  = 4.0L;
static const long double    FIVE  = 5.0L;
static const long double    SIX   = 6.0L;
static const long double    SEVEN = 7.0L;
static const long double    EIGHT = 8.0L;
static const long double    NINE  = 9.0L;
static const long double    TEN   = 10.0L;



#endif // GLOBVAR_H_
