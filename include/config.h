#ifndef CONFIG_H_
#define CONFIG_H_ 



/*
    Solver selection enum 
*/
typedef enum 
{
    SOLVER_UNKNOWN = 0,
    SOLVER_ELL,         // pointwise elliptic grid generator
    SOLVER_ELL_NEWTON,  // newton-based elliptic grid generator
    SOLVER_PSN_TM,      // pointwise Poisson grid generator with 
                        // Thomas-Middlecoff control functions
    SOLVER_PSN_LIM,     // pointwise Poisson grid generator with 
                        // blended Thomas-Middlecoff and Steger-Sorenson 
                        // control functions
    SOLVER_NEWTON,      // newton-based Poisson grid generator with 
                        // blended control functions
} solver_type;



/*
    Run-time configuration parsed from an input 
    file (.in). We track whether optional arguments 
    have been set or not - if not, defaults are 
    used and a warning is displayed 
*/
typedef struct 
{
    // Required arguments 
    char                    mode[8];                // "2D" or "3D"
    long                    niter;                  // max solver iterations
    long                    nblocks;                // number of grid blocks
    long                    *npts_x;                // number of points in x direction
    long                    *npts_y;                // number of points in y direction
    long                    *npts_z;                // number of points in z direction
    char                    init_grid_file[256];    // path to initial grid file
    solver_type             solver;                 // solver type

    // Optional arguments
    // Cap on SS contribution in blended Poisson 
    // control function: default 0.1/0.2 per solver
    long double             lambda_max;
    int                     lambda_max_set;

    // Outer iterations to ramp lambda from 0 to 
    // lambda_max: default 20/500 per solver
    int                     ramp_iters;
    int                     ramp_iters_set; 

    // Frozen (P, Q) sub-iterations per outer 
    // iteration: default 10
    int                     inner_niter; 
    int                     inner_niter_set; 

    // Under-relaxation factor: default 1.0
    long double             omega;
    int                     omega_set; 

    // In loop normalized residual tolerance: 
    // default 1E-8
    long double             tol_step;
    int                     tol_step_set; 

    // Post loop residual warning threshold:
    // default 1E-6
    long double             tol_resid;
    int                     tol_resid_set;
} run_args;



// Parse a gridexec input file into the run_args struct 
void read_config_file (char *, run_args *);

// Fill in defaults for any missing optional fields and 
// show warning 
void apply_solver_defaults (run_args *);

// Release dynamically allocated memory 
void free_run_args (run_args *);



#endif // CONFIG_H_
