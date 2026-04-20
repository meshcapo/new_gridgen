#ifndef AUX_FUNCTIONS_H_
#define AUX_FUNCTIONS_H_
#include "types.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>


/*
   Dynamically allocate 1D arrays
*/
int *allocate_1D_int_array (char *, int);
long *allocate_1D_long_array (char *, int);
long double *allocate_1D_long_double_array (char *, int);
point_2D *allocate_1D_point_2D_array (char *, int);
bh_sol_2D *allocate_1D_bh_sol_2D_array (char *, int);

/*
   Free memory from 1D arrays
*/
void free_1D_int_array (char *, int *);
void free_1D_long_array (char *, long *);
void free_1D_long_double_array (char *, long double *);
void free_1D_point_2D_array (char *, point_2D *);
void free_1D_bh_sol_2D_array (char *, bh_sol_2D *);

/*
   Dynamically allocate 2D arrays
*/
int **allocate_2D_int_array (char *, int, int);
long double **allocate_2D_long_double_array (char *, int, int);
point_2D **allocate_2D_point_2D_array (char *, int, int);
grid_der_2D **allocate_2D_grid_der_2D_array (char *, int, int);
coeffs_1 **allocate_2D_coeffs_1_array (char *, int, int);
grid_dder_2D **allocate_2D_grid_dder_2D_array (char *, int, int);
ell_jacobian_2D **allocate_2D_ell_jacobian_2D_array (char *, int, int);
bh_jacobian_2D **allocate_2D_bh_jacobian_2D_array (char *, int, int);

/*
   Free memory from 2D arrays
*/
void free_2D_int_array (char *, int, int **);
void free_2D_long_double_array (char *, int, long double **);
void free_2D_point_2D_array (char *, int, point_2D **);
void free_2D_grid_der_2D_array (char *, int, grid_der_2D **);
void free_2D_coeffs_1_array (char *, int, coeffs_1 **);
void free_2D_grid_dder_2D_array (char *, int, grid_dder_2D **);
void free_2D_ell_jacobian_2D_array (char *, int, ell_jacobian_2D **);
void free_2D_bh_jacobian_2D_array (char *, int, bh_jacobian_2D **);

/* Unaligned functions */ 
int hash (int, int, int);

/*
   Function to read command line arguments
   The necessary inputs are stored in
   global variables.
*/
long *get_npts_x_arg_2D (char *[], long);
long *get_npts_y_arg_2D (char *[], long);
long *get_npts_x_arg_3D (char *[], long);
long *get_npts_y_arg_3D (char *[], long);
long *get_npts_z_arg_3D (char *[], long);
void read_inputs (int, char *[], char **, long *, long **, long **,
                  long **, char **, long *);

/*
   Arrays for converting either from 1D to
   2D point_2D array or vice versa
*/
point_2D *convert_array_to_column (int, int, point_2D **);
void convert_column_to_array (int, int, point_2D *, point_2D ***);

/* Type-generic macro for assigning a variable as 0 */
point_2D zero_point_2D (void);
grid_der_2D zero_grid_der_2D (void);
coeffs_1 zero_coeffs_1 (void);
grid_dder_2D zero_grid_dder_2D (void);
bh_sol_2D zero_bh_sol_2D (void);
ell_jacobian_2D zero_ell_jacobian_2D (void);
bh_jacobian_2D zero_bh_jacobian_2D (void);

#define zero(T) _Generic((T){0},            \
    point_2D: zero_point_2D,                \
    grid_der_2D: zero_grid_der_2D,          \
    coeffs_1: zero_coeffs_1,                \
    grid_dder_2D: zero_grid_dder_2D,        \
    bh_sol_2D: zero_bh_sol_2D,              \
    ell_jacobian_2D: zero_ell_jacobian_2D,  \
    bh_jacobian_2D: zero_bh_jacobian_2D     \
)()

/* Type-generic macro for checking equality */
int equal_point_2D (point_2D, point_2D);
int equal_grid_der_2D (grid_der_2D, grid_der_2D);
int equal_coeffs_1 (coeffs_1, coeffs_1);
int equal_grid_dder_2D (grid_dder_2D, grid_dder_2D);
int equal_bh_sol_2D (bh_sol_2D, bh_sol_2D);
int equal_ell_jacobian_2D (ell_jacobian_2D, ell_jacobian_2D);
int equal_bh_jacobian_2D (bh_jacobian_2D, bh_jacobian_2D);

#define equal(a, b) _Generic((a),           \
    point_2D: equal_point_2D,               \
    grid_der_2D: equal_grid_der_2D,         \
    coeffs_1: equal_coeffs_1,               \
    grid_dder_2D: equal_grid_dder_2D,       \
    bh_sol_2D: equal_bh_sol_2D,             \
    ell_jacobian_2D: equal_ell_jacobian_2D, \
    bh_jacobian_2D: equal_bh_jacobian_2D    \
)(a, b)


/*
   Array arithmetic functions
*/
point_2D *zero_1D_point_2D_array (char *, int);
point_2D **zero_2D_point_2D_array (char *, int, int);
point_2D *equals_1D_point_2D_array (int, point_2D *);
void equals_2D_point_2D_array (int, int, point_2D **, point_2D ***);
point_2D **sub_2D_point_2D_arrays (int, int, point_2D **, point_2D **);
long double array_max_2D_long_double (int, int, long double **);
point_2D array_max_2D_point_2D (int, int, point_2D **);
long double **zero_2D_long_double_array (char *, int, int);
void equals_zero_2D_long_double_array (int, int, long double ***);
void equals_zero_1D_point_2D_array (int, point_2D **);

/* Newton's method functions */
void convert_ell_jacobian_2D (int, ell_jacobian_2D **, long double ***);
void convert_point_2D (int, point_2D *, long double **);
void convert_to_point_2D (int, long double *, point_2D **);
void convert_bh_jacobian_2D (int, bh_jacobian_2D **, long double ***);
void convert_bh_sol_2D (int, bh_sol_2D *, long double **);
void convert_to_bh_sol_2D (int, long double *, bh_sol_2D **);
void std_matrix_LU (int, long double **, long double ***, long double ***);
void LU_linear_system_solve (int, long double **, long double *, long double **);
void gmres_linear_system_solve (int, int, int, ell_jacobian_2D **, point_2D *, point_2D **);

/*
   GSL functions 
*/
void convert_gsl_vector_2D (int, gsl_vector *, point_2D **);
void convert_to_gsl_vector_2D (int, point_2D *, gsl_vector **);
void convert_to_gsl_spmatrix_2D(int, ell_jacobian_2D **, gsl_spmatrix **);

/* Miscellaneous functions */ 
int hash (int, int, int);
point_2D set_as_equal_point_2D (point_2D);
void negative_point_2D (point_2D *);
point_2D divide_point_2D (point_2D, long double);
ell_jacobian_2D identity_ders_ell ();
int is_zero (bh_jacobian_2D);
void print_as_2D_matrix (int n, bh_jacobian_2D **);


#endif // AUX_FUNCTIONS_H_
