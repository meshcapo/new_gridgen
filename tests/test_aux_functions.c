#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <criterion/parameterized.h>
#include <criterion/logging.h>
#include <libgen.h>
#include <stdlib.h>
#include <unistd.h>
#include "aux_functions.h"
#include "types.h"



/* -------------------------------------------------------- */
/* Test-only 1D array wrappers to unify function signatures */
/* -------------------------------------------------------- */
static void *alloc_int_1D_wrapper (char *name, int n)
{
    return (void *) allocate_1D_int_array (name, n);
}
static void free_int_1D_wrapper (char *name, void *p)
{
    free_1D_int_array (name, (int *) p);
}

static void *alloc_long_1D_wrapper (char *name, int n)
{
    return (void *) allocate_1D_long_array (name, n);
}
static void free_long_1D_wrapper (char *name, void *p)
{
    free_1D_long_array (name, (long *) p);
}

static void *alloc_ld_1D_wrapper (char *name, int n)
{
    return (void *) allocate_1D_long_double_array (name, n);
}
static void free_ld_1D_wrapper (char *name, void *p)
{
    free_1D_long_double_array (name, (long double *) p);
}

static void *alloc_p2D_1D_wrapper (char *name, int n)
{
    return (void *) allocate_1D_point_2D_array (name, n);
}
static void free_p2D_1D_wrapper (char *name, void *p)
{
    free_1D_point_2D_array (name, (point_2D *) p);
}

static void *alloc_bhsol2D_1D_wrapper (char *name, int n)
{
    return (void *) allocate_1D_bh_sol_2D_array (name, n);
}
static void free_bhsol2D_1D_wrapper (char *name, void *p)
{
    free_1D_bh_sol_2D_array (name, (bh_sol_2D *) p);
}


/* Parameterized test for allocation and free routines for 1D arrays */
typedef struct  /* Data type for parameterized test */
{
    char case_name[32];
    int n;
    void *(*alloc_fn)(char *, int);
    void (*free_fn)(char *, void *);
    int type_tag;
} alloc_case_1D_t;

ParameterizedTestParameters (aux_functions_1D_alloc, alloc_free_param_1D)   /* Generate parameter set */
{
    static alloc_case_1D_t cases[] = 
    {
        {"int", 10, alloc_int_1D_wrapper, free_int_1D_wrapper, 1},
        {"long", 10, alloc_long_1D_wrapper, free_long_1D_wrapper, 2},
        {"long double", 10, alloc_ld_1D_wrapper, free_ld_1D_wrapper, 3},
        {"point_2D", 10, alloc_p2D_1D_wrapper, free_p2D_1D_wrapper, 4},
        {"bh_sol_2D", 10, alloc_bhsol2D_1D_wrapper, free_bhsol2D_1D_wrapper, 5},
    };
    return cr_make_param_array (alloc_case_1D_t, cases, sizeof(cases)/sizeof(*cases), NULL);
}

ParameterizedTest (alloc_case_1D_t *param, aux_functions_1D_alloc, alloc_free_param_1D)
{
    // Test allocation
    void *arr = param->alloc_fn ("arr", param->n);
    cr_assert_not_null (arr, "Allocation failed for %s\n", param->case_name);

    // Determine if allocated memory is usable
    if (param->type_tag == 1) 
    {
        int *a = (int *) arr;
        a[0]= 42;
        cr_assert_eq (a[0], 42, "int 1D array write/read check failed for %s\n", param->case_name);
    }
    else if (param->type_tag == 2)
    { 
        long *a = (long *) arr;
        a[0] = 0xDEADBEEF;
        cr_assert_eq (a[0], 0xDEADBEEF, "long 1D array write/read check failed for %s\n", param->case_name);
    }
    else if (param->type_tag == 3)
    { 
        long double *a = (long double *) arr;
        a[0] = 3.1415L;
        cr_assert_float_eq ((double) a[0], 3.1415, 1e-12, "long double 1D array write/read check failed for %s\n", param->case_name);
    }
    else if (param->type_tag == 4)
    {
        point_2D *a = (point_2D *) arr;
        a[0].x = 0.0L;
        a[0].y = 0.0L;
        cr_assert_float_eq ((double) a[0].x, 0.0, 1e-12, "point_2D 1D array write/read check failed for %s.%c\n", param->case_name, 'x');
        cr_assert_float_eq ((double) a[0].y, 0.0, 1e-12, "point_2D 1D array write/read check failed for %s.%c\n", param->case_name, 'y');
    }
    else if (param->type_tag == 5)
    {
        bh_sol_2D *a = (bh_sol_2D *) arr;
        a[0].x = 0.0L;
        a[0].y = 0.0L;
        a[0].p = 0.0L;
        a[0].q = 0.0L;
        cr_assert_float_eq ((double) a[0].x, 0.0, 1e-12, "bh_sol_2D 1D array write/read check failed for %s.%c\n", param->case_name, 'x');
        cr_assert_float_eq ((double) a[0].y, 0.0, 1e-12, "bh_sol_2D 1D array write/read check failed for %s.%c\n", param->case_name, 'y');
        cr_assert_float_eq ((double) a[0].p, 0.0, 1e-12, "bh_sol_2D 1D array write/read check failed for %s.%c\n", param->case_name, 'p');
        cr_assert_float_eq ((double) a[0].q, 0.0, 1e-12, "bh_sol_2D 1D array write/read check failed for %s.%c\n", param->case_name, 'q');
    }

    // Test deallocation
    param->free_fn ("arr", arr); /* Shouldn't crash */
    arr = NULL;
    cr_assert_null (arr);
}



/* -------------------------------------------------------- */
/* Test-only 2D array wrappers to unify function signatures */
/* -------------------------------------------------------- */
static void *alloc_int_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_int_array (name, n1, n2);
}
static void free_int_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_int_array (name, n1, (int **) p);
}

static void *alloc_ld_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_long_double_array (name, n1, n2);
}
static void free_ld_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_long_double_array (name, n1, (long double **) p);
}

static void *alloc_p2D_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_point_2D_array (name, n1, n2);
}
static void free_p2D_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_point_2D_array (name, n1, (point_2D **) p);
}

static void *alloc_gd2D_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_grid_der_2D_array (name, n1, n2);
}
static void free_gd2D_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_grid_der_2D_array (name, n1, (grid_der_2D **) p);
}

static void *alloc_coeffs_1_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_coeffs_1_array (name, n1, n2);
}
static void free_coeffs_1_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_coeffs_1_array (name, n1, (coeffs_1 **) p);
}

static void *alloc_gdd2D_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_grid_dder_2D_array (name, n1, n2);
}
static void free_gdd2D_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_grid_dder_2D_array (name, n1, (grid_dder_2D **) p);
}

static void *alloc_ellj2D_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_ell_jacobian_2D_array (name, n1, n2);
}
static void free_ellj2D_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_ell_jacobian_2D_array (name, n1, (ell_jacobian_2D **) p);
}

static void *alloc_bhj2D_2D_wrapper (char *name, int n1, int n2)
{
    return (void *) allocate_2D_bh_jacobian_2D_array (name, n1, n2);
}
static void free_bhj2D_2D_wrapper (char *name, int n1, void *p)
{
    free_2D_bh_jacobian_2D_array (name, n1, (bh_jacobian_2D **) p);
}


/* Parameterized test for allocation and free routines for 2D arrays */
typedef struct  /* Data type for parameterized test */
{
    char case_name[32];
    int n1;
    int n2;
    void *(*alloc_fn)(char *, int, int);
    void (*free_fn)(char *, int, void *);
    int type_tag;
} alloc_case_2D_t;

ParameterizedTestParameters (aux_functions_2D_alloc, alloc_free_param_2D)   /* Generate parameter set */
{
    static alloc_case_2D_t cases[] = 
    {
        {"int", 10, 10, alloc_int_2D_wrapper, free_int_2D_wrapper, 1},
        {"long double", 10, 10, alloc_ld_2D_wrapper, free_ld_2D_wrapper, 2},
        {"point_2D", 10, 10, alloc_p2D_2D_wrapper, free_p2D_2D_wrapper, 3},
        {"grid_der_2D", 10, 10, alloc_gd2D_2D_wrapper, free_gd2D_2D_wrapper, 4},
        {"coeffs_1", 10, 10, alloc_coeffs_1_2D_wrapper, free_coeffs_1_2D_wrapper, 5},
        {"grid_dder_2D", 10, 10, alloc_gdd2D_2D_wrapper, free_gdd2D_2D_wrapper, 6},
        {"ell_jacobian_2D", 10, 10, alloc_ellj2D_2D_wrapper, free_ellj2D_2D_wrapper, 7},
        {"bh_jacobian_2D", 10, 10, alloc_bhj2D_2D_wrapper, free_bhj2D_2D_wrapper, 8},
    };
    return cr_make_param_array (alloc_case_2D_t, cases, sizeof(cases)/sizeof(*cases), NULL);
}

ParameterizedTest (alloc_case_2D_t *param, aux_functions_2D_alloc, alloc_free_param_2D)
{
    // Test allocation
    void *arr = param->alloc_fn ("arr", param->n1, param->n2);
    cr_assert_not_null (arr, "Allocation failed for %s\n", param->case_name);

    // Determine if allocated memory is usable
    if (param->type_tag == 1) 
    {
        int **a = (int **) arr;
        a[0][0]= 42;
        cr_assert_eq (a[0][0], 42, "int 2D array write/read check failed for %s\n", param->case_name);
    }
    else if (param->type_tag == 2)
    { 
        long double **a = (long double **) arr;
        a[0][0] = 0.0L;
        cr_assert_float_eq ((double) a[0][0], 0.0, 1e-12, "long double 1D array write/read check failed for %s\n", param->case_name);
    }
    else if (param->type_tag == 3)
    { 
        point_2D **a = (point_2D **) arr;
        a[0][0].x = 0.0L;
        a[0][0].y = 0.0L;
        cr_assert_float_eq ((double) a[0][0].x, 0.0, 1e-12, "long double 1D array write/read check failed for %s.%c\n", param->case_name, 'x');
        cr_assert_float_eq ((double) a[0][0].y, 0.0, 1e-12, "long double 1D array write/read check failed for %s.%c\n", param->case_name, 'y');
    }
    else if (param->type_tag == 4)
    {
        grid_der_2D **a = (grid_der_2D **) arr;
        a[0][0].x_der.x = 0.0L;
        a[0][0].x_der.y = 0.0L;
        a[0][0].y_der.x = 0.0L;
        a[0][0].y_der.y = 0.0L;
        cr_assert_float_eq ((double) a[0][0].x_der.x, 0.0, 1e-12, "grid_der_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "x_der", 'x');
        cr_assert_float_eq ((double) a[0][0].x_der.y, 0.0, 1e-12, "grid_der_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "x_der", 'y');
        cr_assert_float_eq ((double) a[0][0].y_der.x, 0.0, 1e-12, "grid_der_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "y_der", 'x');
        cr_assert_float_eq ((double) a[0][0].y_der.y, 0.0, 1e-12, "grid_der_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "y_der", 'y');
    }
    else if (param->type_tag == 5)
    {
        coeffs_1 **a = (coeffs_1 **) arr;
        a[0][0].alpha = 0.0L;
        a[0][0].beta  = 0.0L;
        a[0][0].gamma = 0.0L;
        a[0][0].J     = 0.0L;
        cr_assert_float_eq ((double) a[0][0].alpha, 0.0, 1e-12, "coeffs_1 2D array write/read check failed for %s.%s\n", param->case_name, "alpha");
        cr_assert_float_eq ((double) a[0][0].beta, 0.0, 1e-12, "coeffs_1 2D array write/read check failed for %s.%s\n", param->case_name, "beta");
        cr_assert_float_eq ((double) a[0][0].gamma, 0.0, 1e-12, "coeffs_1 2D array write/read check failed for %s.%s\n", param->case_name, "gamma");
        cr_assert_float_eq ((double) a[0][0].J, 0.0, 1e-12, "coeffs_1 2D array write/read check failed for %s.%s\n", param->case_name, "J");
    }
    else if (param->type_tag == 6)
    {
        grid_dder_2D **a = (grid_dder_2D **) arr;
        a[0][0].x_der.x = 0.0L;
        a[0][0].x_der.y = 0.0L;
        a[0][0].xy_der.x = 0.0L;
        a[0][0].xy_der.y = 0.0L;
        a[0][0].y_der.x = 0.0L;
        a[0][0].y_der.y = 0.0L;
        cr_assert_float_eq ((double) a[0][0].x_der.x, 0.0, 1e-12, "grid_dder_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "x_der", 'x');
        cr_assert_float_eq ((double) a[0][0].x_der.y, 0.0, 1e-12, "grid_dder_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "x_der", 'y');
        cr_assert_float_eq ((double) a[0][0].xy_der.x, 0.0, 1e-12, "grid_dder_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "xy_der", 'x');
        cr_assert_float_eq ((double) a[0][0].xy_der.y, 0.0, 1e-12, "grid_dder_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "xy_der", 'y');
        cr_assert_float_eq ((double) a[0][0].y_der.x, 0.0, 1e-12, "grid_dder_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "y_der", 'x');
        cr_assert_float_eq ((double) a[0][0].y_der.y, 0.0, 1e-12, "grid_dder_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "y_der", 'y');
    }
    else if (param->type_tag == 7)
    {
        ell_jacobian_2D **a = (ell_jacobian_2D **) arr;
        a[0][0].ddx.x = 0.0L;
        a[0][0].ddx.y = 0.0L;
        a[0][0].ddy.x = 0.0L;
        a[0][0].ddy.y = 0.0L;
        cr_assert_float_eq ((double) a[0][0].ddx.x, 0.0, 1e-12, "ell_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddx", 'x');
        cr_assert_float_eq ((double) a[0][0].ddx.y, 0.0, 1e-12, "ell_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddx", 'y');
        cr_assert_float_eq ((double) a[0][0].ddy.x, 0.0, 1e-12, "ell_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddy", 'x');
        cr_assert_float_eq ((double) a[0][0].ddy.y, 0.0, 1e-12, "ell_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddy", 'y');
    }
    else if (param->type_tag == 8)
    {
        bh_jacobian_2D **a = (bh_jacobian_2D **) arr;
        a[0][0].ddx.x = 0.0L;
        a[0][0].ddx.y = 0.0L;
        a[0][0].ddy.x = 0.0L;
        a[0][0].ddy.y = 0.0L;
        a[0][0].ddp.x = 0.0L;
        a[0][0].ddp.y = 0.0L;
        a[0][0].ddq.x = 0.0L;
        a[0][0].ddq.y = 0.0L;
        cr_assert_float_eq ((double) a[0][0].ddx.x, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddx", 'x');
        cr_assert_float_eq ((double) a[0][0].ddx.y, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddx", 'y');
        cr_assert_float_eq ((double) a[0][0].ddy.x, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddy", 'x');
        cr_assert_float_eq ((double) a[0][0].ddy.y, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddy", 'y');
        cr_assert_float_eq ((double) a[0][0].ddp.x, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddp", 'x');
        cr_assert_float_eq ((double) a[0][0].ddp.y, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddp", 'y');
        cr_assert_float_eq ((double) a[0][0].ddq.x, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddq", 'x');
        cr_assert_float_eq ((double) a[0][0].ddq.y, 0.0, 1e-12, "bh_jacobian_2D 2D array write/read check failed for %s.%s.%c\n", param->case_name, "ddq", 'y');
    }

    // Test deallocation
    param->free_fn ("arr", 10, arr); /* Shouldn't crash */
    arr = NULL;
    cr_assert_null (arr);
}


/* Test for read_inputs: 2D mode without initial grid file */
Test(aux_functions, read_inputs_no_file)
{
    // Synthetic argv mimicking: gridexec 2D 1 121 31 100
    char *argv[] = {"gridexec", "2D", "1", "121", "31", "100"};
    int argc = 6;

    // Output variables
    char *mode = "";
    char *filename = "";
    long nblks = 0;
    long *nx = NULL;
    long *ny = NULL;
    long *nz = NULL;
    long niter = 0;

    // Call function
    read_inputs(argc, argv, &mode, &nblks, &nx, &ny, &nz, &filename, &niter);

    // Verify parsed values
    cr_assert_str_eq(mode, "2D", "mode: expected '2D', got '%s'", mode);
    cr_assert_eq(nblks, 1L, "nblks: expected 1, got %ld", nblks);
    cr_assert_not_null(nx, "nx is NULL");
    cr_assert_eq(nx[0], 121L, "nx[0]: expected 121, got %ld", nx[0]);
    cr_assert_not_null(ny, "ny is NULL");
    cr_assert_eq(ny[0], 31L, "ny[0]: expected 31, got %ld", ny[0]);
    cr_assert_null(nz, "nz expected NULL");
    cr_assert_eq(niter, 100L, "niter: expected 100, got %ld", niter);

    // Clean up
    free_1D_long_array("nx", nx);
    free_1D_long_array("ny", ny);
}


/* ------------------------------------------------------ */ 
/* Verifier functions for zero() _Generic macro test      */
/* ------------------------------------------------------ */ 
static int verify_zero_point_2D (void)
{
    point_2D r = zero (point_2D); 
    return r.x == 0.0L && r.y == 0.0L;
}

static int verify_zero_grid_der_2D (void)
{
    grid_der_2D r = zero (grid_der_2D);
    return r.x_der.x == 0.0L && r.x_der.y == 0.0L &&
           r.y_der.x == 0.0L && r.y_der.y == 0.0L;
}

static int verify_zero_coeffs_1 (void)
{
    coeffs_1 r = zero (coeffs_1);
    return r.alpha == 0.0L && r.beta == 0.0L &&
           r.gamma == 0.0L && r.J == 0.0L;
}

static int verify_zero_grid_dder_2D (void)
{
    grid_dder_2D r = zero (grid_dder_2D);
    return r.x_der.x == 0.0L && r.x_der.y == 0.0L &&
           r.xy_der.x == 0.0L && r.xy_der.y == 0.0L &&
           r.y_der.x == 0.0L && r.y_der.y == 0.0L;
}

static int verify_zero_bh_sol_2D (void)
{
    bh_sol_2D r = zero (bh_sol_2D);
    return r.x == 0.0L && r.y == 0.0L &&
           r.p == 0.0L && r.q == 0.0L;
}

static int verify_zero_ell_jacobian_2D (void)
{
    ell_jacobian_2D r = zero (ell_jacobian_2D);
    return r.ddx.x == 0.0L && r.ddx.y == 0.0L &&
           r.ddy.x == 0.0L && r.ddy.y == 0.0L;
}

/* Parameterized test for zero() _Generic macro */
typedef struct
{
    char type_name[32];
    int (*verify)(void);
} zero_case_t;

ParameterizedTestParameters (aux_functions_zero, zero_macro)
{
    static zero_case_t cases[] = 
    {
        {"point_2D", verify_zero_point_2D},
        {"grid_der_2D", verify_zero_grid_der_2D},
        {"coeffs_1", verify_zero_coeffs_1},
        {"grid_dder_2D", verify_zero_grid_dder_2D},
        {"bh_sol_2D", verify_zero_bh_sol_2D},
        {"ell_jacobian_2D", verify_zero_ell_jacobian_2D},
    };

    return cr_make_param_array (zero_case_t, cases, sizeof(cases)/sizeof(*cases), NULL);
}

ParameterizedTest (zero_case_t *param, aux_functions_zero, zero_macro)
{
    cr_assert (param->verify(), "zero(%s) returned non-zero field(s)", param->type_name);
}

/* Test for function is_zero - check whether a bh_jacobian_2D 
   variable is zero-valued or not */
Test(aux_functions, is_zero)
{
    bh_jacobian_2D a = zero (bh_jacobian_2D);   /* Zero-valued input - expect 1 */ 
    bh_jacobian_2D b = zero (bh_jacobian_2D);   /* Non-zero input - expect 0 */
    b.ddx.x = 1.0L;

    cr_assert_eq (is_zero (a), 1, "is_zero: expected 1 for zero input, got 0");
    cr_assert_eq (is_zero (b), 0, "is_zero: expected 0 for non-zero input, got 1");
}


/* ------------------------------------------------------ */ 
/* Verifier functions for equal() _Generic macro test      */
/* ------------------------------------------------------ */ 
static int verify_equal_point_2D (void)
{
    point_2D a = zero (point_2D);
    point_2D b = zero (point_2D);
    point_2D c = {1.0L, 0.0L};
    return equal (a, b) == 1 && equal (a, c) == 0;
}

static int verify_equal_grid_der_2D (void)
{
    grid_der_2D a = zero (grid_der_2D);
    grid_der_2D b = zero (grid_der_2D);
    grid_der_2D c = zero (grid_der_2D);
    c.x_der.x = 1.0L;
    return equal (a, b) == 1 && equal (a, c) == 0;
}

static int verify_equal_coeffs_1 (void)
{
    coeffs_1 a = zero (coeffs_1);
    coeffs_1 b = zero (coeffs_1);
    coeffs_1 c = zero (coeffs_1);
    c.alpha = 1.0L;
    return equal (a, b) == 1 && equal (a, c) == 0;
}

static int verify_equal_grid_dder_2D (void)
{
    grid_dder_2D a = zero (grid_dder_2D);
    grid_dder_2D b = zero (grid_dder_2D);
    grid_dder_2D c = zero (grid_dder_2D);
    c.x_der.x = 1.0L;
    return equal (a, b) == 1 && equal (a, c) == 0;
}

static int verify_equal_bh_sol_2D (void)
{
    bh_sol_2D a = zero (bh_sol_2D);
    bh_sol_2D b = zero (bh_sol_2D);
    bh_sol_2D c = zero (bh_sol_2D);
    c.x = 1.0L;
    return equal (a, b) == 1 && equal (a, c) == 0;
}

static int verify_equal_ell_jacobian_2D (void)
{
    ell_jacobian_2D a = zero (ell_jacobian_2D);
    ell_jacobian_2D b = zero (ell_jacobian_2D);
    ell_jacobian_2D c = zero (ell_jacobian_2D);
    c.ddx.x = 1.0L;
    return equal (a, b) == 1 && equal (a, c) == 0;
}

static int verify_equal_bh_jacobian_2D (void)
{
    bh_jacobian_2D a = zero (bh_jacobian_2D);
    bh_jacobian_2D b = zero (bh_jacobian_2D);
    bh_jacobian_2D c = zero (bh_jacobian_2D);
    c.ddx.x = 1.0L;
    return equal (a, b) == 1 && equal (a, c) == 0;
}

/* Parameterized test for equal() _Generic macro */
typedef struct
{
    char type_name[32];
    int (*verify)(void);
} equal_case_t;

ParameterizedTestParameters (aux_functions_equal, equal_macro)
{
    static equal_case_t cases[] = 
    {
        {"point_2D", verify_equal_point_2D},
        {"grid_der_2D", verify_equal_grid_der_2D},
        {"coeffs_1", verify_equal_coeffs_1},
        {"grid_dder_2D", verify_equal_grid_dder_2D},
        {"bh_sol_2D", verify_equal_bh_sol_2D},
        {"ell_jacobian_2D", verify_equal_ell_jacobian_2D},
        {"bh_jacobian_2D", verify_equal_bh_jacobian_2D},
    };

    return cr_make_param_array (equal_case_t, cases, sizeof(cases)/sizeof(*cases), NULL);
}

ParameterizedTest (equal_case_t *param, aux_functions_equal, equal_macro)
{
    cr_assert (param->verify(), "equal(%s) returned incorrect result", param->type_name);
}


/* ------------------------------------------------------ */ 
/* Tests for miscellaneous functions                      */
/* ------------------------------------------------------ */ 

/* Test for set_as_equal_point_2D */ 
Test(aux_functions, set_as_equal_point_2D)
{
    // Local variable
    point_2D a = {1.5L, -3.5L}; 
    point_2D b = set_as_equal_point_2D (a);

    // Check for equality
    cr_assert_eq (a.x, b.x, "set_as_equal_point_2D failed: .x expected 1.5");
    cr_assert_eq (a.y, b.y, "set_as_equal_point_2D failed: .y expected -3.5");
}

/* Test for divide_point_2D_zero - trigger failure */
Test(aux_functions, divide_point_2D_zero, .exit_code = EXIT_FAILURE)
{
    point_2D a =  {1.0L, 2.0L};
    divide_point_2D (a, 0.0L);
}

Test(aux_functions, divide_point_2D_success)
{
    point_2D a = {10.0L, 20.0L};
    point_2D b = divide_point_2D(a, 2.0L);
    point_2D c = divide_point_2D(a, 5.0L);
    point_2D d = divide_point_2D(a, 10.0L);

    cr_assert_float_eq(b.x, 5.0L, 1e-15L, "divide_point_2D failed: expected 5.0, but got %Lf", b.x);
    cr_assert_float_eq(b.y, 10.0L, 1e-15L, "divide_point_2D failed: expected 10.0, but got %Lf", b.y);
    cr_assert_float_eq(c.x, 2.0L, 1e-15L, "divide_point_2D failed: expected 2.0, but got %Lf", c.x);
    cr_assert_float_eq(c.y, 4.0L, 1e-15L, "divide_point_2D failed: expected 4.0, but got %Lf", c.y);
    cr_assert_float_eq(d.x, 1.0L, 1e-15L, "divide_point_2D failed: expected 1.0, but got %Lf", d.x);
    cr_assert_float_eq(d.y, 2.0L, 1e-15L, "divide_point_2D failed: expected 2.0, but got %Lf", d.y);
}

/* Test for identity_ders_ell */
Test(aux_functions, identity_ders_ell)
{
    ell_jacobian_2D a = identity_ders_ell ();

    cr_assert_eq(a.ddx.x, 1.0L, "identity_ders_ell failed: ddx.x expected 1.0");
    cr_assert_eq(a.ddx.y, 0.0L, "identity_ders_ell failed: ddx.y expected 0.0");
    cr_assert_eq(a.ddy.x, 0.0L, "identity_ders_ell failed: ddy.x expected 0.0");
    cr_assert_eq(a.ddy.y, 1.0L, "identity_ders_ell failed: ddy.y expected 1.0");
}

/* Test for equals_2D_point_2D_array */ 
Test(aux_functions, equals_2D_point_2D_array)
{
    point_2D **a_test, **b;
    int      i, j;

    // Create 1st array
    a_test = allocate_2D_point_2D_array ("a_test", 2, 2);
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            a_test[i][j] = zero (point_2D);
        }
    }
    a_test[1][1].x = 1.0L;
    a_test[1][1].y = 1.0L;

    // Call function
    b = allocate_2D_point_2D_array ("b", 2, 2);
    equals_2D_point_2D_array(2, 2, a_test, &b);

    // Test output array
    cr_assert_float_eq(b[0][0].x, 0.0L, 1e-15L, "b[0][0].x failed: expected 0.0, but got %Lf", b[0][0].x);
    cr_assert_float_eq(b[0][0].y, 0.0L, 1e-15L, "b[0][0].y failed: expected 0.0, but got %Lf", b[0][0].y);
    cr_assert_float_eq(b[0][1].x, 0.0L, 1e-15L, "b[0][1].x failed: expected 0.0, but got %Lf", b[0][1].x);
    cr_assert_float_eq(b[0][1].y, 0.0L, 1e-15L, "b[0][1].y failed: expected 0.0, but got %Lf", b[0][1].y);
    cr_assert_float_eq(b[1][0].x, 0.0L, 1e-15L, "b[1][0].x failed: expected 0.0, but got %Lf", b[1][0].x);
    cr_assert_float_eq(b[1][0].y, 0.0L, 1e-15L, "b[1][0].y failed: expected 0.0, but got %Lf", b[1][0].y);
    cr_assert_float_eq(b[1][1].x, 1.0L, 1e-15L, "b[1][1].x failed: expected 0.0, but got %Lf", b[1][1].x);
    cr_assert_float_eq(b[1][1].y, 1.0L, 1e-15L, "b[1][1].y failed: expected 0.0, but got %Lf", b[1][1].y);

    // Free memory
    free_2D_point_2D_array("a_test", 2, a_test);
    free_2D_point_2D_array("b", 2, b);
}

/* Test for array_max_2D_long_double */
Test(aux_functions, array_max_2D_long_double)
{
    long double         **a_test, a_max; 
    int                 i, j; 

    // Create array
    a_test = allocate_2D_long_double_array ("a_test", 3, 3);
    a_test[0][0] = 31.0L;
    a_test[0][1] = 95.0L;
    a_test[0][2] = 91.0L;
    a_test[1][0] = 18.0L;
    a_test[1][1] = 63.0L;
    a_test[1][2] = 66.0L;
    a_test[2][0] = 70.0L;
    a_test[2][1] = 23.0L;
    a_test[2][2] = 63.0L;

    // Call function
    a_max = array_max_2D_long_double(3, 3, a_test);

    // Test value
    cr_assert_float_eq (a_max, 95.0L, 1e-12L, "2D array max failed: expected 95.0, but got %Lf", a_max);

    // Free memory
    free_2D_long_double_array ("a_test", 3, a_test);
}

/* Test for convert_ell_jacobian_2D */
Test(aux_functions, convert_ell_jacobian_2D)
{
    ell_jacobian_2D **mat;
    long double     **array;

    // Allocate input matrix (2x2) and output array (4x4)
    mat   = allocate_2D_ell_jacobian_2D_array ("mat", 2, 2);
    array = allocate_2D_long_double_array ("array", 4, 4);

    // Fill input matrix
    mat[0][0].ddx.x = 1.0L;  mat[0][0].ddx.y = 2.0L;
    mat[0][0].ddy.x = 3.0L;  mat[0][0].ddy.y = 4.0L;

    mat[0][1].ddx.x = 5.0L;  mat[0][1].ddx.y = 6.0L;
    mat[0][1].ddy.x = 7.0L;  mat[0][1].ddy.y = 8.0L;

    mat[1][0].ddx.x = 9.0L;  mat[1][0].ddx.y = 10.0L;
    mat[1][0].ddy.x = 11.0L; mat[1][0].ddy.y = 12.0L;

    mat[1][1].ddx.x = 13.0L; mat[1][1].ddx.y = 14.0L;
    mat[1][1].ddy.x = 15.0L; mat[1][1].ddy.y = 16.0L;

    // Call function
    convert_ell_jacobian_2D (2, mat, &array);

    // Test output array
    cr_assert_float_eq (array[0][0],  1.0L, 1e-15L, "array[0][0]: expected 1.0, got %Lf",  array[0][0]);
    cr_assert_float_eq (array[0][1],  3.0L, 1e-15L, "array[0][1]: expected 3.0, got %Lf",  array[0][1]);
    cr_assert_float_eq (array[0][2],  5.0L, 1e-15L, "array[0][2]: expected 5.0, got %Lf",  array[0][2]);
    cr_assert_float_eq (array[0][3],  7.0L, 1e-15L, "array[0][3]: expected 7.0, got %Lf",  array[0][3]);

    cr_assert_float_eq (array[1][0],  2.0L, 1e-15L, "array[1][0]: expected 2.0, got %Lf",  array[1][0]);
    cr_assert_float_eq (array[1][1],  4.0L, 1e-15L, "array[1][1]: expected 4.0, got %Lf",  array[1][1]);
    cr_assert_float_eq (array[1][2],  6.0L, 1e-15L, "array[1][2]: expected 6.0, got %Lf",  array[1][2]);
    cr_assert_float_eq (array[1][3],  8.0L, 1e-15L, "array[1][3]: expected 8.0, got %Lf",  array[1][3]);

    cr_assert_float_eq (array[2][0],  9.0L, 1e-15L, "array[2][0]: expected 9.0, got %Lf",  array[2][0]);
    cr_assert_float_eq (array[2][1], 11.0L, 1e-15L, "array[2][1]: expected 11.0, got %Lf", array[2][1]);
    cr_assert_float_eq (array[2][2], 13.0L, 1e-15L, "array[2][2]: expected 13.0, got %Lf", array[2][2]);
    cr_assert_float_eq (array[2][3], 15.0L, 1e-15L, "array[2][3]: expected 15.0, got %Lf", array[2][3]);

    cr_assert_float_eq (array[3][0], 10.0L, 1e-15L, "array[3][0]: expected 10.0, got %Lf", array[3][0]);
    cr_assert_float_eq (array[3][1], 12.0L, 1e-15L, "array[3][1]: expected 12.0, got %Lf", array[3][1]);
    cr_assert_float_eq (array[3][2], 14.0L, 1e-15L, "array[3][2]: expected 14.0, got %Lf", array[3][2]);
    cr_assert_float_eq (array[3][3], 16.0L, 1e-15L, "array[3][3]: expected 16.0, got %Lf", array[3][3]);

    // Free memory
    free_2D_ell_jacobian_2D_array ("mat", 2, mat);
    free_2D_long_double_array ("array", 4, array);
}

/* Test for hash: verify 2D-to-1D array index mapping */
Test(aux_functions, hash)
{
    cr_assert_eq(hash(1, 1, 5), 0, "hash(1, 1, 5): expected 0, got %d", hash(1, 1, 5));
}

/* Test for convert_point_2D and convert_to_point_2D */ 
Test(aux_functions, convert_point_2D)
{
    point_2D                *vec_in, *vec_out; 
    long double             *array; 

    // Allocate inputs and outputs
    vec_in = allocate_1D_point_2D_array ("vec_in", 3);
    vec_out = allocate_1D_point_2D_array ("vec_out", 3);
    array = allocate_1D_long_double_array ("array", 6);

    // Fill input vector
    vec_in[0].x = 1.0L; vec_in[0].y = 2.0L;
    vec_in[1].x = 3.0L; vec_in[1].y = 4.0L;
    vec_in[2].x = 5.0L; vec_in[2].y = 6.0L;

    // Convert point_2D to long double
    convert_point_2D (3, vec_in, &array);

    // Test 1D long double array
    cr_assert_float_eq (array[0], 1.0L, 1e-15L, "long double array[0]: expected 1.0, got %Lf", array[0]);
    cr_assert_float_eq (array[1], 2.0L, 1e-15L, "long double array[1]: expected 2.0, got %Lf", array[1]);
    cr_assert_float_eq (array[2], 3.0L, 1e-15L, "long double array[2]: expected 3.0, got %Lf", array[2]);
    cr_assert_float_eq (array[3], 4.0L, 1e-15L, "long double array[3]: expected 4.0, got %Lf", array[3]);
    cr_assert_float_eq (array[4], 5.0L, 1e-15L, "long double array[4]: expected 5.0, got %Lf", array[4]);
    cr_assert_float_eq (array[5], 6.0L, 1e-15L, "long double array[5]: expected 6.0, got %Lf", array[5]);

    // Convert long double to point_2D
    convert_to_point_2D (3, array, &vec_out);

    // Test point_2D array
    cr_assert_float_eq (vec_out[0].x, 1.0L, 1e-15L, "point_2D vec_out[0].x: expected 1.0, got %Lf", vec_out[0].x);
    cr_assert_float_eq (vec_out[0].y, 2.0L, 1e-15L, "point_2D vec_out[0].y: expected 2.0, got %Lf", vec_out[0].y);
    cr_assert_float_eq (vec_out[1].x, 3.0L, 1e-15L, "point_2D vec_out[1].x: expected 3.0, got %Lf", vec_out[1].x);
    cr_assert_float_eq (vec_out[1].y, 4.0L, 1e-15L, "point_2D vec_out[1].y: expected 4.0, got %Lf", vec_out[1].y);
    cr_assert_float_eq (vec_out[2].x, 5.0L, 1e-15L, "point_2D vec_out[2].x: expected 5.0, got %Lf", vec_out[2].x);
    cr_assert_float_eq (vec_out[2].y, 6.0L, 1e-15L, "point_2D vec_out[2].y: expected 6.0, got %Lf", vec_out[2].y);

    // Free memory
    free_1D_point_2D_array ("vec_in", vec_in);
    free_1D_point_2D_array ("vec_out", vec_out);
    free_1D_long_double_array ("array", array);
}

/* Test for convert_bh_jacobian_2D */
Test(aux_functions, convert_bh_jacobian_2D)
{
    bh_jacobian_2D  **mat;
    long double     **array;

    // Allocate input matrix (1x1) and output array (4x4)
    mat   = allocate_2D_bh_jacobian_2D_array ("mat", 1, 1);
    array = allocate_2D_long_double_array ("array", 4, 4);

    // Fill input matrix
    mat[0][0].ddx.x = 1.0L;  mat[0][0].ddx.y = 2.0L;
    mat[0][0].ddx.p = 3.0L;  mat[0][0].ddx.q = 4.0L;

    mat[0][0].ddy.x = 5.0L;  mat[0][0].ddy.y = 6.0L;
    mat[0][0].ddy.p = 7.0L;  mat[0][0].ddy.q = 8.0L;

    mat[0][0].ddp.x = 9.0L;  mat[0][0].ddp.y = 10.0L;
    mat[0][0].ddp.p = 11.0L; mat[0][0].ddp.q = 12.0L;

    mat[0][0].ddq.x = 13.0L; mat[0][0].ddq.y = 14.0L;
    mat[0][0].ddq.p = 15.0L; mat[0][0].ddq.q = 16.0L;

    // Call function
    convert_bh_jacobian_2D (1, mat, &array);

    // Test output array
    cr_assert_float_eq (array[0][0],  1.0L, 1e-15L, "array[0][0]: expected 1.0, got %Lf",  array[0][0]);
    cr_assert_float_eq (array[0][1],  5.0L, 1e-15L, "array[0][1]: expected 5.0, got %Lf",  array[0][1]);
    cr_assert_float_eq (array[0][2],  9.0L, 1e-15L, "array[0][2]: expected 9.0, got %Lf",  array[0][2]);
    cr_assert_float_eq (array[0][3], 13.0L, 1e-15L, "array[0][3]: expected 13.0, got %Lf", array[0][3]);

    cr_assert_float_eq (array[1][0],  2.0L, 1e-15L, "array[1][0]: expected 2.0, got %Lf",  array[1][0]);
    cr_assert_float_eq (array[1][1],  6.0L, 1e-15L, "array[1][1]: expected 6.0, got %Lf",  array[1][1]);
    cr_assert_float_eq (array[1][2], 10.0L, 1e-15L, "array[1][2]: expected 10.0, got %Lf", array[1][2]);
    cr_assert_float_eq (array[1][3], 14.0L, 1e-15L, "array[1][3]: expected 14.0, got %Lf", array[1][3]);

    cr_assert_float_eq (array[2][0],  3.0L, 1e-15L, "array[2][0]: expected 3.0, got %Lf",  array[2][0]);
    cr_assert_float_eq (array[2][1],  7.0L, 1e-15L, "array[2][1]: expected 7.0, got %Lf",  array[2][1]);
    cr_assert_float_eq (array[2][2], 11.0L, 1e-15L, "array[2][2]: expected 11.0, got %Lf", array[2][2]);
    cr_assert_float_eq (array[2][3], 15.0L, 1e-15L, "array[2][3]: expected 15.0, got %Lf", array[2][3]);

    cr_assert_float_eq (array[3][0],  4.0L, 1e-15L, "array[3][0]: expected 4.0, got %Lf",  array[3][0]);
    cr_assert_float_eq (array[3][1],  8.0L, 1e-15L, "array[3][1]: expected 8.0, got %Lf",  array[3][1]);
    cr_assert_float_eq (array[3][2], 12.0L, 1e-15L, "array[3][2]: expected 12.0, got %Lf", array[3][2]);
    cr_assert_float_eq (array[3][3], 16.0L, 1e-15L, "array[3][3]: expected 16.0, got %Lf", array[3][3]);

    // Free memory
    free_2D_bh_jacobian_2D_array ("mat", 1, mat);
    free_2D_long_double_array ("array", 4, array);
}

/* Test for convert_bh_sol_2D and convert_to_bh_sol_2D (round-trip) */
Test(aux_functions, convert_bh_sol_2D)
{
    bh_sol_2D   *vec_in, *vec_out;
    long double *array;

    // Allocate inputs and outputs
    vec_in  = allocate_1D_bh_sol_2D_array ("vec_in", 2);
    vec_out = allocate_1D_bh_sol_2D_array ("vec_out", 2);
    array   = allocate_1D_long_double_array ("array", 8);

    // Fill input vector
    vec_in[0].x = 1.0L; vec_in[0].y = 2.0L;
    vec_in[0].p = 3.0L; vec_in[0].q = 4.0L;
    vec_in[1].x = 5.0L; vec_in[1].y = 6.0L;
    vec_in[1].p = 7.0L; vec_in[1].q = 8.0L;

    // Convert bh_sol_2D -> long double
    convert_bh_sol_2D (2, vec_in, &array);

    // Test intermediate long double array
    cr_assert_float_eq (array[0], 1.0L, 1e-15L, "array[0]: expected 1.0, got %Lf", array[0]);
    cr_assert_float_eq (array[1], 2.0L, 1e-15L, "array[1]: expected 2.0, got %Lf", array[1]);
    cr_assert_float_eq (array[2], 3.0L, 1e-15L, "array[2]: expected 3.0, got %Lf", array[2]);
    cr_assert_float_eq (array[3], 4.0L, 1e-15L, "array[3]: expected 4.0, got %Lf", array[3]);
    cr_assert_float_eq (array[4], 5.0L, 1e-15L, "array[4]: expected 5.0, got %Lf", array[4]);
    cr_assert_float_eq (array[5], 6.0L, 1e-15L, "array[5]: expected 6.0, got %Lf", array[5]);
    cr_assert_float_eq (array[6], 7.0L, 1e-15L, "array[6]: expected 7.0, got %Lf", array[6]);
    cr_assert_float_eq (array[7], 8.0L, 1e-15L, "array[7]: expected 8.0, got %Lf", array[7]);

    // Convert back: long double -> bh_sol_2D
    convert_to_bh_sol_2D (2, array, &vec_out);

    // Test round-trip output
    cr_assert_float_eq (vec_out[0].x, 1.0L, 1e-15L, "vec_out[0].x: expected 1.0, got %Lf", vec_out[0].x);
    cr_assert_float_eq (vec_out[0].y, 2.0L, 1e-15L, "vec_out[0].y: expected 2.0, got %Lf", vec_out[0].y);
    cr_assert_float_eq (vec_out[0].p, 3.0L, 1e-15L, "vec_out[0].p: expected 3.0, got %Lf", vec_out[0].p);
    cr_assert_float_eq (vec_out[0].q, 4.0L, 1e-15L, "vec_out[0].q: expected 4.0, got %Lf", vec_out[0].q);
    cr_assert_float_eq (vec_out[1].x, 5.0L, 1e-15L, "vec_out[1].x: expected 5.0, got %Lf", vec_out[1].x);
    cr_assert_float_eq (vec_out[1].y, 6.0L, 1e-15L, "vec_out[1].y: expected 6.0, got %Lf", vec_out[1].y);
    cr_assert_float_eq (vec_out[1].p, 7.0L, 1e-15L, "vec_out[1].p: expected 7.0, got %Lf", vec_out[1].p);
    cr_assert_float_eq (vec_out[1].q, 8.0L, 1e-15L, "vec_out[1].q: expected 8.0, got %Lf", vec_out[1].q);

    // Free memory
    free_1D_bh_sol_2D_array ("vec_in", vec_in);
    free_1D_bh_sol_2D_array ("vec_out", vec_out);
    free_1D_long_double_array ("array", array);
}

/* Test LU_linear_system_solve */
Test(aux_functions, LU_solve)
{
    int                 n = 3;
    long double         **M, *b, *x; 

    // Allocate arrays
    M = allocate_2D_long_double_array ("LHS matrix", 3, 3);
    b = allocate_1D_long_double_array ("RHS vector", 3);
    x = allocate_1D_long_double_array ("Solution vector", 3);

    // LHS matrix
    M[0][0] = 2.0L; M[0][1] = 1.0L; M[0][2] = -1.0L;
    M[1][0] = -3.0L; M[1][1] = -1.0L; M[1][2] = 2.0L;
    M[2][0] = -2.0L; M[2][1] = 1.0L; M[2][2] = 2.0L;

    // RHS vector
    b[0] = 8.0L; b[1] = -11.0L; b[2] = -3.0L; 

    // Call function
    LU_linear_system_solve(3, M, b, &x);

    // Check solution
    cr_assert_float_eq (x[0], 2.0L, 1e-12L, "LU solve test failed at x[0]: expected 2.0, but got %Lf", x[0]);
    cr_assert_float_eq (x[1], 3.0L, 1e-12L, "LU solve test failed at x[1]: expected 3.0, but got %Lf", x[1]);
    cr_assert_float_eq (x[2], -1.0L, 1e-12L, "LU solve test failed at x[2]: expected -1.0, but got %Lf", x[2]);

    // Free memory
    free_2D_long_double_array ("LHS matrix", 3, M);
    free_1D_long_double_array ("RHS vector", b);
    free_1D_long_double_array ("Solution vector", x);
}

/* Test for gmres_linear_system_solve */ 
Test(aux_functions, gmres_linear_system_solve)
{
    ell_jacobian_2D         **dFdu; 
    point_2D                *fu, *del_u; 

    // Allocate inputs and output 
    dFdu = allocate_2D_ell_jacobian_2D_array ("dFdu", 1, 1);
    fu = allocate_1D_point_2D_array ("fu", 1);
    del_u = allocate_1D_point_2D_array ("del_u", 1);

    // Fill Jacobian matrix (diagonal: M = [[2, 0], [0, 3]])
    dFdu[0][0].ddx.x = 2.0L; dFdu[0][0].ddx.y = 0.0L;
    dFdu[0][0].ddy.x = 0.0L; dFdu[0][0].ddy.y = 3.0L;

    // Fill RHS vector (b = [4, 9])
    fu[0].x = 4.0L; 
    fu[0].y = 9.0L;

    // Call function 
    gmres_linear_system_solve (0, 1, 1, dFdu, fu, &del_u);

    // Test solution 
    cr_assert_float_eq (del_u[0].x, 2.0L, 1e-5L, "GMRES system solve failed at del_u[0].x: expected 2.0, got %Lf", del_u[0].x);
    cr_assert_float_eq (del_u[0].y, 3.0L, 1e-5L, "GMRES system solve failed at del_u[0].y: expected 3.0, got %Lf", del_u[0].y);

    // Free memory 
    free_2D_ell_jacobian_2D_array ("dFdu", 1, dFdu);
    free_1D_point_2D_array ("fu", fu);
    free_1D_point_2D_array ("del_u", del_u);
}
