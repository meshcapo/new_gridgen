#include <criterion/criterion.h>
#include <criterion/parameterized.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <linux/limits.h>
#include "config.h"



/* Build the path tests/files/config/<basename> relative to this
   source file, so tests can be run from any working directory.
   Returns a pointer to a static buffer; not thread-safe but
   adequate for sequential Criterion test runs */
static const char *config_test_file (const char *basename)
{
    static char             path[PATH_MAX];
    char                    src_path[PATH_MAX];
    char                    *test_dir;

    strncpy (src_path, __FILE__, PATH_MAX - 1);
    src_path[PATH_MAX - 1] = '\0';
    test_dir = dirname (src_path);
    snprintf (path, sizeof(path), "%s/files/config/%s", test_dir, basename);
    return path;
}



/* Write `content` to a fresh temp file under /tmp and return its
   path in a static buffer. Used by the fatal-path tests so each
   bad-input case can carry its own config text inline rather
   than living as a separate file under tests/files/config/.
   Caller does not need to unlink — Criterion forks each test, so
   the file leaks at most once per parameterized run */
static const char *write_temp_config (const char *content)
{
    static char             path[PATH_MAX];
    FILE                    *fp;

    snprintf (path, sizeof(path), "/tmp/gridgen_test_config_%d.in", getpid());
    fp = fopen (path, "w");
    if (fp == NULL) return NULL;
    fputs (content, fp);
    fclose (fp);
    return path;
}



/* Test: required-only 2D config with explicit npts parses correctly.
   Every required field should land in cfg, *_set flags for optionals
   stay 0 (since the file does not set any optionals), and npts_z is
   NULL since mode is 2D */
Test (config, valid_2d_explicit)
{
    run_args                cfg;
    char                    path[PATH_MAX];

    // Resolve fixture path
    strncpy (path, config_test_file ("valid_2d_explicit.in"), PATH_MAX - 1);
    path[PATH_MAX - 1] = '\0';

    // Parse
    read_config_file (path, &cfg);

    // Required fields populated
    cr_assert_str_eq (cfg.mode, "2D");
    cr_assert_eq (cfg.niter, 500);
    cr_assert_eq (cfg.nblocks, 1);
    cr_assert_not_null (cfg.npts_x);
    cr_assert_not_null (cfg.npts_y);
    cr_assert_eq (cfg.npts_x[0], 101);
    cr_assert_eq (cfg.npts_y[0], 31);
    cr_assert_null (cfg.npts_z);
    cr_assert_eq (cfg.solver, SOLVER_PSN_LIM);

    // No optionals set -> *_set flags all zero
    cr_assert_eq (cfg.lambda_max_set, 0);
    cr_assert_eq (cfg.ramp_iters_set, 0);
    cr_assert_eq (cfg.inner_niter_set, 0);
    cr_assert_eq (cfg.omega_set, 0);
    cr_assert_eq (cfg.tol_step_set, 0);
    cr_assert_eq (cfg.tol_resid_set, 0);

    // Cleanup
    free_run_args (&cfg);
}



/* Test: 2D config with init_grid_file set (no explicit npts).
   Parser should leave nblocks=0 and npts_x/y/z = NULL — main fills
   these from the PLOT3D header. init_grid_file string should be
   stored verbatim */
Test (config, valid_2d_init_grid)
{
    run_args                cfg;

    read_config_file ((char *) config_test_file ("valid_2d_init_grid.in"), &cfg);

    cr_assert_str_eq (cfg.mode, "2D");
    cr_assert_eq (cfg.niter, 100);
    cr_assert_str_eq (cfg.init_grid_file, "some_grid.x");
    cr_assert_eq (cfg.solver, SOLVER_ELL_NEWTON);

    // npts arrays left NULL when init_grid_file is provided
    cr_assert_eq (cfg.nblocks, 0);
    cr_assert_null (cfg.npts_x);
    cr_assert_null (cfg.npts_y);
    cr_assert_null (cfg.npts_z);

    free_run_args (&cfg);
}



/* Test: 3D config with explicit npts parses npts_z correctly */
Test (config, valid_3d)
{
    run_args                cfg;

    read_config_file ((char *) config_test_file ("valid_3d.in"), &cfg);

    cr_assert_str_eq (cfg.mode, "3D");
    cr_assert_eq (cfg.niter, 200);
    cr_assert_eq (cfg.nblocks, 1);
    cr_assert_not_null (cfg.npts_x);
    cr_assert_not_null (cfg.npts_y);
    cr_assert_not_null (cfg.npts_z);
    cr_assert_eq (cfg.npts_x[0], 41);
    cr_assert_eq (cfg.npts_y[0], 21);
    cr_assert_eq (cfg.npts_z[0], 11);
    cr_assert_eq (cfg.solver, SOLVER_ELL);

    free_run_args (&cfg);
}



/* Test: every optional field set in the config file is read and
   the corresponding *_set flag is raised. Values should be stored
   verbatim with no defaults applied (defaults are a separate
   function tested elsewhere) */
Test (config, valid_with_all_optionals)
{
    run_args                cfg;

    read_config_file ((char *) config_test_file ("valid_with_all_optionals.in"), &cfg);

    // Solver and required fields
    cr_assert_eq (cfg.solver, SOLVER_NEWTON);
    cr_assert_eq (cfg.niter, 1000);

    // Each optional field landed
    cr_assert_eq (cfg.lambda_max_set, 1);
    cr_assert_eq (cfg.ramp_iters_set, 1);
    cr_assert_eq (cfg.inner_niter_set, 1);
    cr_assert_eq (cfg.omega_set, 1);
    cr_assert_eq (cfg.tol_step_set, 1);
    cr_assert_eq (cfg.tol_resid_set, 1);

    // Optional values match the file
    cr_assert_float_eq (cfg.lambda_max, 0.07L, 1e-15L);
    cr_assert_eq (cfg.ramp_iters, 33);
    cr_assert_eq (cfg.inner_niter, 7);
    cr_assert_float_eq (cfg.omega, 0.85L, 1e-15L);
    cr_assert_float_eq (cfg.tol_step, 1.5e-9L, 1e-20L);
    cr_assert_float_eq (cfg.tol_resid, 2.5e-7L, 1e-18L);

    free_run_args (&cfg);
}



/* Test: comments, blank lines, and erratic whitespace around keys,
   values, and `=` are all handled correctly. Confirms that the
   trim() helper works for keys with leading whitespace (which would
   trigger the buf[j-i] null-terminator path) and that inline
   `# comments` are stripped */
Test (config, comments_and_whitespace)
{
    run_args                cfg;

    read_config_file ((char *) config_test_file ("comments_and_whitespace.in"), &cfg);

    cr_assert_str_eq (cfg.mode, "2D");
    cr_assert_eq (cfg.niter, 250);
    cr_assert_eq (cfg.nblocks, 1);
    cr_assert_not_null (cfg.npts_x);
    cr_assert_not_null (cfg.npts_y);
    cr_assert_eq (cfg.npts_x[0], 61);
    cr_assert_eq (cfg.npts_y[0], 17);
    cr_assert_eq (cfg.solver, SOLVER_PSN_TM);
    cr_assert_eq (cfg.omega_set, 1);
    cr_assert_float_eq (cfg.omega, 0.9L, 1e-15L);

    free_run_args (&cfg);
}



/* Each fatal-path test below loads a config that triggers
   fatal_error() inside read_config_file. fatal_error calls exit(),
   so we use Criterion's .exit_code = EXIT_FAILURE assertion: the
   forked test process should terminate with that code, no body
   assertions needed */

/* Test: nonexistent file -> fopen failure -> fatal_error */
Test (config, fatal_nonexistent_file, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;

    read_config_file ("/tmp/this_path_definitely_does_not_exist_gridgen.in", &cfg);
}


/* Test: missing required field `mode` */
Test (config, fatal_missing_mode, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "niter = 100\n"
        "nblocks = 1\n"
        "npts_x = 11\n"
        "npts_y = 11\n"
        "solver = ell\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: mode set to a value that is neither "2D" nor "3D" */
Test (config, fatal_bad_mode, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = foo\n"
        "niter = 100\n"
        "nblocks = 1\n"
        "npts_x = 11\n"
        "npts_y = 11\n"
        "solver = ell\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: missing required field `niter` (memset gives niter=0,
   validation rejects niter <= 0) */
Test (config, fatal_missing_niter, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = 2D\n"
        "nblocks = 1\n"
        "npts_x = 11\n"
        "npts_y = 11\n"
        "solver = ell\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: missing required field `solver` */
Test (config, fatal_missing_solver, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = 2D\n"
        "niter = 100\n"
        "nblocks = 1\n"
        "npts_x = 11\n"
        "npts_y = 11\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: solver token does not match any known solver. Parser
   warns to stderr then leaves solver = SOLVER_UNKNOWN; validation
   rejects */
Test (config, fatal_bad_solver_token, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = 2D\n"
        "niter = 100\n"
        "nblocks = 1\n"
        "npts_x = 11\n"
        "npts_y = 11\n"
        "solver = not_a_real_solver\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: nblocks omitted (memset gives 0) when no init_grid_file */
Test (config, fatal_nblocks_zero, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = 2D\n"
        "niter = 100\n"
        "npts_x = 11\n"
        "npts_y = 11\n"
        "solver = ell\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: npts_x omitted when no init_grid_file */
Test (config, fatal_missing_npts_x, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = 2D\n"
        "niter = 100\n"
        "nblocks = 1\n"
        "npts_y = 11\n"
        "solver = ell\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: npts_y omitted when no init_grid_file */
Test (config, fatal_missing_npts_y, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = 2D\n"
        "niter = 100\n"
        "nblocks = 1\n"
        "npts_x = 11\n"
        "solver = ell\n");
    read_config_file ((char *) path, &cfg);
}


/* Test: 3D mode requires npts_z when no init_grid_file */
Test (config, fatal_missing_npts_z_3d, .exit_code = EXIT_FAILURE)
{
    run_args                cfg;
    const char              *path;

    path = write_temp_config (
        "mode = 3D\n"
        "niter = 100\n"
        "nblocks = 1\n"
        "npts_x = 11\n"
        "npts_y = 11\n"
        "solver = ell\n");
    read_config_file ((char *) path, &cfg);
}



/* The apply_solver_defaults tests below build cfg directly (no
   parser involvement) so the function is tested in isolation. We
   verify two things: (1) when *_set flags are 0, the function
   writes the documented per-solver defaults; (2) when *_set flags
   are 1, the function leaves the existing values untouched */

/* Test: SOLVER_NEWTON gets lambda_max = 0.1 and ramp_iters = 20.
   Universal defaults (omega, inner_niter, tol_step, tol_resid)
   are also written even though they may not all apply */
Test (config, defaults_newton)
{
    run_args                cfg;

    // Build a cfg with solver = NEWTON and all *_set = 0
    memset (&cfg, 0, sizeof(cfg));
    cfg.solver = SOLVER_NEWTON;

    apply_solver_defaults (&cfg);

    // Solver-specific defaults
    cr_assert_float_eq (cfg.lambda_max, 0.1L, 1e-15L);
    cr_assert_eq (cfg.ramp_iters, 20);

    // Universal defaults
    cr_assert_eq (cfg.inner_niter, 10);
    cr_assert_float_eq (cfg.omega, 1.0L, 1e-15L);
    cr_assert_float_eq (cfg.tol_step, 1.0e-8L, 1e-20L);
    cr_assert_float_eq (cfg.tol_resid, 1.0e-6L, 1e-18L);
}


/* Test: SOLVER_PSN_LIM gets the LIM-flavoured defaults
   (lambda_max = 0.2, ramp_iters = 500) */
Test (config, defaults_psn_lim)
{
    run_args                cfg;

    memset (&cfg, 0, sizeof(cfg));
    cfg.solver = SOLVER_PSN_LIM;

    apply_solver_defaults (&cfg);

    cr_assert_float_eq (cfg.lambda_max, 0.2L, 1e-15L);
    cr_assert_eq (cfg.ramp_iters, 500);
    cr_assert_eq (cfg.inner_niter, 10);
    cr_assert_float_eq (cfg.omega, 1.0L, 1e-15L);
    cr_assert_float_eq (cfg.tol_step, 1.0e-8L, 1e-20L);
    cr_assert_float_eq (cfg.tol_resid, 1.0e-6L, 1e-18L);
}


/* Test: SOLVER_ELL falls into the non-NEWTON branch for
   lambda_max (gets 0.2) and the non-NEWTON branch for ramp_iters
   (gets 500). These values are inert for ELL but the function
   writes them unconditionally */
Test (config, defaults_ell)
{
    run_args                cfg;

    memset (&cfg, 0, sizeof(cfg));
    cfg.solver = SOLVER_ELL;

    apply_solver_defaults (&cfg);

    cr_assert_float_eq (cfg.lambda_max, 0.2L, 1e-15L);
    cr_assert_eq (cfg.ramp_iters, 500);
    cr_assert_float_eq (cfg.omega, 1.0L, 1e-15L);
    cr_assert_eq (cfg.inner_niter, 10);
    cr_assert_float_eq (cfg.tol_step, 1.0e-8L, 1e-20L);
    cr_assert_float_eq (cfg.tol_resid, 1.0e-6L, 1e-18L);
}


/* Test: pre-set values (with their *_set flags raised) are
   preserved by apply_solver_defaults — the function only writes
   when a *_set flag is 0. This is what allows a user to override
   any subset of tunables in the config file */
Test (config, defaults_set_flags_prevent_overwrite)
{
    run_args                cfg;

    // Build a cfg with solver = NEWTON and every optional pre-set
    // to a non-default sentinel value
    memset (&cfg, 0, sizeof(cfg));
    cfg.solver = SOLVER_NEWTON;
    cfg.lambda_max = 0.42L;
    cfg.lambda_max_set = 1;
    cfg.ramp_iters = 99;
    cfg.ramp_iters_set = 1;
    cfg.inner_niter = 13;
    cfg.inner_niter_set = 1;
    cfg.omega = 0.5L;
    cfg.omega_set = 1;
    cfg.tol_step = 3.0e-10L;
    cfg.tol_step_set = 1;
    cfg.tol_resid = 5.0e-9L;
    cfg.tol_resid_set = 1;

    apply_solver_defaults (&cfg);

    // None of the values should have been overwritten
    cr_assert_float_eq (cfg.lambda_max, 0.42L, 1e-15L);
    cr_assert_eq (cfg.ramp_iters, 99);
    cr_assert_eq (cfg.inner_niter, 13);
    cr_assert_float_eq (cfg.omega, 0.5L, 1e-15L);
    cr_assert_float_eq (cfg.tol_step, 3.0e-10L, 1e-20L);
    cr_assert_float_eq (cfg.tol_resid, 5.0e-9L, 1e-18L);
}
