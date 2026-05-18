#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "config.h"
#include "errors.h"



/*
    Strip leading and trailing whitespaces from 
    a string 
*/
static void trim (char *buf)
{
    int                     i, j, len;

    len = strlen (buf);

    // Find first non-whitespace
    for (i = 0; i < len && isspace ((unsigned char) buf[i]); i++);

    // Find last non-whitespace 
    for (j = len; j > i && isspace ((unsigned char) buf[j - 1]); j--);

    // Shift content left if leading whitespace was present
    if (i > 0)
    {
        memmove (buf, buf + i, j - i);
    }
    buf[j - i] = '\0';
}



/*
    Parse a whitespace separated list of long integers 
    from a string into a dynamically allocated array

    Input parameters: str   - string containg the values
                      count - expected number of values 
                      key   - key name (for error messages)
*/
static long *parse_long_list (char *str, long count, const char *key)
{
    long                    *out;
    char                    *tok, *saveptr;
    int                     i;

    if (count <= 0) 
    {
        fatal_error ("parse_long_list: count must be > 0");
    }
    out = malloc (count * sizeof (long));
    if (out == NULL) 
    {
        fatal_error ("parse_long_list: malloc failed");
    }
            
    i = 0;
    tok = strtok_r (str, " \t", &saveptr);
    while (tok != NULL && i < count)
    {
        out[i++] = strtol (tok, NULL, 10);
        tok = strtok_r (NULL, " \t", &saveptr);
    }
    if (i != count)
    {
        fprintf (stderr, "Error: %s expected %ld values, got %d\n", key, count, i);
        fatal_error ("parse_long_list: wrong value count");
    }

    return out;
}



/*
    Map a config file solver token string to the 
    solver_type enum. Returns SOLVER_UNKNOWN if 
    no match is found

    Input parameters: tok       - solver token string 
*/
static solver_type parse_solver_token (const char *tok)
{
    if (strcmp (tok, "ell") == 0) return SOLVER_ELL;
    if (strcmp (tok, "ell_newton") == 0) return SOLVER_ELL_NEWTON;
    if (strcmp (tok, "psn_tm") == 0) return SOLVER_PSN_TM;
    if (strcmp (tok, "psn_lim") == 0) return SOLVER_PSN_LIM;
    if (strcmp (tok, "newton") == 0) return SOLVER_NEWTON;
    return SOLVER_UNKNOWN;
}



/*
    Parse input file into supplied run_args struct.
    Required fields are validated at file parse end:
    missing ones produce a fatal error

    On entry, *cfg is zero-initialized internally, so any
    npts_x/y/z arrays previously held by cfg will be leaked.
    Caller must pass either a fresh struct or one that has
    been released with free_run_args first.

    Input parameters: filename  - path to the input file
                      cfg       - struct to populate
*/
void read_config_file (char *filename, run_args *cfg)
{
    FILE                    *fptr;
    char                    line[512], key[64], value[256];
    char                    *eq, *hash; 
    char                    npts_x_raw[256], npts_y_raw[256], npts_z_raw[256]; 


    fptr = fopen (filename, "r");
    if (fptr == NULL)
    {
        fatal_error ("read_config_file: could not open config file\n");
    }

    // Initialize all fields (memset zeros cfg->solver to SOLVER_UNKNOWN)
    memset (cfg, 0, sizeof (*cfg));
    npts_x_raw[0] = '\0';
    npts_y_raw[0] = '\0';
    npts_z_raw[0] = '\0';

    // Line-by-line parse 
    while (fgets (line, sizeof(line), fptr) != NULL)
    {
        // Strip comments
        hash = strchr (line, '#');
        if (hash != NULL) *hash = '\0';

        trim (line);
        if (strlen (line) == 0) continue; 

        // Split on `=` 
        eq = strchr (line, '=');
        if (eq == NULL)
        {
            fprintf (stderr, "Warning: malformed line (no `=`), skipping: %s\n", line);
            continue;
        }
        *eq = '\0';
        strncpy (key, line, sizeof(key) - 1);
        strncpy (value, eq + 1, sizeof(value) - 1);
        key[sizeof(key) - 1] = '\0';
        value[sizeof(value) - 1] = '\0';
        trim (key);
        trim (value);

        // Required fields
        if (strcmp (key, "mode") == 0) 
        {
            strncpy (cfg->mode, value, sizeof(cfg->mode) - 1);
        } 
        else if (strcmp (key, "niter") == 0)
        {
            cfg->niter = strtol (value, NULL, 10);
        } 
        else if (strcmp (key, "nblocks") == 0)
        {
            cfg->nblocks = strtol (value, NULL, 10);
        } 
        else if (strcmp (key, "init_grid_file") == 0)
        {
            strncpy (cfg->init_grid_file, value, sizeof(cfg->init_grid_file) - 1);
        } 
        else if (strcmp (key, "npts_x") == 0)
        {
            strncpy (npts_x_raw, value, sizeof(npts_x_raw) - 1);
        } 
        else if (strcmp (key, "npts_y") == 0)
        {
            strncpy (npts_y_raw, value, sizeof(npts_y_raw) - 1);
        } 
        else if (strcmp (key, "npts_z") == 0)
        {
            strncpy (npts_z_raw, value, sizeof(npts_z_raw) - 1);
        } 
        else if (strcmp (key, "solver") == 0)
        {
            cfg->solver = parse_solver_token (value);
            if (cfg->solver == SOLVER_UNKNOWN)
            {
                fprintf (stderr, "WARNING: Unrecognized solver token '%s'\n", value);
            }
        } 

        // Optional fields
        else if (strcmp (key, "lambda_max") == 0)
        {
            cfg->lambda_max = strtold (value, NULL); 
            cfg->lambda_max_set = 1;
        } 
        else if (strcmp (key, "ramp_iters") == 0)
        {
            cfg->ramp_iters = (int) strtol (value, NULL, 10); 
            cfg->ramp_iters_set = 1;
        } 
        else if (strcmp (key, "inner_niter") == 0)
        {
            cfg->inner_niter = (int) strtol (value, NULL, 10); 
            cfg->inner_niter_set = 1;
        } 
        else if (strcmp (key, "omega") == 0)
        {
            cfg->omega = strtold (value, NULL); 
            cfg->omega_set = 1;
        } 
        else if (strcmp (key, "tol_step") == 0)
        {
            cfg->tol_step = strtold (value, NULL); 
            cfg->tol_step_set = 1;
        } 
        else if (strcmp (key, "tol_resid") == 0)
        {
            cfg->tol_resid = strtold (value, NULL); 
            cfg->tol_resid_set = 1;
        } 
        else 
        {
            fprintf (stderr, "WARNING: Unrecognized key '%s', ignoring\n", key);
        }
    }

    // Close file
    fclose (fptr);


    // Validate required fields 
    if (strlen (cfg->mode) == 0) 
    {
        fatal_error ("Required field 'mode' missing");
    }
    if (strcmp (cfg->mode, "2D") != 0 && strcmp (cfg->mode, "3D") != 0)
    {
        fprintf (stderr, "Error: mode = '%s' (must be \"2D\" or \"3D\")\n", cfg->mode);
        fatal_error ("Field 'mode' has an invalid value");
    }
    if (cfg->niter <= 0) 
    {
        fatal_error ("Required field 'niter' missing or invalid");
    }
    if (cfg->solver == SOLVER_UNKNOWN) 
    {
        fatal_error ("Required field 'solver' missing or unrecognized");
    }

    // Initial grid input
    if (strlen (cfg->init_grid_file) > 0)
    {
        // Input read from PLOT3D file 
        cfg->nblocks = 0;
        cfg->npts_x = NULL;
        cfg->npts_y = NULL;
        cfg->npts_z = NULL;
    }
    else 
    {
        if (cfg->nblocks <= 0)
        {
            fatal_error ("nblocks > 0 required when an initial grid file was not provided");
        }
        if (strlen (npts_x_raw) == 0)
        {
            fatal_error ("npts_x required when an initial grid file was not provided");
        }
        if (strlen (npts_y_raw) == 0)
        {
            fatal_error ("npts_y required when an initial grid file was not provided");
        }

        cfg->npts_x = parse_long_list (npts_x_raw, cfg->nblocks, "npts_x");
        cfg->npts_y = parse_long_list (npts_y_raw, cfg->nblocks, "npts_y");
        if (strcmp (cfg->mode, "3D") == 0)
        {
            if (strlen (npts_z_raw) == 0)
            {
                fatal_error ("npts_z required for 3D mode");
            }
            cfg->npts_z = parse_long_list (npts_z_raw, cfg->nblocks, "npts_z");
        }
        else 
        {
            cfg->npts_z = NULL;
        }
    }
}



/*
    Fill in defaults for any optional fields not set in
    the config file. The defaults are:
        SOLVER_NEWTON       lambda_max  = 0.1
                            ramp_iters  = 20
        SOLVER_PSN_LIM      lambda_max  = 0.2
                            ramp_iters  = 500
                            inner_niter = 10
    Universal defaults:     omega       = 1.0
                            tol_step    = 1E-8
                            tol_resid   = 1E-6

    Issues a warning when a default is being applied to a
    field that is applicable to the chosen solver, so the
    user knows they are getting a default

    Input parameters: cfg   - config struct to fill in
*/
void apply_solver_defaults (run_args *cfg)
{
    int                     blend_solver, pointwise_solver, lim_solver;


    // Determine which solver category is in use
    blend_solver = (cfg->solver == SOLVER_NEWTON) || (cfg->solver == SOLVER_PSN_LIM);
    pointwise_solver = (cfg->solver == SOLVER_ELL) || (cfg->solver == SOLVER_PSN_TM);
    lim_solver = (cfg->solver == SOLVER_PSN_LIM);

    // lambda_max: applies to newton and psn_lim
    if (!cfg->lambda_max_set)
    {
        cfg->lambda_max = (cfg->solver == SOLVER_NEWTON) ? 0.1L : 0.2L;
        if (blend_solver)
        {
            fprintf (stderr, "Warning: lambda_max not specified, using default %.4Le\n", cfg->lambda_max);
        }
    }

    // ramp_iters: applies to newton and psn_lim
    if (!cfg->ramp_iters_set)
    {
        cfg->ramp_iters = (cfg->solver == SOLVER_NEWTON) ? 20 : 500;
        if (blend_solver)
        {
            fprintf (stderr, "Warning: ramp_iters not specified, using default %d\n", cfg->ramp_iters);
        }
    }

    // inner_niter: applies only to psn_lim
    if (!cfg->inner_niter_set)
    {
        cfg->inner_niter = 10;
        if (lim_solver)
        {
            fprintf (stderr, "Warning: inner_niter not specified, using default %d\n", cfg->inner_niter);
        }
    }

    // omega: applies to ell and psn_tm
    if (!cfg->omega_set)
    {
        cfg->omega = 1.0L;
        if (pointwise_solver)
        {
            fprintf (stderr, "Warning: omega not specified, using default %.4Le\n", cfg->omega);
        }
    }

    // tol_step and tol_resid apply to all solvers; warn unconditionally
    if (!cfg->tol_step_set)
    {
        cfg->tol_step = 1.0E-8L;
        fprintf (stderr, "Warning: tol_step not specified, using default %.4Le\n", cfg->tol_step);
    }
    if (!cfg->tol_resid_set)
    {
        cfg->tol_resid = 1.0E-6L;
        fprintf (stderr, "Warning: tol_resid not specified, using default %.4Le\n", cfg->tol_resid);
    }
}



/*
    Release dynamically allocated members of cfg (the
    npts_x/y/z arrays if non-NULL)

    Input parameters: cfg   - config struct to clean up
*/
void free_run_args (run_args *cfg)
{
    if (cfg->npts_x != NULL)
    {
        free (cfg->npts_x);
        cfg->npts_x = NULL;
    }
    if (cfg->npts_y != NULL)
    {
        free (cfg->npts_y);
        cfg->npts_y = NULL;
    }
    if (cfg->npts_z != NULL)
    {
        free (cfg->npts_z);
        cfg->npts_z = NULL;
    }
}
