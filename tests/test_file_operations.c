#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <criterion/parameterized.h>
#include <stdio.h>
#include <libgen.h>
#include <unistd.h>
#include <linux/limits.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#ifdef HAVE_LIBXML2
#include <libxml/parser.h>
#include <libxml/tree.h>
#endif
#include "file_operations.h"
#include "aux_functions.h"


/* Helper function to create a uniform grid on a 
   square grid [0, 1] x [0, 1] domain */
static point_2D **make_unit_square_grid(int nx, int ny)
{
    point_2D **grid = allocate_2D_point_2D_array("grid", nx, ny);
    if (grid == NULL) return NULL;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            grid[i][j].x = (long double) i/(long double)(nx - 1);
            grid[i][j].y = (long double) j/(long double)(ny - 1);
        }
    }
    return grid;
}


/* Parameterized test for function get_file_extension 
   No failure test cases since this function calls 
   fatal_error which is tested elsewhere */
typedef struct  /* Data type for parameterized test */
{
    char filename[64];
    char expected_extension[16];
} ext_case_t;

ParameterizedTestParameters (file_operations, get_file_extension_param) /* Generate parameter set */
{
    static ext_case_t cases[] =
    {
        {"grid.x", ".x"},
        {"grid.vts", ".vts"},
    };
    return cr_make_param_array (ext_case_t, cases, sizeof(cases)/sizeof(*cases), NULL);
}

ParameterizedTest (ext_case_t *param, file_operations, get_file_extension_param)
{
    char *extension = get_file_extension (param->filename);
    cr_assert_not_null (extension, "Expected non-null extension pointer for filename '%s'", param->filename);
    cr_assert_str_eq (extension, param-> expected_extension, "For input '%s', expected extensiom '%s', but got '%s'", 
                      param->filename, param->expected_extension, extension);   
}


/* Test for plot3D file reading functions */
Test (file_operations, read_plot3D)
{
    // Compute relative path from test source file
    char src_path[PATH_MAX];
    strncpy (src_path, __FILE__, PATH_MAX - 1);
    src_path[PATH_MAX - 1] = '\0';
    char *test_dir = dirname (src_path);

    char test_file[PATH_MAX];
    snprintf (test_file, sizeof(test_file), "%s/files/smoothbump/initial_grid.x", test_dir);

    // Check that file exists
    if (access (test_file, R_OK) != 0) 
    {
        cr_skip_test ("Test file not found at %s", test_file);
    }

    // Read quantities from plot3D file
    long nblks = read_nblocks_from_plot3D (test_file);
    long *npts_x = read_npts_x_from_plot3D ("2D", test_file, nblks);
    long *npts_y = read_npts_y_from_plot3D ("2D", test_file, nblks);

    // Test number of blocks, number of points in
    // X-direction and Y-direction
    cr_assert_eq (nblks, 1L, "Expected 1 block but got %ld", nblks);
    cr_assert_not_null (npts_x, "read_npts_x_from_plot3D returned NULL");
    cr_assert_eq (npts_x[0], 121L, "Expected 121 points in X direction but got %ld", npts_x[0]);
    cr_assert_not_null (npts_y, "read_npts_y_from_plot3D returned NULL");
    cr_assert_eq (npts_y[0], 31L, "Expected 31 points in Y direction but got %ld", npts_y[0]);

    free_1D_long_array ("npts_x", npts_x);
    free_1D_long_array ("npts_y", npts_y);
}


/* Test for (x, y) boundary information file reading functions */
Test (file_operations, read_xy_2D)
{
    // Compute relative path from test source file
    char src_path[PATH_MAX];    
    strncpy (src_path, __FILE__, PATH_MAX - 1);
    src_path[PATH_MAX - 1] = '\0';
    char *test_dir = dirname (src_path);

    // Path to smoothbump X-boundaries file
    char sbump_dir[PATH_MAX];
    snprintf (sbump_dir, sizeof(sbump_dir), "%s/files/smoothbump", test_dir);

    // Save current cwd and chdir into smoothbump directory so test finds "x_boundaries.dat"
    char oldcwd[PATH_MAX];
    cr_assert_not_null (getcwd (oldcwd, sizeof(oldcwd)), "getcwd failed");
    cr_assert_eq (chdir(sbump_dir), 0, "chdir to %s failed", sbump_dir);
   
    
    // X-boundaries
    char sbump_x_file[PATH_MAX];
    snprintf (sbump_x_file, sizeof(sbump_x_file), "%s/x_boundaries.dat", sbump_dir);

    // Check that file exists
    if (access (sbump_x_file, R_OK) != 0)
    {
        cr_skip_test ("Test file not found at %s", sbump_x_file);
    }

    // Call function for reading X-boundaries
    point_2D **x_bounds = read_x_bounds_2D (121);
    cr_assert_not_null (x_bounds, "read_x_bounds_2D returned NULL");

    // Check values
    long double tol = 1e-12L;
    int i11[9] = {0, 0, 0, 0, 1, 1, 1, 1};
    int i12[9] = {0, 0, 60, 60, 60, 60, 120, 120};
    long double ref_values_1[8] = {-1.5, 0.0, 0.0, 0.15, 0.0, 0.8, 1.5, 0.8};
    int i1, j1, k1;
    for (i1 = 0; i1 < 4; i1++)
    {
        j1 = 2 * i1;
        k1 = j1 + 1;
        cr_assert_float_eq (x_bounds[i11[j1]][i12[j1]].x, ref_values_1[j1], tol, "mismatch in x_bounds[%d][%d].x: expected %LF, but got %LF", 
                            i11[j1], i12[j1], ref_values_1[j1], x_bounds[i11[j1]][i12[j1]].x);
        cr_assert_float_eq (x_bounds[i11[k1]][i12[k1]].y, ref_values_1[k1], tol, "mismatch in x_bounds[%d][%d].y: expected %LF, but got %LF", 
                            i11[k1], i12[k1], ref_values_1[k1], x_bounds[i11[k1]][i12[k1]].y);
    }


    // Y-boundaries
    char sbump_y_file[PATH_MAX];
    snprintf (sbump_y_file, sizeof(sbump_y_file), "%s/y_boundaries.dat", sbump_dir);

    // Check that file exists
    if (access (sbump_y_file, R_OK) != 0)
    {
        cr_skip_test ("Test file not found at %s", sbump_y_file);
    }

    // Call function for reading Y-boundaries
    point_2D **y_bounds = read_y_bounds_2D (31);
    cr_assert_not_null (y_bounds, "read_y_bounds_2D returned NULL");

    // Check values
    int i21[9] = {0, 0, 0, 0, 1, 1, 1, 1};
    int i22[9] = {0, 0, 15, 15, 15, 15, 30, 30};
    long double ref_values_2[8] = {-1.5, 0.0, -1.5, 0.4, 1.5, 0.4, 1.5, 0.8};
    int i2, j2, k2;
    for (i2 = 0; i2 < 4; i2++)
    {
        j2 = 2 * i2;
        k2 = j2 + 1;
        cr_assert_float_eq (y_bounds[i21[j2]][i22[j2]].x, ref_values_2[j2], tol, "mismatch in y_bounds[%d][%d].x: expected %LF, but got %LF", 
                            i21[j2], i22[j2], ref_values_2[j2], y_bounds[i21[j2]][i22[j2]].x);
        cr_assert_float_eq (y_bounds[i21[k2]][i22[k2]].y, ref_values_2[k2], tol, "mismatch in y_bounds[%d][%d].y: expected %LF, but got %LF", 
                            i21[k2], i22[k2], ref_values_2[k2], y_bounds[i21[k2]][i22[k2]].y);
    }


    // X-normals
    char sbump_xn_file[PATH_MAX];
    snprintf (sbump_xn_file, sizeof(sbump_xn_file), "%s/xbound_normals.dat", sbump_dir);

    // Check that file exists
    if (access (sbump_xn_file, R_OK) != 0)
    {
        cr_skip_test ("Test file not found at %s", sbump_xn_file);
    }

    // Call function for reading X-normals
    point_2D **x_normals = read_xbound_normals_2D (121);
    cr_assert_not_null (x_normals, "read_xbound_normals_2D returned NULL");

    // Check values
    int i31[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};
    int i32[13] = {0, 0, 30, 30, 60, 60, 90, 90, 120, 120, 0, 0};
    long double ref_values_3[13] = {0.0, -1.0, 0.0045881492887020, -0.9999894743876582,
                                    0.0000000000000036, -1.0, -0.0045881492887021, -0.9999894743876582,
                                    0.0, -1.0, 0.0, 1.0};
    int i3, j3, k3;
    for (i3 = 0; i3 < 6; i3++)
    {
        j3 = 2 * i3;
        k3 = j3 + 1;
        cr_assert_float_eq (x_normals[i31[j3]][i32[j3]].x, ref_values_3[j3], tol, "mismatch in x_normals[%d][%d].x: expected %LF, but got %LF", 
                            i31[j3], i32[j3], ref_values_3[j3], x_normals[i31[j3]][i32[j3]].x);
        cr_assert_float_eq (x_normals[i31[k3]][i32[k3]].y, ref_values_3[k3], tol, "mismatch in x_normals[%d][%d].y: expected %LF, but got %LF", 
                            i31[k3], i32[k3], ref_values_3[k3], x_normals[i31[k3]][i32[k3]].y);
    }


    // Y-normals
    char sbump_yn_file[PATH_MAX];
    snprintf (sbump_yn_file, sizeof(sbump_yn_file), "%s/ybound_normals.dat", sbump_dir);

    // Check that file exists
    if (access (sbump_yn_file, R_OK) != 0)
    {
        cr_skip_test ("Test file not found at %s", sbump_yn_file);
    }

    // Call function for reading Y-normals
    point_2D **y_normals = read_ybound_normals_2D (31);
    cr_assert_not_null (y_normals, "read_ybound_normals_2D returned NULL");

    // Check values
    int i41[5] = {0, 0, 1, 1};
    long double ref_values_4[5] = {-1.0, 0.0, 1.0, 0.0};
    int i4, j4, k4;
    for (i4 = 0; i4 < 2; i4++)
    {
        j4 = 2 * i4;
        k4 = j4 + 1;
        cr_assert_float_eq (y_normals[i41[j4]][0].x, ref_values_4[j4], tol, "mismatch in y_normals[%d][%d].x: expected %LF, but got %LF",
                            i41[j4], 0, ref_values_4[j4], y_normals[i41[j4]][0].x);
        cr_assert_float_eq (y_normals[i41[k4]][0].y, ref_values_4[k4], tol, "mismatch in y_normals[%d][%d].y: expected %LF, but got %LF",
                            i41[k4], 0, ref_values_4[k4], y_normals[i41[k4]][0].y);
    }


    // Free memory and restore cwd
    free_2D_point_2D_array ("x_bounds", 2, x_bounds);
    free_2D_point_2D_array ("y_bounds", 2, y_bounds);
    free_2D_point_2D_array ("x_normals", 2, x_normals);
    free_2D_point_2D_array ("y_normals", 2, y_normals);
    cr_assert_eq (chdir (oldcwd), 0, "chdir back failed");
}


/* Test for x boundary information file writing function */
Test(file_operations, write_bounds_2D_x)
{
    // Compute test directory relative to this source file
    char src[PATH_MAX];
    strncpy(src, __FILE__, PATH_MAX - 1);
    src[PATH_MAX - 1] = '\0';
    char *tests_dir = dirname(src);

    // Build path to reference file
    char ref_file[PATH_MAX];
    snprintf(ref_file, sizeof(ref_file), "%s/files/smoothbump/x_boundaries.dat", tests_dir);

    if (access(ref_file, R_OK) != 0) {
        cr_skip_test("Reference file not found at %s", ref_file);
    }
   
    // Open reference file
    FILE *ref_f = fopen(ref_file, "r");
    cr_assert_not_null(ref_f, "Cannot open reference file %s", ref_file);
    rewind (ref_f);

    // Generate boundary array data
    int nx = 121;
    point_2D **bounds = allocate_2D_point_2D_array ("bounds", 2, nx);
    cr_assert_not_null (bounds, "Failed to allocate bounds array");

    // Read reference file into bounds array
    for (int i = 0; i < nx; i++)
    {
        int ret = fscanf(ref_f, "%Lf %Lf %Lf %Lf\n",
                        &bounds[0][i].x, &bounds[0][i].y,
                        &bounds[1][i].x, &bounds[1][i].y);
        cr_assert_eq(ret, 4, "Failed to read line %d from reference file %s", i, ref_file);
    }
    fclose (ref_f);

    // Create temporary directory for output
    char tmpdir_template[] = "/tmp/gridgen_write_bounds_XXXXXX";
    char *tmpdir = mkdtemp(tmpdir_template);
    cr_assert_not_null(tmpdir, "mkdtemp failed");

    char out_file[PATH_MAX];
    snprintf(out_file, sizeof(out_file), "%s/x_boundaries_out.dat", tmpdir);

    // Call write_bounds_2D (writes to relative filename in current dir)
    write_bounds_2D(out_file, nx, bounds);
    cr_assert_eq(access(out_file, R_OK), 0, "write_bounds_2D output file %s not created", out_file);

    // Compare output file with reference file line-by-line
    FILE *out_f = fopen(out_file, "r");
    cr_assert_not_null(out_f, "Cannot open output file %s", out_file);
    rewind (ref_f);
    ref_f = fopen(ref_file, "r");
    cr_assert_not_null (ref_f, "Cannot reopen reference file %s", ref_file);

    // Test output
    char ref_line[512], out_line[512];
    int line_num = 0;
    long double tol = 1E-12L;
    while (fgets(ref_line, sizeof(ref_line), ref_f) != NULL &&
           fgets(out_line, sizeof(out_line), out_f) != NULL) 
    {
        // Parse and compare each numeric value with tolerance for floating-point rounding
        long double ref_x0, ref_y0, ref_x1, ref_y1;
        long double out_x0, out_y0, out_x1, out_y1;
        
        int ref_ret = sscanf(ref_line, "%Lf %Lf %Lf %Lf",
                            &ref_x0, &ref_y0, &ref_x1, &ref_y1);
        int out_ret = sscanf(out_line, "%Lf %Lf %Lf %Lf",
                            &out_x0, &out_y0, &out_x1, &out_y1);
        
        cr_assert_eq(ref_ret, 4, "Reference line %d parse error", line_num);
        cr_assert_eq(out_ret, 4, "Output line %d parse error", line_num);

        cr_assert(fabsl(out_x0 - ref_x0) < tol, "Line %d, bounds[0][%d].x: Expected %Lf but got %Lf", line_num, line_num, ref_x0, out_x0);
        cr_assert(fabsl(out_y0 - ref_y0) < tol, "Line %d, bounds[0][%d].y: Expected %Lf but got %Lf", line_num, line_num, ref_y0, out_y0);
        cr_assert(fabsl(out_x1 - ref_x1) < tol, "Line %d, bounds[1][%d].x: Expected %Lf but got %Lf", line_num, line_num, ref_x1, out_x1);
        cr_assert(fabsl(out_y1 - ref_y1) < tol, "Line %d, bounds[1][%d].y: Expected %Lf but got %Lf", line_num, line_num, ref_y1, out_y1);

        ++line_num;
    }

    // Ensure both files read the same number of lines
    cr_assert_eq(line_num, nx, "Output file has %d lines, expected %d", line_num, nx);

    // Cleanup
    fclose(out_f);
    fclose(ref_f);
    free_2D_point_2D_array("bounds", 2, bounds);
    unlink(out_file);
    rmdir(tmpdir);
}


/* Test for y boundary information file writing function */
Test(file_operations, write_bounds_2D_y)
{
    // Compute test directory relative to this source file
    char src[PATH_MAX];
    strncpy(src, __FILE__, PATH_MAX - 1);
    src[PATH_MAX - 1] = '\0';
    char *tests_dir = dirname(src);

    // Build path to reference file
    char ref_file[PATH_MAX];
    snprintf(ref_file, sizeof(ref_file), "%s/files/smoothbump/y_boundaries.dat", tests_dir);

    if (access(ref_file, R_OK) != 0) {
        cr_skip_test("Reference file not found at %s", ref_file);
    }
   
    // Open reference file
    FILE *ref_f = fopen(ref_file, "r");
    cr_assert_not_null(ref_f, "Cannot open reference file %s", ref_file);
    rewind (ref_f);

    // Generate boundary array data
    int ny = 31;
    point_2D **bounds = allocate_2D_point_2D_array ("bounds", 2, ny);
    cr_assert_not_null (bounds, "Failed to allocate bounds array");

    // Read reference file into bounds array
    for (int i = 0; i < ny; i++)
    {
        int ret = fscanf(ref_f, "%Lf %Lf %Lf %Lf\n",
                        &bounds[0][i].x, &bounds[0][i].y,
                        &bounds[1][i].x, &bounds[1][i].y);
        cr_assert_eq(ret, 4, "Failed to read line %d from reference file %s", i, ref_file);
    }
    fclose (ref_f);

    // Create temporary directory for output
    char tmpdir_template[] = "/tmp/gridgen_write_bounds_y_XXXXXX";
    char *tmpdir = mkdtemp(tmpdir_template);
    cr_assert_not_null(tmpdir, "mkdtemp failed");

    char out_file[PATH_MAX];
    snprintf(out_file, sizeof(out_file), "%s/y_boundaries_out.dat", tmpdir);

    // Call write_bounds_2D (writes to relative filename in current dir)
    write_bounds_2D(out_file, ny, bounds);
    cr_assert_eq(access(out_file, R_OK), 0, "write_bounds_2D output file %s not created", out_file);

    // Open files
    FILE *out_f = fopen(out_file, "r");
    cr_assert_not_null(out_f, "Cannot open output file %s", out_file);
    rewind (ref_f);
    ref_f = fopen(ref_file, "r");
    cr_assert_not_null (ref_f, "Cannot reopen reference file %s", ref_file);

    // Test output
    char ref_line[512], out_line[512];
    int line_num = 0;
    long double tol = 1E-12L;
    while (fgets(ref_line, sizeof(ref_line), ref_f) != NULL &&
           fgets(out_line, sizeof(out_line), out_f) != NULL) 
    {
        // Parse and compare each numeric value with tolerance for floating-point rounding
        long double ref_x0, ref_y0, ref_x1, ref_y1;
        long double out_x0, out_y0, out_x1, out_y1;
        
        int ref_ret = sscanf(ref_line, "%Lf %Lf %Lf %Lf",
                            &ref_x0, &ref_y0, &ref_x1, &ref_y1);
        int out_ret = sscanf(out_line, "%Lf %Lf %Lf %Lf",
                            &out_x0, &out_y0, &out_x1, &out_y1);
        
        cr_assert_eq (ref_ret, 4, "Reference line %d parse error", line_num);
        cr_assert_eq (out_ret, 4, "Output line %d parse error", line_num);

        cr_assert (fabsl(out_x0 - ref_x0) < tol, "Line %d, bounds[0][%d].x: Expected %Lf but got %Lf", line_num, line_num, ref_x0, out_x0);
        cr_assert (fabsl(out_y0 - ref_y0) < tol, "Line %d, bounds[0][%d].y: Expected %Lf but got %Lf", line_num, line_num, ref_y0, out_y0);
        cr_assert (fabsl(out_x1 - ref_x1) < tol, "Line %d, bounds[1][%d].x: Expected %Lf but got %Lf", line_num, line_num, ref_x1, out_x1);
        cr_assert (fabsl(out_y1 - ref_y1) < tol, "Line %d, bounds[1][%d].y: Expected %Lf but got %Lf", line_num, line_num, ref_y1, out_y1);

        ++line_num;
    }

    // Ensure both files read the same number of lines
    cr_assert_eq (line_num, ny, "Output file has %d lines, expected %d", line_num, ny);

    // Cleanup
    fclose(out_f);
    fclose(ref_f);
    free_2D_point_2D_array("bounds", 2, bounds);
    unlink(out_file);
    rmdir(tmpdir);
}


/* Test for x normal information file writing function */
Test(file_operations, write_normals_2D_x)
{
    // Compute test directory relative to this source file
    char src[PATH_MAX];
    strncpy(src, __FILE__, PATH_MAX - 1);
    src[PATH_MAX - 1] = '\0';
    char *tests_dir = dirname(src);

    // Build path to reference file
    char ref_file[PATH_MAX];
    snprintf(ref_file, sizeof(ref_file), "%s/files/smoothbump/x_boundaries.dat", tests_dir);

    if (access(ref_file, R_OK) != 0) {
        cr_skip_test("Reference file not found at %s", ref_file);
    }
   
    // Open reference file
    FILE *ref_f = fopen(ref_file, "r");
    cr_assert_not_null(ref_f, "Cannot open reference file %s", ref_file);
    rewind (ref_f);

    // Generate boundary array data
    int nx = 121;
    point_2D **bounds = allocate_2D_point_2D_array ("bounds", 2, nx);
    cr_assert_not_null (bounds, "Failed to allocate bounds array");

    // Read reference file into bounds array
    for (int i = 0; i < nx; i++)
    {
        int ret = fscanf(ref_f, "%Lf %Lf %Lf %Lf\n",
                        &bounds[0][i].x, &bounds[0][i].y,
                        &bounds[1][i].x, &bounds[1][i].y);
        cr_assert_eq(ret, 4, "Failed to read line %d from reference file %s", i, ref_file);
    }
    fclose (ref_f);

    // Create temporary directory for output
    char tmpdir_template[] = "/tmp/gridgen_write_bounds_XXXXXX";
    char *tmpdir = mkdtemp(tmpdir_template);
    cr_assert_not_null(tmpdir, "mkdtemp failed");

    char out_file[PATH_MAX];
    snprintf(out_file, sizeof(out_file), "%s/x_boundaries_out.dat", tmpdir);

    // Call write_bounds_2D (writes to relative filename in current dir)
    write_bounds_2D(out_file, nx, bounds);
    cr_assert_eq(access(out_file, R_OK), 0, "write_bounds_2D output file %s not created", out_file);

    // Compare output file with reference file line-by-line
    FILE *out_f = fopen(out_file, "r");
    cr_assert_not_null(out_f, "Cannot open output file %s", out_file);
    rewind (ref_f);
    ref_f = fopen(ref_file, "r");
    cr_assert_not_null (ref_f, "Cannot reopen reference file %s", ref_file);

    // Test output
    char ref_line[512], out_line[512];
    int line_num = 0;
    long double tol = 1E-12L;
    while (fgets(ref_line, sizeof(ref_line), ref_f) != NULL &&
           fgets(out_line, sizeof(out_line), out_f) != NULL) 
    {
        // Parse and compare each numeric value with tolerance for floating-point rounding
        long double ref_x0, ref_y0, ref_x1, ref_y1;
        long double out_x0, out_y0, out_x1, out_y1;
        
        int ref_ret = sscanf(ref_line, "%Lf %Lf %Lf %Lf",
                            &ref_x0, &ref_y0, &ref_x1, &ref_y1);
        int out_ret = sscanf(out_line, "%Lf %Lf %Lf %Lf",
                            &out_x0, &out_y0, &out_x1, &out_y1);
        
        cr_assert_eq(ref_ret, 4, "Reference line %d parse error", line_num);
        cr_assert_eq(out_ret, 4, "Output line %d parse error", line_num);

        cr_assert(fabsl(out_x0 - ref_x0) < tol, "Line %d, bounds[0][%d].x: Expected %Lf but got %Lf", line_num, line_num, ref_x0, out_x0);
        cr_assert(fabsl(out_y0 - ref_y0) < tol, "Line %d, bounds[0][%d].y: Expected %Lf but got %Lf", line_num, line_num, ref_y0, out_y0);
        cr_assert(fabsl(out_x1 - ref_x1) < tol, "Line %d, bounds[1][%d].x: Expected %Lf but got %Lf", line_num, line_num, ref_x1, out_x1);
        cr_assert(fabsl(out_y1 - ref_y1) < tol, "Line %d, bounds[1][%d].y: Expected %Lf but got %Lf", line_num, line_num, ref_y1, out_y1);

        ++line_num;
    }

    // Ensure both files read the same number of lines
    cr_assert_eq(line_num, nx, "Output file has %d lines, expected %d", line_num, nx);

    // Cleanup
    fclose(out_f);
    fclose(ref_f);
    free_2D_point_2D_array("bounds", 2, bounds);
    unlink(out_file);
    rmdir(tmpdir);
}


/* Test for function reading a single block plot3D file 
   containing a 2D grid */
Test (file_operations, read_2D_singleblock_plot3D)
{
    // Compute path to test file relative to this source file
    char src_path[PATH_MAX];
    strncpy (src_path, __FILE__, PATH_MAX - 1);
    src_path[PATH_MAX - 1] = '\0';
    char *test_dir = dirname (src_path);

    // Determine if file exists
    char test_file[PATH_MAX];
    snprintf (test_file, sizeof(test_file), "%s/files/smoothbump/initial_grid.x", test_dir);
    if (access(test_file, R_OK) != 0)
    {
        cr_skip_test("Test file not found at %s", test_file);
    }

    // Read grid and check if it was allocated successfully
    point_2D **grid = read_2D_singleblock_plot3D(test_file);
    cr_assert_not_null(grid, "read_2D_singleblock_plot3D returned NULL");

    // Test output
    long double tol = 1E-6L;
    cr_assert_float_eq(grid[0][0].x, -1.5, tol, "grid[0][0].x: Expected -1.5, but got %LF", grid[0][0].x);
    cr_assert_float_eq(grid[120][0].x, 1.5, tol, "grid[120][0].x: Expected 1.5, but got %LF", grid[120][0].x);
    cr_assert_float_eq(grid[0][0].y, 0.0, tol, "grid[0][0].y: Expected 0.0, but got %LF", grid[0][0].y);
    cr_assert_float_eq(grid[120][0].y, 0.0, tol, "grid[120][0].y: Expected 0.0, but got %LF", grid[120][0].y);
    cr_assert_float_eq(grid[120][30].y, 0.8, tol, "grid[120][30].y: Expected 0.8, but got %LF", grid[120][30].y);

    // Free allocated memory
    free_2D_point_2D_array("grid", 121, grid);
}


/* Test for function writing a single block plot3D file 
   Write a synthetic 3x3 unit-square grid to a plot3D file, 
   read it back, and verify all values are preserved */
Test(file_operations, write_2D_singleblock_plot3D)
{
    // Build a 3x3 uniform grid on a [0, 1] x [0, 1] domain
    // with spacing 0.5
    int i, j, nx = 3, ny = 3;
    point_2D **grid = make_unit_square_grid(nx, ny);
    cr_assert_not_null(grid, "allocate_2D_point_2D_array returned NULL");

    // Write grid to a temporary file
    char tmpdir_template[] = "/tmp/gridgen_plot3D_XXXXXX";
    char *tmpdir = mkdtemp(tmpdir_template);
    cr_assert_not_null(tmpdir, "mkdtemp failed");

    char out_file[PATH_MAX];
    snprintf(out_file, sizeof(out_file), "%s/square_grid.x", tmpdir);

    write_2D_singleblock_plot3D(out_file, nx, ny, grid);
    cr_assert_eq(access(out_file, R_OK), 0, "Output file not created: %s", out_file);

    // Read the file back
    point_2D **grid_read = read_2D_singleblock_plot3D(out_file);
    cr_assert_not_null(grid_read, "read_2D_singleblock_plot3D returned NULL");

    // Compare every grid point
    long double tol = 1E-12L;
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            cr_assert_float_eq(grid_read[i][j].x, grid[i][j].x, tol, "grid[%d][%d].x: expected %LF, but got %LF", i, j, grid[i][j].x, grid_read[i][j].x);
            cr_assert_float_eq(grid_read[i][j].y, grid[i][j].y, tol, "grid[%d][%d].y: expected %LF, but got %LF", i, j, grid[i][j].y, grid_read[i][j].y);
        }
    }

    // Clean up
    free_2D_point_2D_array("grid", nx, grid);
    free_2D_point_2D_array("grid_read", nx, grid_read);
    unlink(out_file);
    rmdir(tmpdir);
}


/* Test for write_2D_point_2D: Write a unit-square grid to 
   two separate .dat files (x and y), read them back, and 
   verify all values are preserved */ 
Test(file_operations, write_2D_point_2D)
{
    // Local variables
    int         i, j, nx = 21, ny = 21;
    long double val_x, val_y, tol = 1E-15L;

    // Build a nx x ny uniform grid on a [0, 1] x [0, 1] domain
    point_2D **grid = make_unit_square_grid(nx, ny);
    cr_assert_not_null(grid, "allocate_2D_point_2D_array returned NULL");

    // Create a temporary directory for output files
    char tmpdir_template[] = "/tmp/gridgen_point2D_XXXXXX";
    char *tmpdir = mkdtemp(tmpdir_template);
    cr_assert_not_null(tmpdir, "mkdtemp failed");

    char x_file[PATH_MAX], y_file[PATH_MAX];
    snprintf(x_file, sizeof(x_file), "%s/grid_x.dat", tmpdir);
    snprintf(y_file, sizeof(y_file), "%s/grid_y.dat", tmpdir);

    // Write grid components to files
    write_2D_point_2D(x_file, y_file, nx, ny, grid);
    cr_assert_eq(access(x_file, R_OK), 0, "x output file not created: %s", x_file);
    cr_assert_eq(access(y_file, R_OK), 0, "y output file not created: %s", y_file);

    // Read back and compare
    FILE *fx = fopen(x_file, "r");
    FILE *fy = fopen(y_file, "r");
    cr_assert_not_null(fx, "Cannot open x output file %s", x_file);
    cr_assert_not_null(fy, "Cannot open y output file %s", y_file);

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            cr_assert_eq(fscanf(fx, "%LF", &val_x), 1, "fscanf failed at [%d][%d] in x file", i, j);
            cr_assert_eq(fscanf(fy, "%LF", &val_y), 1, "fscanf failed at [%d][%d] in y file", i, j);
            cr_assert_float_eq(val_x, grid[i][j].x, tol, "grid[%d][%d].x: expected %LF, but got %LF", i, j, grid[i][j].x, val_x);
            cr_assert_float_eq(val_y, grid[i][j].y, tol, "grid[%d][%d].y: expected %LF, but got %LF", i, j, grid[i][j].y, val_y);
        }
    }

    // Clean up
    fclose(fx);
    fclose(fy);
    free_2D_point_2D_array("grid", nx, grid);
    unlink(x_file);
    unlink(y_file);
    rmdir(tmpdir);
}

#ifdef HAVE_LIBXML2
/* Find the first child element node with the given tag name
   Helper function for write_2D_singleblock_vts test */
static xmlNodePtr find_child(xmlNodePtr parent, const char *name)
{
    xmlNodePtr node;
    for (node = parent->children; node != NULL; node = node->next)
    {
        if (node->type == XML_ELEMENT_NODE &&
            xmlStrcmp(node->name, (const xmlChar *)name) == 0)
            return node;
    }
    return NULL;
}

/* Test for write_2D_singleblock_vts: write a unit-square grid 
   to a .vts file, parse with libxml2 and verify coordinate values */
Test(file_operations, write_2D_singleblock_vts)
{
    // Local variables
    int                 i, j, nx = 21, ny = 21;
    long double         rx, ry, tol = 1E-12L;

    // Build a nx x ny uniform grid on [0, 1] x [0, 1] domain
    point_2D **grid = make_unit_square_grid(nx, ny);
    cr_assert_not_null(grid, "allocate_2D_point_2D_array returned NULL");

    // Create a temporary directory and output file path
    char tmpdir_template[] = "/tmp/gridgen_vts_XXXXXX";
    char *tmpdir = mkdtemp(tmpdir_template);
    cr_assert_not_null(tmpdir, "mkdtemp failed");

    char vts_file[PATH_MAX];
    snprintf(vts_file, sizeof(vts_file), "%s/grid.vts", tmpdir);

    // Write the grid to a .vts file
    write_2D_singleblock_vts(vts_file, nx, ny, grid, 1);
    cr_assert_eq(access(vts_file, R_OK), 0, "Output file not created: %s", vts_file);

    // Parse the XML output file
    xmlDoc *doc = xmlReadFile(vts_file, NULL, 0);
    cr_assert_not_null(doc, "xmlReadFile failed to parse %s", vts_file);

    // Navigate to the DataArray node:
    // VTKFile -> StructuredGrid -> Piece -> Points -> Dataarray
    xmlNodePtr root = xmlDocGetRootElement(doc);
    xmlNodePtr struct_grid = find_child(root, "StructuredGrid");
    cr_assert_not_null(struct_grid, "StructuredGrid node not found");

    xmlNodePtr piece = find_child(struct_grid, "Piece");
    cr_assert_not_null(piece, "Piece node not found");

    xmlNodePtr points_node = find_child(piece, "Points");
    cr_assert_not_null(points_node, "Points node not found");

    xmlNodePtr data_array = find_child(points_node, "DataArray");
    cr_assert_not_null(data_array, "DataArray node not found");

    // Verify WholeExtent and Piece Extent attributes
    char expected_extent[64];
    snprintf(expected_extent, sizeof(expected_extent), "0 %d 0 %d 0 0", nx - 1, ny - 1);

    xmlChar *whole_extent = xmlGetProp(struct_grid, (const xmlChar *)"WholeExtent");
    xmlChar *piece_extent = xmlGetProp(piece, (const xmlChar *)"Extent");

    cr_assert_str_eq((char *)whole_extent, expected_extent, "WholeExtent: expected '%s', got '%s'", 
                     expected_extent, (char *)whole_extent);
    cr_assert_str_eq((char *)piece_extent, expected_extent, "Piece Extent: expected '%s', got '%s'", 
                     expected_extent, (char *)piece_extent);

    xmlFree(whole_extent);
    xmlFree(piece_extent);

    // Parse coordinate values from DataArray and compare against
    // original grid values
    xmlChar *content = xmlNodeGetContent(data_array);
    cr_assert_not_null(content, "DataArray has no text content");

    char *ptr = (char *)content;
    char *end;

    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            rx = strtold(ptr, &end); ptr = end;
            ry = strtold(ptr, &end); ptr = end;
            strtold(ptr, &end); ptr = end; // skip z = 0.0 
        
            cr_assert_float_eq(rx, grid[i][j].x, tol, "grid[%d][%d].x: expected %LF, but got %LF", i, j, grid[i][j].x, rx);
            cr_assert_float_eq(ry, grid[i][j].y, tol, "grid[%d][%d].y: expected %LF, but got %LF", i, j, grid[i][j].y, ry);
        }
    }

    // Clean up
    xmlFree(content);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    free_2D_point_2D_array("grid", nx, grid);
    unlink(vts_file);
    rmdir(tmpdir);
}


/* Test for write_point_2D_to_vts: write a unit-square grid with a 
   manufactured vector field to a .vts file, parse with libxml2, and 
   veify the PointData vector values */
Test(file_operations, write_point_2D_to_vts)
{
    // Local variables
    int                     i, j, nx = 5, ny = 5;
    long double             tol = 1E-12L;
    const long double       pi = acosl(-1.0L);

    // Build a nx x ny uniform grid on [0, 1] x [0, 1] domain
    point_2D **grid = make_unit_square_grid(nx, ny);
    cr_assert_not_null(grid, "allocate_2D_point_2D_array returned NULL");

    // Build a manufactured vector field: f(x, y) = (sin(pi*x), cos(pi*y))
    point_2D **array = allocate_2D_point_2D_array("array", nx, ny);
    cr_assert_not_null(array, "allocate_2D_point_2D_array returned NULL for array");

    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            array[i][j].x = sinl(pi * grid[i][j].x);
            array[i][j].y = cosl(pi * grid[i][j].y);
        }
    }

    // Create a temporary director and output file path
    char tmpdir_template[] = "/tmp/gridgen_point2D_vts_XXXXXX";
    char *tmpdir = mkdtemp(tmpdir_template);
    cr_assert_not_null(tmpdir, "mkdtemp failed");

    char vts_file[PATH_MAX];
    snprintf(vts_file, sizeof(vts_file), "%s/vector.vts", tmpdir);

    // Write vector field to .vts file
    write_point_2D_to_vts(vts_file, nx, ny, grid, "f", array);
    cr_assert_eq(access(vts_file, R_OK), 0, "Output file not created: %s", vts_file);

    // Parse the XML output file
    xmlDoc *doc = xmlReadFile(vts_file, NULL, 0);
    cr_assert_not_null(doc, "xmlReadFile failed to parse %s", vts_file);

    // Navigate to the PointData DataArray node:
    // VTKFile -> StructuredGrid -> Piece -> PointData -> DataArray
    xmlNodePtr root = xmlDocGetRootElement(doc);
    xmlNodePtr struct_grid = find_child(root, "StructuredGrid");
    cr_assert_not_null(struct_grid, "StructuredGrid node not found");

    xmlNodePtr piece = find_child(struct_grid, "Piece");
    cr_assert_not_null(piece, "Piece node not found");

    xmlNodePtr point_data = find_child(piece, "PointData");
    cr_assert_not_null(point_data, "PointData node not found");

    // Verify the Vectors attribute contains the correct field name
    xmlChar *vectors_attr = xmlGetProp(point_data, (const xmlChar *)"Vectors");
    cr_assert_not_null(vectors_attr, "Vectors attribute not found on PointData node");
    cr_assert_str_eq((char *)vectors_attr, "f", 
                     "Vectors attribute: expected 'f', got '%s'", (char *)vectors_attr);
    xmlFree(vectors_attr);

    xmlNodePtr data_array = find_child(point_data, "DataArray");
    cr_assert_not_null(data_array, "DataArray node not found in PointData");

    // Parse vector values from DataArray and compare against
    // manufactured field
    xmlChar *content = xmlNodeGetContent(data_array);
    cr_assert_not_null(content, "DataArray has no text content");

    char *ptr = (char *)content;
    char *end;
    long double rx, ry, rz;

    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            rx = strtold(ptr, &end); ptr = end;
            ry = strtold(ptr, &end); ptr = end;
            rz = strtold(ptr, &end); ptr = end;

            cr_assert_float_eq(rx, array[i][j].x, tol,
                               "array[%d][%d].x: expected %LF, but got %LF", i, j, array[i][j].x, rx);
            cr_assert_float_eq(ry, array[i][j].y, tol,
                               "array[%d][%d].y: expected %LF, but got %LF", i, j, array[i][j].y, ry);
            cr_assert_float_eq(rz, 0.0L, tol,
                               "array[%d][%d].z: expected 0.0, but got %LF", i, j, rz);
        }
    }

    // Clean-up
    xmlFree(content);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    free_2D_point_2D_array("grid", nx, grid);
    free_2D_point_2D_array("array", nx, array);
    unlink(vts_file);
    rmdir(tmpdir);
}
#endif
