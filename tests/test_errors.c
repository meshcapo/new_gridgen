#include <criterion/criterion.h>
#include <criterion/parameterized.h>
#include <criterion/redirect.h>
#include "errors.h"


/* Test for function fatal_error - captures exit status and prints 
   error message to screen when running test suite */
Test (errors, fatal_error_exit_cond, .exit_code = EXIT_FAILURE)
{
    fatal_error ("This is a test fatal error\n");
}


/* Tests for function file_extension_error */
Test (errors, file_extension_error_valid_x, .exit_code = EXIT_SUCCESS)
{
    file_extension_error (".x");
}

Test (errors, file_extension_error_valid_vts, .exit_code = EXIT_SUCCESS)
{
    file_extension_error (".vts");
}

Test (errors, file_extension_error_invalid, .exit_code = EXIT_FAILURE)
{
    file_extension_error (".dat");
}