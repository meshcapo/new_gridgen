#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/*
   Function to print an error message and
   exit with failure

   Input parameters: error_message  - the error string
*/
void fatal_error (char *error_message)
{
    printf ("FATAL ERROR: %s\n", error_message);
    exit (EXIT_FAILURE);
}






/*
   Print an error message if array dynamic memory
   allocation fails

   Input parameters: array_name - array descriptor
                     irow       - row index (for 2D arrays)
*/
void array_allocation_msg (char *array_name, int irow)
{
    if (irow == 0) {
        printf ("Memory allocation failed for %s\n", array_name);
    } else {
        printf ("Memory allocation failed for %s[%d]\n", array_name, irow);
    }
}






/*
   Print an error message if the extension of the
   file containing the initial grid is unsupported

   Input parameters: ext    - file extension
*/
void file_extension_error (char *ext)
{
    if (strcmp(ext, ".x") != 0 && strcmp(ext, ".vts") != 0)
    {
        fatal_error ("Initial grids can only be provided in PLOT3D/VTK format with extensions .x/.vts.");
    }
}
