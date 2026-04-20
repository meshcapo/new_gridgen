#include <stdio.h>
#include <string.h>
#include "types.h"
#include "errors.h"
#include "aux_functions.h"



/*
   Determine the extension of a file - used for reading
   the initial grid, if provided. Raise a fatal error if
   the file has no extension

   Input parameters: filename   - name of the file containing the initial grid
*/
char *get_file_extension (char *filename)
{
    /* Return the file extension */
    char                    *extension;

    /* Local variables */
    char                    dot = '.', *char_ptr;


    /* Initialization */
    extension                           = "";

    /* Find location of '.' in the filename */
    char_ptr                            = strchr(filename, dot);

    /* Raise fatal error if '.' is not found in the filename */
    if (char_ptr == NULL)
    {
        fatal_error ("No extension found for the file containing the initial grid.");
    }
    else
    {
        extension                           = char_ptr; // char_ptr contains the part of the filename starting from '.'
    }


    return extension;
}






/*
   Read the number of blocks from a plot3D file

   Input parameters: filename   - name of the initial grid file
*/
long read_nblocks_from_plot3D (char *filename)
{
    /* Return nblks */
    long                    nblks;

    /* Local variables */
    FILE                    *fptr;


    /* Read only the first line from the incoming
       .x file to get the number of blocks in
       structured grid */
    fptr                                = fopen (filename, "r");
    fscanf (fptr, "%ld", &nblks);
    fclose (fptr);


    return nblks;
}






/*
   Read number of x points in each block of a structured
   grid from a plot3D file

   Input parameters: run_mode   - either 2D/3D
                     filename   - name of the initial grid file
                     nblks      - number of blocks in structured grid
*/
long *read_npts_x_from_plot3D (char *run_mode, char *filename, long nblks)
{
    /* Return nx */
    long                    *nx;

    /* Local variables */
    FILE                    *fptr;
    long                    temp;
    int                     i;


    /* Allocate memory to array containing number
       of x points for each block */
    nx                                  = allocate_1D_long_array ("npts_x", (int)nblks);

    fptr                                = fopen (filename, "r");

    fscanf (fptr, "%ld\n", &temp);  // Read unnecessary input to temporary variable

    /* Read number of x points for each block for 2D grids */
    if (strcmp(run_mode, "2D") == 0)
    {
        for (i = 0; i < 2 * (int)nblks; i++)
        {
            fscanf (fptr, "%ld ", &temp);
            if (i%2 == 0)
                nx[i/2]                 = temp;
        }
    }

    /* Read number of x points for each block for 3D grids */
    else
    {
        for (i = 0; i < 3 * (int)nblks; i++)
        {
            fscanf (fptr, "%ld ", &temp);
            if (i%3 == 0)
                nx[i/3]                 = temp;
        }
    }

    fclose (fptr);


    return nx;
}






/*
   Read number of y points in each block of a structured
   grid from a plot3D file

   Input parameters: run_mode   - either 2D/3D
                     filename   - name of the initial grid file
                     nblks      - number of blocks in structured grid
*/
long *read_npts_y_from_plot3D (char *run_mode, char *filename, long nblks)
{
    /* Return nx */
    long                    *ny;

    /* Local variables */
    FILE                    *fptr;
    long                    temp;
    int                     i;


    /* Allocate memory to array containing number
       of y points for each block */
    ny                                  = allocate_1D_long_array ("npts_y", (int)nblks);

    fptr                                = fopen (filename, "r");

    fscanf (fptr, "%ld\n", &temp);  // Read unnecessary input to temporary variable

    /* Read number of y points for each block for 2D grids */
    if (strcmp(run_mode, "2D") == 0)
    {
        for (i = 0; i < 2 * (int)nblks; i++)
        {
            fscanf (fptr, "%ld ", &temp);
            if (i%2 == 1)
                ny[(i - 1)/2]           = temp;
        }
    }

    /* Read number of y points for each block for 3D grids */
    else
    {
        for (i = 0; i < 3 * (int)nblks; i++)
        {
            fscanf (fptr, "%ld ", &temp);
            if (i%3 == 1)
                ny[(i - 1)/3]           = temp;
        }
    }

    fclose (fptr);


    return ny;
}






/*
   Read number of z points in each block of a structured
   grid from a plot3D file

   Input parameters: filename   - name of the initial grid file
                     nblks      - number of blocks in structured grid
*/
long *read_npts_z_from_plot3D (char *filename, long nblks)
{
    /* Return nz */
    long                    *nz;

    /* Local variables */
    FILE                    *fptr;
    long                    temp;
    int                     i;


    /* Allocate memory to array containing number
       of z points for each block */
    nz                                  = allocate_1D_long_array ("npts_z", (int)nblks);

    fptr                                = fopen (filename, "r");

    fscanf (fptr, "%ld\n", &temp);  // Read unnecessary input to temporary variable

    /* Read number of z points for each block for 3D grids */
    for (i = 0; i < 3 * (int)nblks; i++)
    {
        fscanf (fptr, "%ld ", &temp);
        if (i%3 == 2)
        {
            nz[(i - 2)/3]           = temp;
        }
    }

    fclose (fptr);


    return nz;
}






/*
   Read the (x, y) coordinates along the boundaries
   which map to eta = 0 and eta = 1 from a file

   Input parameters: nx - number of boundary points
*/
point_2D **read_x_bounds_2D (int nx)
{
    /* Return bounds */
    point_2D                **bounds;

    /* Local variables */
    FILE                    *fptr;
    int                     i;


    /* Initialization */
    bounds                              = allocate_2D_point_2D_array ("x_bounds", 2, nx);

    /* Read x boundary coordinates from file */
    fptr                                = fopen ("x_boundaries.dat", "r");
    for (i = 0; i < nx; i++)
    {
        fscanf (fptr, "%LF    %LF    %LF    %LF\n",
                &bounds[0][i].x, &bounds[0][i].y, &bounds[1][i].x, &bounds[1][i].y);
    }
    fclose (fptr);


    return bounds;
}






/*
   Read the (x, y) coordinates along the boundaries
   which map to xi = 0 and xi = 1 from a file

   Input parameters: ny - number of boundary points
*/
point_2D **read_y_bounds_2D (int ny)
{
    /* Return bounds */
    point_2D                **bounds;

    /* Local variables */
    FILE                    *fptr;
    int                     i;


    /* Initialization */
    bounds                              = allocate_2D_point_2D_array ("x_bounds", 2, ny);

    /* Read x boundary coordinates from file */
    fptr                                = fopen ("y_boundaries.dat", "r");
    for (i = 0; i < ny; i++)
    {
        fscanf (fptr, "%LF    %LF    %LF    %LF\n",
                &bounds[0][i].x, &bounds[0][i].y, &bounds[1][i].x, &bounds[1][i].y);
    }
    fclose (fptr);


    return bounds;
}






/*
   Write domain boundary coordinates to a file

   Input parameters: filename   - name of the boundary coordinate file
                     n          - number of boundary points
                     bounds     - array containing boundary coordinates
*/
void write_bounds_2D (char *filename, int n, point_2D **bounds)
{
    /* Local variables */
    FILE                    *fptr;
    int                     j;


    fptr                                = fopen (filename, "w");
    for (j = 0; j < n; j++)
    {
        fprintf (fptr, "%22.16LF    %22.16LF    %22.16LF    %22.16LF\n",
                 bounds[0][j].x, bounds[0][j].y, bounds[1][j].x, bounds[1][j].y);
    }
    fclose (fptr);
}






/*
   Read the normals to the boundaries which map
   to eta = 0 and eta = 1 from a file

   Input parameters: nx - number of boundary points
*/
point_2D **read_xbound_normals_2D (int nx)
{
    /* Return nx_bounds */
    point_2D                **nx_bounds;

    /* Local variables */
    FILE                    *fptr;
    int                     i;


    /* Initialization */
    nx_bounds                           = allocate_2D_point_2D_array ("nx_bounds", 2, nx);

    /* Read x boundary coordinates from file */
    fptr                                = fopen ("xbound_normals.dat", "r");
    for (i = 0; i < nx; i++)
    {
        fscanf (fptr, "%LF    %LF    %LF    %LF\n",
                &nx_bounds[0][i].x, &nx_bounds[0][i].y, &nx_bounds[1][i].x, &nx_bounds[1][i].y);
    }
    fclose (fptr);


    return nx_bounds;
}






/*
   Read the normals to the boundaries which map
   to xi = 0 and xi = 1 from a file

   Input parameters: ny - number of boundary points
*/
point_2D **read_ybound_normals_2D (int ny)
{
    /* Return ny_bounds */
    point_2D                **ny_bounds;

    /* Local variables */
    FILE                    *fptr;
    int                     i;


    /* Initialization */
    ny_bounds                           = allocate_2D_point_2D_array ("ny_bounds", 2, ny);

    /* Read x boundary coordinates from file */
    fptr                                = fopen ("ybound_normals.dat", "r");
    for (i = 0; i < ny; i++)
    {
        fscanf (fptr, "%LF    %LF    %LF    %LF\n",
                &ny_bounds[0][i].x, &ny_bounds[0][i].y, &ny_bounds[1][i].x, &ny_bounds[1][i].y);
    }
    fclose (fptr);


    return ny_bounds;
}






/*
   Read a single block formatted 2D plot3D file
   The grid coordinates are stored in a 2D array of
   type point_2D
*/
point_2D **read_2D_singleblock_plot3D (char *filename)
{
    /* Return 2D grid */
    point_2D                **grid;

    /* Local variables */
    FILE                    *fptr;
    int                     nblks, nx, ny, i, j;


    /* Initialization */
    grid                                = NULL;

    fptr                                = fopen (filename, "r");

    /* Read number of blocks */
    fscanf (fptr, "%d\n", &nblks);

    /* Read number of points in x and y */
    fscanf (fptr, "%d    %d\n", &nx, &ny);

    /* Allocate point_2D array */
    grid                                = allocate_2D_point_2D_array ("grid", nx, ny);

    /* Read x and y coordinates */
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fscanf (fptr, "%LF ", &grid[i][j].x);
        }
    }
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fscanf (fptr, "%LF ", &grid[i][j].y);
        }
    }

    fclose (fptr);


    return grid;
}






/*
   Write a formatted 2D plot3D file with a .x extension

   Input parameters: filename   - name of the output file (should have .x extension)
                     nx         - number of x points
                     ny         - number of y points
                     pts        - 2D array of type point_2D with size (nx, ny)
                                  containing x and y coordinates
*/
void write_2D_singleblock_plot3D (char *filename, int nx, int ny, point_2D **grid)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Write grid data to a file */
    fptr                                = fopen (filename, "w");

    /* Number of blocks */
    fprintf (fptr, "1\n");

    /* Number of points in x and y directions */
    fprintf (fptr, "%d    %d\n", nx, ny);

    /* x coordinates */
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "%LF    ", grid[i][j].x);
        }
    }

    /* y coordinates */
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "%LF    ", grid[i][j].y);
        }
    }

    /* Close file */
    fclose (fptr);
}






/*
   Write a 2D grid to a XML structured grid
   file with a .vts extension

   Input parameters: filename   - name of the output file (should have .x extension)
                     nx         - number of x points
                     ny         - number of y points
                     grid       - 2D array of type point_2D with size (nx, ny)
                                  containing x and y coordinates
*/
void write_2D_singleblock_vts (char *filename, int nx, int ny, point_2D **grid, int write_data)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Write grid data to a .vts file */
    fptr                                = fopen (filename, "w");

    fprintf (fptr, "<?xml version=\"1.0\"?>\n");    // XML identifier
    fprintf (fptr, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");  // Structured grid file
    fprintf (fptr, "    <StructuredGrid WholeExtent=\"0 %d 0 %d 0 0\">\n", nx - 1, ny - 1); // Number of points in grid
    fprintf (fptr, "        <Piece Extent=\"0 %d 0 %d 0 0\">\n", nx - 1, ny - 1);   // Singleblock so piece extent = whole extent
    fprintf (fptr, "            <Points>\n");
    fprintf (fptr, "                <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

    /* Write points such that x increases quickest
       followed by y */
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF    %22.16LF    0.0\n", grid[i][j].x, grid[i][j].y);
        }
    }

    fprintf (fptr, "                </DataArray>\n");
    fprintf (fptr, "            </Points>\n");
    if (write_data == 0)
    {
        /* Close file without ending piece field to
           enable writing scalar/vector data */
        fclose (fptr);
    }
    else
    {
        fprintf (fptr, "        </Piece>\n");
        fprintf (fptr, "    </StructuredGrid>\n");
        fprintf (fptr, "</VTKFile>");

        /* Close file */
        fclose (fptr);
    }
}






/*
   Write a scalar 2D array to a .vts file

   Input parameters: filename   - name of the output file (should have .x extension)
                     nx         - number of x points
                     ny         - number of y points
                     grid       - 2D array of type point_2D with size (nx, ny)
                                  containing x and y coordinates
                     array      - scalar array to be written to file
*/
void write_long_double_to_vts (char *filename, int nx, int ny, point_2D **grid,
                               char *fieldname, long double **array)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Write grid point data to file */
    write_2D_singleblock_vts (filename, nx, ny, grid, 0);

    /* Open file to append */
    fptr                                = fopen (filename, "a");

    /* Write scalar data to file */
    fprintf (fptr, "            <PointData Scalars=\"%s\">\n", fieldname);
    fprintf (fptr, "                <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"ascii\">\n", fieldname);
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF\n", array[i][j]);
        }
    }
    fprintf (fptr, "                </DataArray>\n");
    fprintf (fptr, "            </PointData>\n");

    /* End piece field */
    fprintf (fptr, "        </Piece>\n");
    fprintf (fptr, "    </StructuredGrid>\n");
    fprintf (fptr, "</VTKFile>");

    /* Close file */
    fclose (fptr);
}






/*
   Write a point_2D 2D array to a .vts file

   Input parameters: filename   - name of the output file (should have .x extension)
                     nx         - number of x points
                     ny         - number of y points
                     grid       - 2D array of type point_2D with size (nx, ny)
                                  containing x and y coordinates
                     array      - point_2D array to be written to file
*/
void write_point_2D_to_vts (char *filename, int nx, int ny, point_2D **grid,
                            char *fieldname, point_2D **array)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Write grid point data to file */
    write_2D_singleblock_vts (filename, nx, ny, grid, 0);

    /* Open file to append */
    fptr                                = fopen (filename, "a");

    /* Write scalar data to file */
    fprintf (fptr, "            <PointData Vectors=\"%s\">\n", fieldname);
    fprintf (fptr, "                <DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", fieldname);
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF    %22.16LF    %22.16LF\n", array[i][j].x, array[i][j].y, (long double) 0.0);
        }
    }
    fprintf (fptr, "                </DataArray>\n");
    fprintf (fptr, "            </PointData>\n");

    /* End piece field */
    fprintf (fptr, "        </Piece>\n");
    fprintf (fptr, "    </StructuredGrid>\n");
    fprintf (fptr, "</VTKFile>");

    /* Close file */
    fclose (fptr);
}






/*
   Write grid derivatives for a 2D grid to a
   file with a .vts extension

   Input parameters: filename   - name of the output file (should have .x extension)
                     nx         - number of x points
                     ny         - number of y points
                     grid       - 2D array of type point_2D with size (nx, ny)
                                  containing x and y coordinates
                     dgrid      - 2D array of type grid_der_2D with size (nx, ny)
                                  containing the derivatives of x and y with respect to
                                  xi and eta
*/
void write_grid_der_2D_to_vts (char *filename, int nx, int ny, point_2D **grid,
                               grid_der_2D **dgrid)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Write grid point data to file */
    write_2D_singleblock_vts (filename, nx, ny, grid, 0);

    /* Open file to append */
    fptr                                = fopen (filename, "a");

    /* Write xi derivatives to file */
    fprintf (fptr, "            <PointData Vectors=\"xi_ders eta_ders\">\n");
    fprintf (fptr, "                <DataArray type=\"Float32\" Name=\"xi_ders\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF    %22.16LF    0.0\n", dgrid[i][j].x_der.x, dgrid[i][j].x_der.y);
        }
    }
    fprintf (fptr, "                </DataArray>\n");

    /* Write eta derivatives to file */
    fprintf (fptr, "                <DataArray type=\"Float32\" Name=\"eta_ders\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF    %22.16LF    0.0\n", dgrid[i][j].y_der.x, dgrid[i][j].y_der.y);
        }
    }
    fprintf (fptr, "                </DataArray>\n");
    fprintf (fptr, "            </PointData>\n");

    /* End piece field */
    fprintf (fptr, "        </Piece>\n");
    fprintf (fptr, "    </StructuredGrid>\n");
    fprintf (fptr, "</VTKFile>");

    /* Close file */
    fclose (fptr);
}






/*
   Write grid derivatives for a 2D grid to a
   file with a .vts extension

   Input parameters: filename   - name of the output file (should have .x extension)
                     nx         - number of x points
                     ny         - number of y points
                     grid       - 2D array of type point_2D with size (nx, ny)
                                  containing x and y coordinates
                     dgrid      - 2D array of type grid_der_2D with size (nx, ny)
                                  containing the derivatives of x and y with respect to
                                  xi and eta
*/
void write_grid_dder_2D_to_vts (char *filename, int nx, int ny, point_2D **grid,
                                grid_dder_2D **d2grid)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Write grid point data to file */
    write_2D_singleblock_vts (filename, nx, ny, grid, 0);

    /* Open file to append */
    fptr                                = fopen (filename, "a");

    /* Write xi 2nd derivatives to file */
    fprintf (fptr, "            <PointData Vectors=\"ddxi2 ddxideta ddeta2\">\n");
    fprintf (fptr, "                <DataArray type=\"Float32\" Name=\"ddxi2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF    %22.16LF    0.0\n", d2grid[i][j].x_der.x, d2grid[i][j].x_der.y);
        }
    }
    fprintf (fptr, "                </DataArray>\n");

    /* Write eta 2nd derivatives to file */
    fprintf (fptr, "                <DataArray type=\"Float32\" Name=\"ddeta2\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF    %22.16LF    0.0\n", d2grid[i][j].y_der.x, d2grid[i][j].y_der.y);
        }
    }
    fprintf (fptr, "                </DataArray>\n");

    /* Write mixed derivatives to file */
    fprintf (fptr, "                <DataArray type=\"Float32\" Name=\"dxideta\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (j = 0; j < ny; j++)
    {
        for (i = 0; i < nx; i++)
        {
            fprintf (fptr, "                    %22.16LF    %22.16LF    0.0\n", d2grid[i][j].xy_der.x, d2grid[i][j].xy_der.y);
        }
    }
    fprintf (fptr, "                </DataArray>\n");
    fprintf (fptr, "            </PointData>\n");

    /* End piece field */
    fprintf (fptr, "        </Piece>\n");
    fprintf (fptr, "    </StructuredGrid>\n");
    fprintf (fptr, "</VTKFile>");

    /* Close file */
    fclose (fptr);
}






/********************************************************************/
/*                                                                  */
/*                     ARRAY WRITING FUNCTIONS                      */
/*                                                                  */
/********************************************************************/



/*
   Write a 1D long double array to .dat file
   such that each row is on a separate line

   Input parameters: filename   - name of output file
                     n          - size of the array
                     array      - array to be written to file
*/
void write_long_double_1D (char *filename, int n, long double *array)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i;


    fptr                                = fopen (filename, "w");
    for (i = 0; i < n; i++)
    {
        fprintf (fptr, "%22.16LF\n", array[i]);
    }
    fclose (fptr);
}



/*
   Write a 2D long double array to .dat file
   such that each row is on a separate line

   Input parameters: filename   - name of output file
                     nx         - size of 1st array dimension
                     ny         - size of 2nd array dimension
                     array      - array to be written to file
*/
void write_long_double_2D (char *filename, int nx, int ny, long double **array)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    fptr                                = fopen (filename, "w");
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", array[i][j]);
        }
        fprintf (fptr, "\n");
    }
    fclose (fptr);
}



/*
   Write a 1D long double array to .dat file
   such that each row is on a separate line

   Input parameters: filename   - name of output file
                     n          - size of the array
                     array      - array to be written to file
*/
void write_1D_point_2D (char *filename, int n, point_2D *array)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i;


    fptr                                = fopen (filename, "w");
    for (i = 0; i < n; i++)
    {
        fprintf (fptr, "%22.16LF    %22.16LF\n", array[i].x, array[i].y);
    }
    fclose (fptr);
}



/*
   Write x and y components of a 2D point_2D array
   to separate .dat files such that each row is on
   a separate line

   Input parameters: filename1  - name of x component output file
                     filename2  - name of y component output file
                     nx         - size of 1st array dimension
                     ny         - size of 2nd array dimension
                     array      - array to be written to file
*/
void write_2D_point_2D (char *filename1, char *filename2, int nx, int ny, point_2D **array)
{
    /* Local variables */
    FILE                    *fptr1, *fptr2;
    int                     i, j;


    fptr1                               = fopen (filename1, "w");
    fptr2                               = fopen (filename2, "w");
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr1, "%22.16LF    ", array[i][j].x);
            fprintf (fptr2, "%22.16LF    ", array[i][j].y);
        }
        fprintf (fptr1, "\n");
        fprintf (fptr2, "\n");
    }
    fclose (fptr1);
    fclose (fptr2);
}



/*
   Write derivatives of all (x, y) points in a
   2D grid to separate files (derivatives with
   respect to xi and eta)

   Input parameters: filename1  - name of output file containing x component of x_der
                     filename2  - name of output file containing x component of y_der
                     filename3  - name of output file containing y component of x_der
                     filename4  - name of output file containing y component of y_der
                     nx         - size of 1st array dimension
                     ny         - size of 2nd array dimension
                     array      - array to be written to file
*/
void write_grid_der_2D (char *filename1, char *filename2, char *filename3,
                        char *filename4, int nx, int ny, grid_der_2D **array)
{
    /* Local variables */
    FILE                    *fptr1, *fptr2, *fptr3, *fptr4;
    int                     i, j;


    fptr1                               = fopen (filename1, "w");
    fptr2                               = fopen (filename2, "w");
    fptr3                               = fopen (filename3, "w");
    fptr4                               = fopen (filename4, "w");
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr1, "%22.16LF    ", array[i][j].x_der.x);
            fprintf (fptr2, "%22.16LF    ", array[i][j].y_der.x);
            fprintf (fptr3, "%22.16LF    ", array[i][j].x_der.y);
            fprintf (fptr4, "%22.16LF    ", array[i][j].y_der.y);
        }
        fprintf (fptr1, "\n");
        fprintf (fptr2, "\n");
        fprintf (fptr3, "\n");
        fprintf (fptr4, "\n");
    }
    fclose (fptr1);
    fclose (fptr2);
    fclose (fptr3);
    fclose (fptr4);
}



/*
   Write second derivatives of all (x, y) points
   in a 2D grid to separate files (derivatives with
   respect to xi and eta)

   Input parameters: nx     - size of 1st array dimension
                     ny     - size of 2nd array dimension
                     array  - array to be written to file
*/
void write_grid_dder_2D (int nx, int ny, grid_dder_2D **array)
{
    /* Local variables */
    FILE                    *fptr1, *fptr2, *fptr3, *fptr4, *fptr5, *fptr6;
    int                     i, j;


    fptr1                               = fopen ("d2xdxi2.dat", "w");
    fptr2                               = fopen ("d2ydxi2.dat", "w");
    fptr3                               = fopen ("d2xdxideta.dat", "w");
    fptr4                               = fopen ("d2ydxideta.dat", "w");
    fptr5                               = fopen ("d2xdeta2.dat", "w");
    fptr6                               = fopen ("d2ydeta2.dat", "w");
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr1, "%22.16LF    ", array[i][j].x_der.x);
            fprintf (fptr2, "%22.16LF    ", array[i][j].x_der.y);
            fprintf (fptr3, "%22.16LF    ", array[i][j].xy_der.x);
            fprintf (fptr4, "%22.16LF    ", array[i][j].xy_der.y);
            fprintf (fptr5, "%22.16LF    ", array[i][j].y_der.x);
            fprintf (fptr6, "%22.16LF    ", array[i][j].y_der.y);
        }
        fprintf (fptr1, "\n");
        fprintf (fptr2, "\n");
        fprintf (fptr3, "\n");
        fprintf (fptr4, "\n");
        fprintf (fptr5, "\n");
        fprintf (fptr6, "\n");
    }
    fclose (fptr1);
    fclose (fptr2);
    fclose (fptr3);
    fclose (fptr4);
    fclose (fptr5);
    fclose (fptr6);
}



/*
   Write the values of (x, y, p, q) for a
   given grid to a single file. Each variable
   is written as a matrix of size (nx, ny)

   Input parameters: filename   - name of the output file
                     nx         - number of x points
                     ny         - number of y points
                     grid       - array containing (x, y) values
                     pq         - array containing (p, q) values
*/
void write_x_to_file_2D (char *filename, int nx, int ny, point_2D **grid, point_2D **pq)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Open file */
    fptr                                = fopen (filename, "w");

    /* Write x values to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", grid[i][j].x);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write y values to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", grid[i][j].y);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write p values to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", pq[i][j].x);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write q values to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", pq[i][j].y);
        }
        fprintf (fptr, "\n");
    }

    /* Close file */
    fclose(fptr);
}



/*
   Write all derivatives dx/dxi, dy/dxi,
   dx/deta and dy/deta to a single file.
   Each derivative is written to the file
   as a matrix of size (nx, ny)

   Input parameters: filename   - name of the output file
                     nx         - number of x points
                     ny         - number of y points
                     dgrid      - array of type grid_der_2D containing
                                  all derivatives
*/
void write_dgrid_to_file_2D (char *filename, int nx, int ny, grid_der_2D **dgrid)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Open file */
    fptr                                = fopen (filename, "w");

    /* Write dx/dxi to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", dgrid[i][j].x_der.x);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write dy/dxi to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", dgrid[i][j].x_der.y);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write dx/deta to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", dgrid[i][j].y_der.x);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write dy/deta to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", dgrid[i][j].y_der.y);
        }
        fprintf (fptr, "\n");
    }

    /* Close file */
    fclose (fptr);
}



/*
   Write all first derivative coefficients
   to a single file. Each coefficient is
   written to the file as a matrix of size (nx, ny)

   Input parameters: filename   - name of the output file
                     nx         - number of x points
                     ny         - number of y points
                     coeffs     - array of type coeffs_1 containing
                                  all coefficients
*/
void write_coeffs_to_file_2D (char *filename, int nx, int ny, coeffs_1 **coeffs)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Open file */
    fptr                                = fopen (filename, "w");

    /* Write alpha to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", coeffs[i][j].alpha);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write beta to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", coeffs[i][j].beta);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write gamma to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", coeffs[i][j].gamma);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write J to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", coeffs[i][j].J);
        }
        fprintf (fptr, "\n");
    }

    /* Close file */
    fclose (fptr);
}



/*
   Write all second order derivatives to
   a single file. Each derivative is written
   to the file as a matrix of size (nx, ny)

   Input parameters: filename   - name of the output file
                     nx         - number of x points
                     ny         - number of y points
                     d2grid     - array of type grid_dder_2D containing
                                  all derivatives
*/
void write_d2grid_to_file_2D (char *filename, int nx, int ny, grid_dder_2D **d2grid)
{
    /* Local variables */
    FILE                    *fptr;
    int                     i, j;


    /* Open file */
    fptr                                = fopen (filename, "w");

    /* Write d2x/dxi2 to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", d2grid[i][j].x_der.x);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write d2y/dxi2 to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", d2grid[i][j].x_der.y);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write d2x/dxideta to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", d2grid[i][j].xy_der.x);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write d2y/dxideta to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", d2grid[i][j].xy_der.y);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write d2x/deta2 to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", d2grid[i][j].y_der.x);
        }
        fprintf (fptr, "\n");
    }
    fprintf (fptr, "\n");

    /* Write d2y/deta2 to file */
    for (i = 0; i < nx; i++)
    {
        for (j = 0; j < ny; j++)
        {
            fprintf (fptr, "%22.16LF    ", d2grid[i][j].y_der.y);
        }
        fprintf (fptr, "\n");
    }

    /* Close file */
    fclose (fptr);
}






/*
   Write elements of a 1D newton_sol_2D
   array to a file

   Input parameters: filename - name of output file
                     n        - number of array elements
                     array    - the array to be written out
*/
void write_1D_bh_sol_2D (char *filename, int nx, int ny, bh_sol_2D *array)
{
    /* Open file */
    FILE *fptr                          = fopen (filename, "w");

    /* Write array elements to a file */
    for (int i = 0; i < nx * ny; i++)
    {
        fprintf (fptr, "%d    %d    %10.6LF    %10.6LF    %10.6LF    %10.6LF\n",
                 i%nx, i/nx, array[i].x, array[i].y, array[i].p, array[i].q);
    }

    /* Close file */
    fclose (fptr);
}
