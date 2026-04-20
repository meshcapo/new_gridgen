#ifndef ERRORS_H_
#define ERRORS_H_



/*
   Function to print an error message and
   exit with failure
*/
void fatal_error (char* );

/*
   Print a message to screen if dynamic memory
   allocation for a array fails
*/
void array_allocation_msg (char *, int);

void file_extension_error (char *);


#endif // ERRORS_H_
