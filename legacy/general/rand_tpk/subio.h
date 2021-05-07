#include <stdio.h>
/*Header file for io subroutines*/

FILE *file_open(const char *fname, const char *mode);
int file_close(FILE *f);
void errquit(const char *message);
