#include <stdio.h>
#include <stdlib.h>
#include "subio.h"

FILE *file_open(const char *fname, const char *mode)
{
  FILE *f;
  char s[80];
  if (!(f=fopen(fname,mode)))
  {
	 sprintf(s,"Can't open file %s",fname);
	 perror(s);
	 printf("Program terminated.\n");
	 exit(EXIT_FAILURE);
  }
  return f;
}

int file_close(FILE *f)
{
  return fclose(f);
}

void errquit(const char *s)
{
  fputs(s,stderr); exit(EXIT_FAILURE);
}
