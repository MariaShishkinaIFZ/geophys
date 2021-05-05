#include <cstdio>
#include <cstdlib>
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

FILE *dir_file_open(const char *dir, const char *fname, const char *mode)
{
  FILE *f;
  char s[256], full_name[256];
  sprintf(full_name,"%s/%s",dir,fname);
  if (!(f=fopen(full_name,mode)))
  {
	 sprintf(s,"Can't open file %s",full_name);
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
  puts(s); exit(EXIT_FAILURE);
}
