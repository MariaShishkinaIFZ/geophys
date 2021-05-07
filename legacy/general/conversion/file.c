#include "file.h"

FILE *f_open(const char *path, const char *name, const char *mode)
{
  char path_i[50];
  FILE *fi;
  strcpy(path_i,path);
  strcat(path_i,name);
  fi=fopen(path_i,mode);
  if (fi==NULL)
  {
    printf("Incorrect path or filename! Program terminated."); exit(0);
  }
  return fi;
}

void putintxy(int x, int y, int n, int a)
{
  int i;
  for (i=0;i<n;i++)
    printf(" ");
  printf("%*d",n,a);
}

