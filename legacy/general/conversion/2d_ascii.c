#include <stdio.h>

#include "submem.h"
#include "subio.h"

int main(int narg, char **argv)
{
  int Nx, Ny;   /*Dimentions of grid*/
  int i, j;
  float **M;  /*Data matrix*/

  FILE *fi, *fo;

  if (narg < 5) errquit("Insufficient parameters.\n USAGE: 2d_ascii <file[.2d]> <out ascii file> Nx Ny\n");

  sscanf(argv[3],"%d",&Nx);
  sscanf(argv[4],"%d",&Ny);

  alloc_float_matrix(&M,Ny,Nx);

  fi=file_open(argv[1],"rb");
  for (i=0;i<Ny;i++)
    fread(M[i],sizeof(float),Nx,fi);
  fclose(fi);

  fo=file_open(argv[2],"wt");
  for (i=0;i<Ny;i++)
  {
    for (j=0;j<Nx;j++)
      fprintf(fo,"%10g ",M[i][j]);
    fprintf(fo,"\n");
  }
  fclose(fo);

  free_float_matrix(M,Ny,Nx);
  printf("Success.\n");
}
