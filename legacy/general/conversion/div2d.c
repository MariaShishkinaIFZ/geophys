#include <stdio.h>

#include "subio.h"
#include "submem.h"

int main(int narg, char **argv)
{
  int Nx, Ny;
  int i, j;
  float **D1, **D2;
  FILE *fi, *fo;

  if (narg < 6) errquit("Insufficient parameters.\n USAGE: div2d <file 1[.2d]> <file 2[.2d]> <out file[.2d]> Nx Ny\n");

  sscanf(argv[4],"%d",&Nx);
  sscanf(argv[5],"%d",&Ny);

  alloc_float_matrix(&D1,Ny,Nx);
  alloc_float_matrix(&D2,Ny,Nx);

  fi=file_open(argv[1],"rb");
  for (i=0;i<Ny;i++)
    fread(D1[i],sizeof(float),Nx,fi);
  fclose(fi);
  fi=file_open(argv[2],"rb");
  for (i=0;i<Ny;i++)
    fread(D2[i],sizeof(float),Nx,fi);
  fclose(fi);

  fo=file_open(argv[3],"wb");
  for (i=0;i<Ny;i++)
  {
    for (j=0;j<Nx;j++) D1[i][j]/=D2[i][j];
    fwrite(D1[i],sizeof(float),Nx,fo);
  }
  fclose(fo);

  free_float_matrix(D2,Ny,Nx);
  free_float_matrix(D1,Ny,Nx);
  printf("Success.\n");
  return 0;
}
