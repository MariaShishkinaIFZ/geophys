#include <stdio.h>
#include "subio.h"
#include "iogras.h"

void read_3D_volume(float ***Buf, long Nz, long Ny, long Nx, char *filename)
{
  long j,k;
  FILE *fi;
  
  fi=file_open(filename,"rb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fread(Buf[k][j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_3D_volume(float ***Buf, long Nz, long Ny, long Nx, char *filename)
{
  long j,k;
  FILE *fo;
  
  fo=file_open(filename,"wb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fwrite(Buf[k][j],sizeof(float),Nx,fo);
  fclose(fo);
}

void read_2D_surface(float **Buf, long Ny, long Nx, char *filename)
{
  long j;
  FILE *fi;
  
  fi=file_open(filename,"rb");
  for (j=0;j<Ny;j++)
    fread(Buf[j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_2D_surface(float **Buf, long Ny, long Nx, char *filename)
{
  long j;
  FILE *fo;
  
  fo=file_open(filename,"wb");
  for (j=0;j<Ny;j++)
    fwrite(Buf[j],sizeof(float),Nx,fo);
  fclose(fo);
}
