#include <cstdio>
#include <cmath>
#include "subio.h"

void read_3D_volume(float ***Buf, long Nz, long Ny, long Nx, const char *filename)
{
  long j,k;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fread(Buf[k][j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_3D_volume(float ***Buf, long Nz, long Ny, long Nx, const char *dir, const char *filename)
{
  long j,k;
  FILE *fo;

  fo=dir_file_open(dir,filename,"wb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fwrite(Buf[k][j],sizeof(float),Nx,fo);
  fclose(fo);
}

void read_2D_surface(float **Buf, long Ny, long Nx, const char *filename)
{
  long j;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (j=0;j<Ny;j++)
    fread(Buf[j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_2D_surface(float **Buf, long Ny, long Nx, const char *dir, const char *filename)
{
  long j;
  FILE *fo;

  fo=dir_file_open(dir,filename,"wb");
  for (j=0;j<Ny;j++)
    fwrite(Buf[j],sizeof(float),Nx,fo);
  fclose(fo);
}
