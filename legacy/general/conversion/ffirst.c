#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "subio.h"
#include "iogras.h"
#include "submem.h"

#define GEOMFILE "model-geometry.gras"

#define DUMMY_LIM 1e9
#define DUMMY_VALUE 1e10

int main(int narg, char **argv)
{
  long Nx, Ny, Nz;
  long i,j,k;
  float x0,y0,z0;
  float h, V;
  float ***data;
  float **res;
  FILE *fp;
  
  if (narg<4)
  {
    fprintf(stderr,"Insufficient parameters. Use:\n3lin <data file [.3d]> <out file [.2d]> V\n");
    exit(EXIT_FAILURE);
  }
  
  sscanf(argv[3],"%g",&V);
  
  fp = file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
  
  alloc_float_3d(&data,Nz,Ny,Nx);
  alloc_float_matrix(&res,Ny,Nx);
  
  read_3D_volume(data,Nz,Ny,Nx,argv[1]);
  
  for (j=0;j<Ny;j++)
    for (i=0;i<Nx;i++)
    {
      for (k=1;k<Nz;k++)
      {
	if ((data[k-1][j][i]<V) && (data[k][j][i]>=V))
	{
	  res[j][i]=z0+h*k+(V-data[k-1][j][i])/(data[k][j][i]-data[k-1][j][i])*h;
	  break;
	}
      }
      if (k==Nz) res[j][i]=DUMMY_VALUE;
    }

  write_2D_surface(res,Ny,Nx,argv[2]);
                 
  free_float_matrix(res,Ny,Nx);
  free_float_3d(data,Nz,Ny,Nx);
}
