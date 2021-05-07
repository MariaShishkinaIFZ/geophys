#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define DZFILE "zgrad.3d"
#define DZMAXFILE "zgradmax.3d"

#define DUMMY_LIM 1e9
#define DUMMY_VALUE 1e10

int main(int narg, char **argv)
{
  long Nx, Ny, Nz;
  long i,j,k;
  float x0,y0,z0;
  float h;
  float ***data;
  float ***res;
  FILE *fp;
  
  if (narg<2)
  {
    fprintf(stderr,"Insufficient parameters. Use:\n3lin <data file [.3d]>\n");
    exit(EXIT_FAILURE);
  }
  
  fp = file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
  
  alloc_float_3d(&data,Nz,Ny,Nx);
  alloc_float_3d(&res,Nz,Ny,Nx);
  
  read_3D_volume(data,Nz,Ny,Nx,argv[1]);
  
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      for (i=0;i<Nx;i++)
      {
        if ((k>0) && (k<Nz-1))
	{
	  if ((data[k+1][j][i]<DUMMY_LIM) && (data[k-1][j][i]<DUMMY_LIM))
	    res[k][j][i]=(data[k+1][j][i]-data[k-1][j][i])/(2.0*h);
	}
	else
	  res[k][j][i]=DUMMY_VALUE;
      }

  write_3D_volume(res,Nz,Ny,Nx,DZFILE);
  
  for (k=2;k<Nz-2;k++)
    for (j=0;j<Ny;j++)
      for (i=0;i<Nx;i++)
      {
        if ((res[k][j][i]>res[k-1][j][i]) && (res[k][j][i]>res[k+1][j][i]))
	  data[k][j][i]=res[k][j][i];
	else
	  data[k][j][i]=DUMMY_VALUE;
      }
      
  for (j=0;j<Ny;j++)
    for (i=0;i<Nx;i++)
      data[0][j][i]=data[1][j][i]=data[Nz-1][j][i]=data[Nz-2][j][i]=DUMMY_VALUE;
  
  write_3D_volume(data,Nz,Ny,Nx,DZMAXFILE);
          
  free_float_3d(res,Nz,Ny,Nx);
  free_float_3d(data,Nz,Ny,Nx);
}
