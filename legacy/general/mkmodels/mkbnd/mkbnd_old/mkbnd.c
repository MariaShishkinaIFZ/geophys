#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iogras.h"
#include "submem.h"
#include "subio.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 80 	/*Maximum length of the file name*/


int main(int narg, char **argv)
{
  int Nx, Ny, Nz, i, j, k;
  float x0, y0, z0, h;
  float a1, a2, a3;
  float nx, ny, nz;
  float alfa, delta;
  float **Int;
  char out_file[FLNLN];
  FILE *fp;
 
  if (narg<2)
  {
    fprintf(stderr,"Use: mkbnd <parameters file>\n");
    exit(0);
  }
 
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
 
  fp=file_open(argv[1],"rt");
  fscanf(fp,"%g%g%g",&a1,&a2,&a3);
  fscanf(fp,"%g%g",&alfa, &delta);
  fscanf(fp,"%s",out_file);
  fclose(fp);

  if (delta>=90)
  {
    fprintf(stderr,"Hade should be  in interval between 0 and 90 degrees!!!\n");
    exit(0);
  }

  if (delta<0)
  {
    fprintf(stderr,"Hade should be  in interval between 0 and 90 degrees!!!\n");
    exit(0);
  }

  alfa =(alfa*(M_PI/180.0));
  delta=(delta*(M_PI/180.0));

  nx=sin(delta)*sin(alfa);
  ny=sin(delta)*cos(alfa);
  nz=-cos(delta);

  alloc_float_matrix(&Int,Ny,Nx);
 
  for (i=0; i<Ny; i++){
    for (j=0; j<Nx; j++){ 
	Int[i][j]=(nz*a3-nx*(j*h-a1)-ny*(i*h-a2))/nz;
    }
  }

  write_2D_surface(Int,Ny,Nx,out_file);

  free_float_matrix(Int,Ny,Nx);
}
