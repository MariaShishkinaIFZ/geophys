#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 80 	/*Maximum length of the file name*/

int main(int narg, char **argv)
{
  int Nx, Ny, Nz, i, j, k;
  float x0, y0, z0, h;
  float x, y, z;
  float ***Vu, ***Vl, **Int;
  char u_file[FLNLN], l_file[FLNLN], int_file[FLNLN];
  char out_file[FLNLN];
  FILE *fp;
  
  if (narg<2)
  {
    fprintf(stderr,"Use: stacklays <parameters file>\n");
    exit(0);
  }
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
 
  fp=file_open(argv[1],"rt");
  fscanf(fp,"%s",u_file);
  fscanf(fp,"%s",l_file);
  fscanf(fp,"%s",int_file);
  fscanf(fp,"%s",out_file);
  fclose(fp);
 
  alloc_float_3d(&Vu,Nz,Ny,Nx);
  alloc_float_3d(&Vl,Nz,Ny,Nx);
  alloc_float_matrix(&Int,Ny,Nx);
  
  read_3D_volume(Vu,Nz,Ny,Nx,u_file);
  read_3D_volume(Vl,Nz,Ny,Nx,l_file);
  read_2D_surface(Int,Ny,Nx,int_file);
  
  for (i=0;i<Ny;i++)
    for (j=0;j<Nx;j++)
      for (k=(int)floor((Int[i][j]-z0)/h);k<Nz;k++)
        Vu[k][i][j]=Vl[k][i][j];
      
  write_3D_volume(Vu,Nz,Ny,Nx,out_file);
  
  free_float_matrix(Int,Ny,Nx);
  free_float_3d(Vl,Nz,Ny,Nx);
  free_float_3d(Vu,Nz,Ny,Nx);
}
