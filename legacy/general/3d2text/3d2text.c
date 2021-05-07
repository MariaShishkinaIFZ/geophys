#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/

int main(int narg, char **argv)
{
  long int Nx, Ny, Nz;
  int i, j, k;
  float x0, y0, z0, h;
  float ***Vu;
  char u_file[FLNLN];
  FILE *fp;
  
  if (narg < 2) errquit("Use: 3d2text <velocyty_model.3d>.\n");
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
 
  
  if (narg < 2) errquit("Use: 3d2text <velocyty_model.3d>.\n");
  
//   fscanf(argv[1],"%s",u_file);
 
  alloc_float_3d(&Vu,Nz,Ny,Nx);
  read_3D_volume(Vu,Nz,Ny,Nx,argv[1]);
  
  for (i=0;i<Ny;i++)
    for (j=0;j<Nx;j++)
      for (k=0;k<Nz;k++)
        fprintf(stdout,"%10g %10g %10g %15g \n", x0+j*h, y0+i*h, z0+k*h, Vu[k][i][j]);
      
  
  free_float_3d(Vu,Nz,Ny,Nx);
  exit(1);
}
