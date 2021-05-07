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
  
  if (narg < 2) errquit("Use: text23d <velocyty_model.txt> <velocyty_model.3d>\n");
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
 
  // According to general agreements:
  //    * Z axes has 0 on the ground surface, and direction toward Earth center
  alloc_float_3d(&Vu,Nz,Ny,Nx);
 
  fp=file_open(argv[1],"rt");
  
  for (i=0;i<Ny;i++)
    for (j=0;j<Nx;j++)
      for (k=0;k<Nz;k++) {
        int x,y,z;
        fscanf(fp,"%g %g %g %g \n", &x0, &y0,&z0,&Vu[k][i][j]);
      }


  write_3D_volume(Vu,Nz,Ny,Nx,argv[2]);
      
  
  free_float_3d(Vu,Nz,Ny,Nx);
  exit(1);
}
