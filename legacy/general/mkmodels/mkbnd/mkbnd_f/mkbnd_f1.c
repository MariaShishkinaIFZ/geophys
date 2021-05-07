#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "submem.h"
#include "subio.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/


int main(int narg, char **argv)
{
  int Nx, Ny, Nz, i, j, k;
  int ask;
  float x0, y0, z0, h;
  float a1, a2, a3;
  float nx, ny, nz;
  float b1, b2, b3;
  float c1, c2, c3;
  float **Int;
  char out_file[FLNLN];
  FILE *fp;
 
 if (narg < 2) errquit("Insufficient parameters.\n");
 
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%d%d%d",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
 
  if (strcmp(argv[1],"-f")==0) /*Use parameters file with the name given as 2nd parameter*/
  {
    if (narg < 3) errquit("Parameters file name not given.\n");
    fp=file_open(argv[2],"rt");
    ask=0;
  }
  else
  {
    if (strcmp(argv[1],"-s")==0) /*Read parameters from standart input*/
    {
      fp=stdin;
      ask=1;
    }
    else
    {
      printf("Unrecognized flag %s.\n",argv[1]);
      exit(1);
    }
  }
  
  if (ask) fprintf(stdout,"Input x y z, - coordinates of 1st point \n");
  fscanf(fp,"%g%g%g",&a1,&a2,&a3);
  if (ask) fprintf(stdout,"Input x y z - coordinates of 2nd point \n");
  fscanf(fp,"%g%g%g",&b1, &b2, &b3);
  if (ask) fprintf(stdout,"Input x y z - coordinates of 3d point \n");
  fscanf(fp,"%g%g%g",&c1, &c2, %c3);
  if (ask) fprintf(stdout,"Input boundary model filename.2d (%d chars max) \n",FLNLN);
  fscanf(fp,"%s",out_file);
  if (!ask) fclose(fp);
  
  nx=(b2-a2)*(c3-a3)-(b3-a3)*(b2-a2);
  
  ny=(b1-a1)*(c3-a3)-(b3-a3)*(c1-a1);
  
  nz=(b1-a1)*(c2-c3-c1)-(b2-b1)*(c1-a1);
 

  alloc_float_matrix(&Int,Ny,Nx);
 
  for (i=0; i<Ny; i++){
    for (j=0; j<Nx; j++){ 
	Int[i][j]=(nz*a3-nx*(j*h-a1)-ny*(i*h-a2))/nz;
    }
  }

  write_2D_surface(Int,Ny,Nx,out_file);

  free_float_matrix(Int,Ny,Nx);
}
