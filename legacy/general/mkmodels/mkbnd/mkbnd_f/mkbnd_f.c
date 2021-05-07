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
  float alfa, delta;
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
  
  if (ask) fprintf(stdout,"Input x y z, - coordinates \n");
  fscanf(fp,"%g%g%g",&a1,&a2,&a3);
  if (ask) fprintf(stdout,"Input hade and strech angles \n");
  fscanf(fp,"%g%g",&alfa, &delta);
  if (ask) fprintf(stdout,"Input boundary model filename.2d (%d chars max) \n",FLNLN);
  fscanf(fp,"%s",out_file);
  if (!ask) fclose(fp);
 
  if (delta>=90)
  {
    fprintf(stderr,"Hade should be  in interval between 0 and 90 degrees!\n");
    exit(0);
  }

  if (delta<0)
  {
    fprintf(stderr,"Hade should be  in interval between 0 and 90 degrees!\n");
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
