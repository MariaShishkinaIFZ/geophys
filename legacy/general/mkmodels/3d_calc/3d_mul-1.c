#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 80 	/*Maximum length of the file name*/

int main(int narg, char **argv)
{
  long int Nx, Ny, Nz;
  int i, j, k;
  int ask;
  float x0, y0, z0, h;
  float ***Vu, ***Vl;
  char u_file[FLNLN], l_file[FLNLN];
  char out_file[FLNLN];
  FILE *fp;
  
  if (narg < 2) errquit("Insufficient parameters.\n");
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
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
      exit(0);
    }
  }
  
  if (ask) fprintf(stdout,"Input first velocity model filename (%d chars max) > ",FLNLN);
  fscanf(fp,"%s",u_file);
  if (ask) fprintf(stdout,"Input second velocity model filename (%d chars max) > ",FLNLN);
  fscanf(fp,"%s",l_file);  
  if (ask) fprintf(stdout,"Input resulting velocity model filename (%d chars max) > ",FLNLN);
  fscanf(fp,"%s",out_file);
  if (!ask) fclose(fp);
 
  alloc_float_3d(&Vu,Nz,Ny,Nx);
  alloc_float_3d(&Vl,Nz,Ny,Nx);
  
  read_3D_volume(Vu,Nz,Ny,Nx,u_file);
  read_3D_volume(Vl,Nz,Ny,Nx,l_file);
  
  for (i=0;i<Ny;i++)
    for (j=0;j<Nx;j++)
      for (k=0;k<Nz;k++)
        Vu[k][i][j]*=Vl[k][i][j];
      
  write_3D_volume(Vu,Nz,Ny,Nx,out_file);
  
  free_float_3d(Vl,Nz,Ny,Nx);
  free_float_3d(Vu,Nz,Ny,Nx);
  exit(1);
}
