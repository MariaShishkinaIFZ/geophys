#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.coarse"
#define FLNLN 80 	/*Maximum length of the file name*/

int main(int narg, char **argv)
{
  long int Nx, Ny, Nz;
  int i, j, k;
  int ask;
  float rc_par;
  float ***Vu, ***Vl;
  char u_file[FLNLN], l_file[FLNLN];
  char out_file[FLNLN];
  FILE *fp;
  
  printf("This program is desighned to convert 3D corse grid file and ray coverage to text table.\n");
  if (narg < 2) errquit("Insufficient parameters.\n");
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
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
  if (ask) fprintf(stdout,"Input ray coverage filename (%d chars max) > ",FLNLN);
  fscanf(fp,"%s",l_file);  
  if (ask) fprintf(stdout,"Input ray covearge cutoff parametr > ");
  fscanf(fp,"%g",&rc_par);
  if (ask) fprintf(stdout,"Input outpu table filename (%d chars max) > ",FLNLN);
  fscanf(fp,"%s",out_file);  
  if (!ask) fclose(fp);
 
  
  alloc_float_3d(&Vu,Nz,Ny,Nx);
  alloc_float_3d(&Vl,Nz,Ny,Nx);
  
  fp=file_open(out_file,"wt");
  
  read_3D_volume(Vu,Nz,Ny,Nx,u_file);
  read_3D_volume(Vl,Nz,Ny,Nx,l_file);
 
  
  for (i=0;i<Ny;i++)
    for (j=0;j<Nx;j++)
      for (k=0;k<Nz;k++)
      {
	if (Vl[k][i][j]>rc_par)	  fprintf(fp,"%10g %10g \n",Vl[k][i][j],fabs(Vu[k][i][j]));  
      }
  
  fclose(fp);
  free_float_3d(Vl,Nz,Ny,Nx);
  free_float_3d(Vu,Nz,Ny,Nx);
  exit(1);
}
