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
  long Nx, Ny, Nz, i, j, k;
  int ask;
  float x0, y0, z0, h;
  float ***Vu, ***Vl, **Int;
  char u_file[FLNLN], l_file[FLNLN], int_file[FLNLN];
  char out_file[FLNLN];
  FILE *fp;
  
  
  if (narg < 2) errquit("Insufficient parameters.\n");
  
//  fprintf(stdout,"chekpoint (%d chars max) > \n",FLNLN);
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
 
//  fprintf(stderr,"long: %lu int: %lu\n",sizeof(long int),sizeof(int));
 
//  fprintf(stderr,"%ld %ld %ld\n",Nz,Ny,Nx);
 
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
  
//  fprintf(stdout,"chekpoint (%d chars max) > \n",FLNLN);
  
  if (ask) fprintf(stdout,"Input upper velocity model filename (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",u_file);
  if (ask) fprintf(stdout,"Input lower velocity model filename (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",l_file);  
  if (ask) fprintf(stdout,"Input boundary filename (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",int_file);
  if (ask) fprintf(stdout,"Input resulting velocity model filename (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",out_file);
  if (!ask) fclose(fp);
 
  printf("check 1\n");
  
  alloc_float_3d(&Vu,Nz,Ny,Nx);
  alloc_float_3d(&Vl,Nz,Ny,Nx);
  alloc_float_matrix(&Int,Nz,Nx); //Ny-> Nz
  
  read_3D_volume(Vu,Nz,Ny,Nx,u_file);
  read_3D_volume(Vl,Nz,Ny,Nx,l_file);
  read_2D_surface(Int,Nz,Nx,int_file); //Ny-> Nz
  
    printf("check 1\n");
  
  for (i=0;i<Nz;i++) //Ny-> Nz
    for (j=0;j<Nx;j++)
      for (k=(int)ceil((Int[i][j]-y0)/h);k<Ny;k++) //Nz-> Ny, z0-> y0
        Vl[i][k][j]=Vu[i][k][j];
      
//  fprintf(stderr,"%ld %ld %ld\n",Nz,Ny,Nx);
	
	  printf("check 1\n");
  
  write_3D_volume(Vl,Nz,Ny,Nx,out_file);
  
  free_float_matrix(Int,Nz,Nx); //Ny-> Nz
  free_float_3d(Vl,Nz,Ny,Nx);
  free_float_3d(Vu,Nz,Ny,Nx);
  
  return 0;
}
