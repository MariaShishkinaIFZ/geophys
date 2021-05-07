#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/

#define EMPTY_VALUE 1e10
#define EMPTY_LIMIT 1e9

/*x0, y0 - центр, a,b - полуоси, mag - амплитуда, t - полуширина переходной зоны, x,y - текущие координаты*/
float ellips(float x0, float y0, float a, float b, float t, float x, float y)
{
  float x1 = (x-x0);
  float y1 = (y-y0);
  
  float q = (x1*x1/a+y1*y1/b);
  
  if (q > 1+t)
    return EMPTY_VALUE;
  else
  {
    if (q < 1-t) return 1;
    else return sin(M_PI/2.0*((q-1+t)/(2*t)+1));
  }
}


int main(int narg, char **argv)
{
  long int Nx, Ny, Nz; 
  int i, j, k, ask;
  float x0, y0, z0, h, ct;
  float cx0, cy0, cz0, a, b;
  float abs_vel;
  float ***V;
  float incl_C;
  char fout[FLNLN], vm_file[FLNLN];
  char cs;
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

  if (ask) fprintf(stdout,"Please enter velocity model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s%*c",vm_file);
  if (ask) fprintf(stdout,"Please enter axis to include cylinder along X or Y or Z > \n");
  fscanf(fp,"%c",&cs);
    switch (cs)
    {
      case 'X':   
	if (ask) fprintf(stdout,"Please enter y0, zO, a, b and t - wide of transit  > \n");
	fscanf(fp,"%g%g%g%g%g",&cy0,&cz0,&a,&b,&ct); cx0=0; break;
      case 'Y': 
	if (ask) fprintf(stdout,"Please enter x0, zO, a, b, and t - wide of transit  > \n");
	fscanf(fp,"%g%g%g%g%g",&cx0,&cz0,&a,&b,&ct); cy0=0; break;
      case 'Z': 
	if (ask) fprintf(stdout,"Please enter x0, yO, a, b, and t - wide of transit  > \n");
	fscanf(fp,"%g%g%g%g%g",&cx0,&cy0,&a,&b,&ct); cz0=0; break;
      default: errquit("Incorrect label for slice position. MUST be either X, Y or Z.\n");
    }
  if (ask) fprintf(stdout,"Please enter the absolute velocity of the inclusion > \n");
  fscanf(fp,"%g",&abs_vel);
  if (ask) fprintf(stdout,"Please enter chessed model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",fout);
  fclose(fp);
    
  alloc_float_3d(&V,Nz,Ny,Nx);
  read_3D_volume(V,Nz,Ny,Nx,vm_file);
  
  for (k=0;k<Nz;k++)
    for (i=0;i<Ny;i++)
      for (j=0;j<Nx;j++)
      {
	switch (cs)
	{
	  
	 case 'X': incl_C=ellips(cy0,cz0,a,b,ct,y0+h*i,z0+h*k);
	 case 'Y': incl_C=ellips(cx0,cz0,a,b,ct,x0+h*j,z0+h*k);
	 case 'Z': incl_C=ellips(cx0,cy0,a,b,ct,x0+h*j,y0+h*i);
	}
        if (incl_C<EMPTY_LIMIT) V[k][i][j]+=(abs_vel-V[k][i][j])*(incl_C);
      }
      
  write_3D_volume(V,Nz,Ny,Nx,fout);
  free_float_3d(V,Nz,Ny,Nx);
  exit (1);
}
