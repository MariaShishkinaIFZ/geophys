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

/*x0 - начало отсчёта, w - ширина клетки, t - полуширина переходной зоны*/
float ch1d(float x0, float y0, float z0, float R, float t,float x, float y, float z)
{
  long a = 1;
  float r=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
  
  if (r<R+t)
  {
    if (r<=R-t) return a;
    else return a*(sin(M_PI*(r+2*t-R)/(2.0*t))/2+0.5);
  }
  else
    return EMPTY_VALUE;
}



int main(int narg, char **argv)
{
  long int Nx, Ny, Nz; 
  int i, j, k, ask;
  float x0, y0, z0, h, cR, ct;
  float cx0, cy0, cz0;
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
	if (ask) fprintf(stdout,"Please enter y0, zO, R and t - wide of transit  > \n");
	fscanf(fp,"%g%g%g%g",&cy0,&cz0,&cR,&ct); cx0=0; 
	if (ask) fprintf(stdout,"Axis %c \n",cs); break;
	
      case 'Y': 
	if (ask) fprintf(stdout,"Please enter x0, zO, R and t - wide of transit  > \n");
	fscanf(fp,"%g%g%g%g",&cx0,&cz0,&cR,&ct); cy0=0; break;
      case 'Z': 
	if (ask) fprintf(stdout,"Please enter x0, yO, R and t - wide of transit  > \n");
	fscanf(fp,"%g%g%g%g",&cx0,&cy0,&cR,&ct); cz0=0; break;
      default: errquit("Incorrect label for slice position. MUST be either X, Y or Z.\n");
    }
  if (ask) fprintf(stdout,"Please enter the absolute velocity of the inclusion > \n");
  fscanf(fp,"%g",&abs_vel);
  if (ask) fprintf(stdout,"Please enter chessed model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",fout);
  fclose(fp);
  
  
  if (ask) fprintf(stdout,"x0 %g \n",cx0);
  if (ask) fprintf(stdout, "y0 %g \n",cy0);
  if (ask) fprintf(stdout, "z0 %g \n",cz0);
  if (ask) fprintf(stdout, "R %g \n",cR);
  if (ask) fprintf(stdout, "t %g \n",ct);
    
  alloc_float_3d(&V,Nz,Ny,Nx);
  read_3D_volume(V,Nz,Ny,Nx,vm_file);
  
  for (k=0;k<Nz;k++)
    for (i=0;i<Ny;i++)
      for (j=0;j<Nx;j++)
      {
	switch (cs)
	{
	  
	 case 'X': incl_C=ch1d(x0+h*j,cy0,cz0,cR,ct,x0+h*j,y0+h*i,z0+h*k); break;
	 case 'Y': incl_C=ch1d(cx0,y0+h*i,cz0,cR,ct,x0+h*j,y0+h*i,z0+h*k); break;
	 case 'Z': incl_C=ch1d(cx0,cy0,z0+h*k,cR,ct,x0+h*j,y0+h*i,z0+h*k); break;
	}
        if (incl_C<EMPTY_LIMIT) V[k][i][j]+=(abs_vel-V[k][i][j])*(incl_C);
      }
      
  write_3D_volume(V,Nz,Ny,Nx,fout);
  free_float_3d(V,Nz,Ny,Nx);
  exit (1);
}
