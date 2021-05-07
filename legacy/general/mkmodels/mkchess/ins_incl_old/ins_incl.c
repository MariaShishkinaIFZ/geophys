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
float ch1d(float x, float x0, float w, float t)
{
  long n = floor((x-x0)/w); /*номер клетки*/
  long a = 1;
  
  if (n==0)
  {
    if ((x-x0>=n*w+t) && (x-x0<=(n+1)*w-t)) return a;
    else
    {
      if (x-x0<n*w+t) return a*sin(M_PI*(x-x0-n*w)/(2.0*t));
      else return a*sin(M_PI*((n+1)*w-(x-x0))/(2.0*t));
    }
  }
  else
    return EMPTY_VALUE;
}

int main(int narg, char **argv)
{
  long int Nx, Ny, Nz; 
  int i, j, k, ask;
  float x0, y0, z0, h;
  float cx0, cxw, cxt, cy0, cyw, cyt, cz0, czw, czt, cx1, cy1, cz1;
  float abs_vel;
  float ***V;
  float incl_X, incl_Y, incl_Z;
  char fout[FLNLN], vm_file[FLNLN];
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
  fscanf(fp,"%s",vm_file);
  if (ask) fprintf(stdout,"Please enter x0, x1, X wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cx0,&cx1,&cxt); cxw=cx1-cx0;
  if (ask) fprintf(stdout,"Please enter y0, y1, Y wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cy0,&cy1,&cyt); cyw=cy1-cy0;
  if (ask) fprintf(stdout,"Please enter z0, z1, Z wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cz0,&cz1,&czt); czw=cz1-cz0;
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
        incl_Z=ch1d(z0+h*k,cz0,czw,czt);
        incl_Y=ch1d(y0+h*i,cy0,cyw,cyt);
        incl_X=ch1d(x0+h*j,cx0,cxw,cxt);
        if (incl_X<EMPTY_LIMIT && incl_Y<EMPTY_LIMIT && incl_Z<EMPTY_LIMIT)
        V[k][i][j]+=(abs_vel-V[k][i][j])*incl_X*incl_Y*incl_Z;
      }
      
  write_3D_volume(V,Nz,Ny,Nx,fout);
  free_float_3d(V,Nz,Ny,Nx);
  exit (1);
}
