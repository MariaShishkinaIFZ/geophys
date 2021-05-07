#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 80 	/*Maximum length of the file name*/

/*x0 - начало отсчёта, w - ширина клетки, t - полуширина переходной зоны*/
double ch1d(double x, double x0, double w, double t)
{
  long n = floor((x-x0)/w); /*номер клетки*/
  long a = pow(-1.0,n);
  
  if ((x-x0>=n*w+t) && (x-x0<=(n+1)*w-t)) return a;
  else
  {
    if (x-x0<n*w+t) return a*sin(M_PI*(x-x0-n*w)/(2.0*t));
    else return a*sin(M_PI*((n+1)*w-(x-x0))/(2.0*t));
  }
}

int main(int narg, char **argv)
{
  int Nx, Ny, Nz, i, j;
  float x0, y0, z0, h;
  float cx0, cxw, cxt, cy0, cyw, cyt;
  float x, y, z;
  float **V;
  char fout[FLNLN];
  FILE *fp;
  
  if (narg<2)
  {
    fprintf(stderr,"Use: mkchess2d <parameters file name>\n");
    exit(0);
  }
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
 
  fp=file_open(argv[1],"rt");
  fscanf(fp,"%g%g%g",&cx0,&cxw,&cxt);
  fscanf(fp,"%g%g%g",&cy0,&cyw,&cyt);
  fscanf(fp,"%s",fout);
  fclose(fp);
  
  alloc_float_matrix(&V,Ny,Nx);
  
  for (i=0;i<Ny;i++)
    for (j=0;j<Nx;j++)
      V[i][j]=ch1d(h*i,cy0,cyw,cyt)*ch1d(h*j,cx0,cxw,cxt);
      
  write_2D_surface(V,Ny,Nx,fout);
  free_float_matrix(V,Ny,Nx);
}
