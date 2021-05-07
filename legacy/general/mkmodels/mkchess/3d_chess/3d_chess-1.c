#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/

/*x0 - начало отсчёта, w - ширина клетки, t - полуширина переходной зоны*/
double ch1d(double x, double x0, double w, double t)
{
  long n = floor((x-x0)/w); /*номер клетки*/
  long a = pow(-1.0,n);
  
  if (w<=0) return 1.0;
  if ((x-x0>=n*w+t) && (x-x0<=(n+1)*w-t)) return a;
  else
  {
    if (x-x0<n*w+t) return a*sin(M_PI*(x-x0-n*w)/(2.0*t));
    else return a*sin(M_PI*((n+1)*w-(x-x0))/(2.0*t));
  }
}

int main(int narg, char **argv)
{
  long int Nx, Ny, Nz; 
  int i, j, k, ask;
  float x0, y0, z0, h;
  float cx0, cxw, cxt, cy0, cyw, cyt, cz0, czw, czt;
  float rel_anom;
  float ***V;
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
  if (ask) fprintf(stdout,"Please enter x0, wide of cell in X dimension, X wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cx0,&cxw,&cxt);
  if (ask) fprintf(stdout,"Please enter y0, wide of cell in Y dimension, Y wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cy0,&cyw,&cyt);
  if (ask) fprintf(stdout,"Please enter z0, wide of cell in Z dimension, Z wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cz0,&czw,&czt);
  if (ask) fprintf(stdout,"Please enter relative amplitude of anomaly in part of one > \n");
  fscanf(fp,"%g",&rel_anom);
  if (ask) fprintf(stdout,"Please enter chessed model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",fout);
  fclose(fp);
  
  alloc_float_3d(&V,Nz,Ny,Nx);
  read_3D_volume(V,Nz,Ny,Nx,vm_file);
  
  for (k=0;k<Nz;k++)
    for (i=0;i<Ny;i++)
      for (j=0;j<Nx;j++)
        V[k][i][j]=V[k][i][j]*
                    (1+rel_anom*ch1d(z0+h*k,cz0,czw,czt)*ch1d(y0+h*i,cy0,cyw,cyt)*ch1d(x0+h*j,cx0,cxw,cxt));
      
  write_3D_volume(V,Nz,Ny,Nx,fout);
  free_float_3d(V,Nz,Ny,Nx);
  exit (1);
}
