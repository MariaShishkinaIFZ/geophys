#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "submem.h"
#include "subio.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/

void write_3D_volume(float ***Buf, long Nz, long Ny, long Nx, char *filename)
{
  long j,k;
  FILE *fo;
  
  fo=file_open(filename,"wb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fwrite(Buf[k][j],sizeof(float),Nx,fo);
  fclose(fo);
}

int main (int narg, char **argv)
{
  int Nx, Ny, Nz, i, j, k;
  int ask;
  float ***Vmod;
  float xc, yc, zc, V0, gradV;
  float x0, y0, z0, h;
  float x, y, z;
  char outfile[FLNLN];
  FILE *fo, *fp;
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
  
  if (narg < 2) errquit("Insufficient parameters.\n");
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
  
  if (ask) fprintf(stdout,"Input V0 gradV > \n");
  fscanf(fp,"%g%g",&V0,&gradV);
  if (ask) fprintf(stdout,"Input velocity model filename (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",outfile);
  if (!ask) fclose(fp);
  
   fprintf(stdout,"\n V0= %g >  \n",V0);
   fprintf(stdout,"\n gradV= %g >  \n",gradV);


  alloc_float_3d(&Vmod,Nz,Ny,Nx);
  
  for (k=0;k<Nz;k++)
  {
    z=z0+k*h;
    for (j=0;j<Ny;j++)
    {

      for (i=0;i<Nx;i++)
      {
	Vmod[k][j][i]=V0+z*gradV;
      }
    }
  }
write_3D_volume(Vmod,Nz,Ny,Nx,outfile);  
  free_float_3d(Vmod,Nz,Ny,Nx);
}
