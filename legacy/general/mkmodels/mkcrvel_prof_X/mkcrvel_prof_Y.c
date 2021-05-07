#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "submem.h"
#include "subio.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN  256	/*Maximum length of the file name*/

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
  int Nx, Ny, Nz, i, j, k, l, Nl;
  int ask;
  float ***Vmod;
  float xc, yc, zc, V0, gradV;
  float x0, y0, z0, h;
  float x, y, z;
  float *VP, *H;
  char outfile[FLNLN];
  FILE *fo, *fp;
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%d%d%d",&Nx,&Ny,&Nz);
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
  
  if (ask) fprintf(stdout,"Input number of levels > \n");
  fscanf(fp,"%d",&Nl);
  VP=fmemalloc(Nl,sizeof(float));
  H=fmemalloc(Nl,sizeof(float));
  for (i=0;i<Nl;i++)
  {
    if (ask) fprintf(stdout,"Input H VP > \n");
    fscanf(fp,"%g%g",H+i,VP+i);
    fprintf(stderr,"%d %g %g",i,H[i],VP[i]);
  }
  if (ask) fprintf(stdout,"Input velocity model filename (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",outfile);
  if (!ask) fclose(fp);
  
  alloc_float_3d(&Vmod,Nz,Ny,Nx);
  
  for (k=0;k<Nz;k++)
  {
    z=z0+k*h;
    for (j=0;j<Ny;j++)
    {
      y=y0+j*h;
      for (i=0;i<Nx;i++)
      {
	if (y<=H[0])
	  Vmod[k][j][i]=VP[0];
	else
	{
	  if (y>H[Nl-1])
	    Vmod[k][j][i]=VP[Nl-1];
	  else
	  {
	    for (l=0;l<Nl-1;l++)
	      if (y<=H[l+1]) break;
	    if (i==0 && j==0) fprintf(stderr,"%d",l);
	    Vmod[k][j][i]=VP[l]+(VP[l+1]-VP[l])/(H[l+1]-H[l])*(y-H[l]);
	  }
        }
      }
    }
  }

  fmemfree(H);
  fmemfree(VP);

  write_3D_volume(Vmod,Nz,Ny,Nx,outfile);  
  free_float_3d(Vmod,Nz,Ny,Nx);
  fprintf(stderr,"Success.\n");
}
