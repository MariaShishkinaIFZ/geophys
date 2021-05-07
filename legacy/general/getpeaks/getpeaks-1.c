#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "submem.h"
#include "subio.h"

#define GEOMFILE "model-geometry.gras"
#define STRLN 256
#define FLNLN 256 	/*Maximum length of the file name*/

typedef struct 
	{ 
	  int nprof;
	  long long np;
	  float xp,yp,zp;
	}
	string;

int gras_round(float x)

{
  return (int)(x-floor(x) >= 0.5 ? ceil(x) : floor(x));
}

int main(int narg, char **argv)
{
  int Np;
  int Nx, Ny, Nz;
  int i, j, k, l, ind;
  int ask;
  float x0,y0,z0;
  float *TF;
  string *rec_data;
  float t;
  float h;
  char tfln[FLNLN], outfln[FLNLN], parfln[FLNLN];
  FILE *fi, *fo;

  if (narg<2)
  {
    printf("Insufficient parameters.\nUSAGE:\ngetpeaks -f <parameters file>\n");
    exit(0);
  }
  /*--------------------Read parameters-------------------*/
  fi=file_open(GEOMFILE,"rt");
  fscanf(fi,"%d%d%d",&Nx,&Ny,&Nz);
  fscanf(fi,"%g%g%g",&x0,&y0,&z0);
  fscanf(fi,"%g",&h);
  fclose(fi);

  if (strcmp(argv[1],"-f")==0) /*Use parameters file with the name given as 2nd parameter*/
  {
    if (narg < 3) errquit("Parameters file name not given.\n");
    fi=file_open(argv[2],"rt");
    ask=0;
  }
  else
  {
    if (strcmp(argv[1],"-s")==0) /*Read parameters from standart input*/
    {
      fi=stdin;
      ask=1;
    }
    else
    {
      printf("Unrecognized flag %s.\n",argv[1]);
      exit(1);
    }
  }
  
  if (ask) fprintf(stdout,"Input model filename.3d (%d chars max) > \n ",FLNLN);
  fscanf(fi,"%s",tfln);
  if (ask) fprintf(stdout,"Input peaks filename.dat (%d chars max) > \n ",FLNLN);
  fscanf(fi,"%s",outfln);
  if (ask) fprintf(stdout,"Input points position data filename.dat (%d chars max) > \n ",FLNLN);
  fscanf(fi,"%s",parfln);
  if (ask) fi=file_open(parfln,"rt");

  fscanf(fi,"%d",&Np);
  rec_data=(string*)fmemalloc(Np,sizeof(string));
  
  for (i=0;i<Np;i++)
    fscanf(fi,"%d%Ld%g%g%g",&rec_data[i].nprof,&rec_data[i].np, &rec_data[i].xp, &rec_data[i].yp, &rec_data[i].zp);
  fclose(fi);
  printf("Parameters read.\n");
  /*------------------------------------------------------*/
  printf("Time field array size %dx%dx%d = %d elements\n",Nx,Ny,Nz,Nx*Ny*Nz);

  TF=fmemalloc(Nx*Ny*Nz,sizeof(float));
  printf("Memory allocated\n");
  fi=file_open(tfln,"rb");
  for (i=0;i<Nz;i++)
  {
    fread(&TF[i*Nx*Ny],sizeof(float),Nx*Ny,fi);
  }
  fclose(fi);
  printf("Timefield read...\n");

  fo=fopen(outfln,"wt");
  fprintf(fo,"%d \n",Np);
  for (l=0;l<Np;l++)
  {
    i=gras_round((rec_data[l].yp-y0)/h);
    j=gras_round((rec_data[l].xp-x0)/h); 
    k=gras_round((rec_data[l].zp-z0)/h);
    if ((i>=0) && (i<Ny) && (j>=0) && (j<Nx) && (k>=0) && (k<Nz))
    {
      ind=k*Nx*Ny+i*Nx+j;
      t=TF[ind];
      fprintf(fo,"%d %13Ld %10g %10g %10g %5d %5d %5d %10g\n",rec_data[l].nprof,rec_data[l].np,rec_data[l].xp,rec_data[l].yp,rec_data[l].zp,i,j,k,t);
    }
    else
    {
      fprintf(fo,"%d %13Ld %5d %5d %5d Out of volume\n",rec_data[l].nprof,rec_data[l].np,j,i,k);
    }
  }
  fclose(fo);

  fmemfree(TF);
  fmemfree(rec_data);

  return 1;
}
