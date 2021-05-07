#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "submem.h"

#define STRLN 80

int gras_round(float x)
{
  return (int)(x-floor(x) >= 0.5 ? ceil(x) : floor(x));
}

int main(int narg, char **argv)
{
  int Np;
  int Nx, Ny, Nz;
  int i, j, k, l, ind;
  float x0,y0,z0;
  float *TF;
  float *xp, *yp, *zp;
  int *np;
  float t;
  float h;
  char tfln[STRLN], outfln[STRLN];
  FILE *fi, *fo;

  if (narg<2)
  {
    printf("Insufficient parameters.\nUSAGE:\ngetpeaks <parameters file>\n");
    exit(1);
  }
  /*--------------------Read parameters-------------------*/
  fi=fopen(argv[1],"rt");
  fscanf(fi,"%s",tfln);
  fscanf(fi,"%s",outfln);
  fscanf(fi,"%d%d%d",&Nx,&Ny,&Nz);
  fscanf(fi,"%g%g%g",&x0,&y0,&z0);
  fscanf(fi,"%g",&h);
  fscanf(fi,"%d",&Np);
  np=(int*)fmemalloc(Np,sizeof(int));
  xp=(float*)fmemalloc(Np,sizeof(float));
  yp=(float*)fmemalloc(Np,sizeof(float));
  zp=(float*)fmemalloc(Np,sizeof(float));
  for (i=0;i<Np;i++)
    fscanf(fi,"%d%g%g%g",np+i,xp+i,yp+i,zp+i);
  fclose(fi);
  printf("Parameters read.\n");
  /*------------------------------------------------------*/
  printf("Time field array size %dx%dx%d = %d elements",Nx,Ny,Nz,Nx*Ny*Nz);

  TF=fmemalloc(Nx*Ny*Nz,sizeof(float));
  printf("Memory allocated\n");
  fi=fopen(tfln,"rb");
  for (i=0;i<Nz;i++)
  {
    fread(&TF[i*Nx*Ny],sizeof(float),Nx*Ny,fi);
  }
  fclose(fi);
  printf("Timefield read...");

  fo=fopen(outfln,"wt");
  for (l=0;l<Np;l++)
  {
    i=gras_round((yp[l]-y0)/h);
    j=gras_round((xp[l]-x0)/h);
    k=gras_round((zp[l]-z0)/h);
    if ((i>=0) && (i<Ny) && (j>=0) && (j<Nx) && (k>=0) && (k<Nz))
    {
      ind=k*Nx*Ny+i*Nx+j;
      t=TF[ind];
      fprintf(fo,"%5d %10g %10g %10g %5d %5d %5d %10g\n",np[l],xp[l],yp[l],zp[l],i,j,k,t);
    }
    else
    {
      fprintf(fo,"%5d %5d %5d %5d Out of volume\n",np[l],j,i,k);
    }
  }
  fclose(fo);

  fmemfree(TF);
  fmemfree(zp);
  fmemfree(yp);
  fmemfree(xp);
  fmemfree(np);

  return 1;
}
