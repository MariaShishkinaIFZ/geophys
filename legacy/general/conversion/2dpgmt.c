#include <stdio.h>
#include <math.h>

#include "file.h"
#include "submem.h"
#include "gmt_nan.h"

#define STR_LN 256

int main(int narg, char **argv)
{
  int M, N;
  int i,j;
  float **V;
  float long_0, long_E, latt_0, latt_E, d_long, d_latt, min_z, max_z;
  struct
  {
    int nx;
    int ny;
    int node_offset;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
    double x_inc;
    double y_inc;
    double z_scale_factor;
    double z_add_offset;
    char x_units[80];
    char y_units[80];
    char z_units[80];
    char title[80];
    char command[320];
    char remark[160];
  } ncg_h;

  FILE *fi, *fo;
  int verbose=0;
  int centered=0;
  int no_conv=0;
  char label[STR_LN];

  float GMT_f_NaN;

  GMT_make_fnan(GMT_f_NaN);

  if (narg<8)
  {
    printf("Insufficient parameters. USAGE:\n2dgmt <2d data file> <NCG out file> Nx Ny x0 y0 h\n");
    exit(1);
  }
  if (narg>3)
  {
    for (i=3;i<narg;i++)
    {
      if (strcmp(argv[i],"-v")==0) verbose=1;
    }
  }

  sscanf(argv[3],"%d",&N);
  sscanf(argv[4],"%d",&M);
  sscanf(argv[5],"%g",&long_0);
  sscanf(argv[6],"%g",&latt_0);
  sscanf(argv[7],"%g",&d_latt);
  d_long=d_latt;

  latt_E=latt_0+M*d_latt;
  long_E=long_0+N*d_long;

  if (verbose)
  {
      printf("Data info:\n");
      printf("Dimensions: %d rows x %d cols\n",M,N);
      printf("Latt: %g - %g by %g units.\n",latt_0,latt_E,d_latt);
      printf("Long: %g - %g by %g units.\n",long_0,long_E,d_long);
      printf("Values are treated as pixels.\n");
      printf("\n");
  }

  alloc_float_matrix(&V,M,N);

  fi=f_open("./",argv[1],"rb");
  for (i=0;i<M;i++)
  {
    fread(V[i],sizeof(float),N,fi);
  }
  fclose(fi);

  ncg_h.nx=N;
  ncg_h.ny=M;
  ncg_h.node_offset=1;

  ncg_h.x_min=(double)long_0;
  ncg_h.x_max=(double)long_E;
  ncg_h.y_min=(double)latt_0;
  ncg_h.y_max=(double)latt_E;
  ncg_h.x_inc=(double)d_long;
  ncg_h.y_inc=(double)d_latt;
  ncg_h.z_scale_factor=1.0;
  ncg_h.z_add_offset=0.0;
  strcpy(ncg_h.x_units,"unknown");
  strcpy(ncg_h.y_units,"unknown");
  strcpy(ncg_h.z_units,"unknown");
  strcpy(ncg_h.title,"Grid from 2d file");
  strcpy(ncg_h.command,argv[0]);
  for (i=1;i<narg;i++) {strcat(ncg_h.command," "); strcat(ncg_h.command,argv[i]);}
  strcpy(ncg_h.remark,"");

  ncg_h.z_min=1.7e38;
  ncg_h.z_max=-1.7e38;
  for (i=0;i<M;i++)
    for (j=0;j<N;j++)
    {
      if (V[i][j]>1e9) GMT_make_fnan(V[i][j]);
      else
      {
        if (ncg_h.z_min>V[i][j]) ncg_h.z_min=V[i][j];
        if (ncg_h.z_max<V[i][j]) ncg_h.z_max=V[i][j];
      }
    }

  fo=f_open("./",argv[2],"wb");
  fwrite(&ncg_h,1,sizeof(ncg_h),fo);
  for (i=ncg_h.ny-1;i>=0;i--)
    fwrite(V[i],ncg_h.nx,sizeof(float),fo);
  fclose(fo);

  free_float_matrix(V,M,N);
  if (verbose) printf("Success.\n");
}
