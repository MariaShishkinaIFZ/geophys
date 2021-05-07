#include <iostream> 
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstring>
#include <limits>

#include "subio.h"
#include "submem.h"
#include "nrutil.h"

#define GEOMFILE "model-geometry.gras"
#define COARSE_UL_GEOMFILE "model-geometry.coarse"

#define FLNLN 80 	/*Maximum length of the file name*/
#define DUMMY_LIM 1e9
#define DUMMY_VALUE 1e10
#define STRLN 256

#define ACC_CONST 1e-5

int av3d_ins_lay(float ***src,float **top,float **bot, float ***dest,
                 long nlay, long nrow, long ncol, long kav,
		 long kmin, long kz, long NL,
		 float z0, float h)
{
  long i, j, k, k1, k2, ii, jj, kk;
  long M, N;
  float w, z;

  N=(ncol-1)/kav;
  M=(nrow-1)/kav;

  if (((N*kav)!=(ncol-1)) || ((M*kav)!=(nrow-1)))
    return -1; /*Error*/

  for (kk=1;kk<=NL;kk++)
  {
   k1=kmin+(kk-1)*kz-1;
   k2=kmin+kk*kz-2;
   //(fprintf(stderr,"kk=%ld k1=%ld k2=%ld\n",kk,k1,k2));
   for (jj=1;jj<=M;jj++)
    for (ii=1;ii<=N;ii++)
    {
      dest[kk][jj][ii]=0.0; w=0.0;
      /*Internal layers*/
      for (k=k1+1;k<=k2-1;k++)
      {
       z=z0+k*h;
       /*Internal cell points*/
       for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
        for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
        {
	      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
	      {
	          fprintf(stderr,"1: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	      }
	      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h))
	      {
            dest[kk][jj][ii]+=src[k][j][i];
	        w+=1.0;
	      }
        }
       /*Front cell boundary*/
       j=(jj-1)*kav;
       i=(ii-1)*kav;
       if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"2: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       i=ii*kav-1;
       if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"3: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
       /*Rear cell boundary*/
       j=jj*kav-1;
       i=(ii-1)*kav;
       if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"4: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       i=ii*kav-1;
       if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"5: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
       /*Left cell boundary*/
       i=(ii-1)*kav;
       if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"6: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
       for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
       /*Right cell boundary*/
       i=ii*kav-1;
       if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"7: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
       for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
      }

      k=k1; /*Cell top*/
      z=z0+k*h;
      /*Internal cell points*/
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
       {
	     if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
         {
	          fprintf(stderr,"8: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	     }
	     if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h))
	     {
           dest[kk][jj][ii]+=0.5*src[k][j][i];
	       w+=0.5;
	     }
       }
      /*Front cell boundary*/
      j=(jj-1)*kav;
      i=(ii-1)*kav;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"9: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav-1;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"10: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Rear cell boundary*/
      j=jj*kav-1;
      i=(ii-1)*kav;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"11: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav-1;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"12: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Left cell boundary*/
      i=(ii-1)*kav;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"13: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Right cell boundary*/
      i=ii*kav-1;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"14: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}

      k=k2; /*Cell bottom*/
      z=z0+k*h;
      /*Internal cell points*/
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
       {
	     if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
         {
	          fprintf(stderr,"15: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	     }
	     if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h))
	     {
           dest[kk][jj][ii]+=0.5*src[k][j][i];
	       w+=0.5;
	     }
       }
      /*Front cell boundary*/
      j=(jj-1)*kav;
      i=(ii-1)*kav;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"16: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav-1;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"17: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Rear cell boundary*/
      j=jj*kav-1;
      i=(ii-1)*kav;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"18: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"19: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Left cell boundary*/
      i=(ii-1)*kav;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"20: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Right cell boundary*/
      i=ii*kav-1;
      if ((i<0) || (i>=ncol) || (j<0) || (j>=nrow))
       {
	          fprintf(stderr,"21: jj=%ld ii=%ld kav=%ld i=%ld j=%ld nrow=%ld ncol=%ld",jj,ii,kav,i,j,nrow,ncol);
	          return -1;
	   }
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}

      if (w>0)
        dest[kk][jj][ii]/=w;
      else
        dest[kk][jj][ii]=DUMMY_VALUE;
    }
   } /*for kk*/
   return 1;
}

long roun(double x)
{
  long i;
  i = floor(x);
  if (x-i >= 0.5) i++;
  return i;
}

void read_3D_volume(float ***Buf, long Nz, long Ny, long Nx, const char *filename)
{
  long j,k;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fread(Buf[k][j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_3D_volume_from_1(float ***Buf, long Nz, long Ny, long Nx, const char *filename)
{
  long j,k;
  FILE *fo;

  fo=file_open(filename,"wb");
  for (k=1;k<=Nz;k++)
    for (j=1;j<=Ny;j++)
      fwrite(Buf[k][j]+1,sizeof(float),Nx,fo);
  fclose(fo);
}

void read_2D_surface(float **Buf, long Ny, long Nx, const char *filename)
{
  long j;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (j=0;j<Ny;j++)
    fread(Buf[j],sizeof(float),Nx,fi);
  fclose(fi);
}

int main(int narg, char **argv)
{
  long Nx, Ny, Nz, NNx, NNy, NNz, N;
  long k_gen, k_gen_v;
  long uk_min, uk_max;
  float x0,y0,z0;
  float h,dd,ddz;
  float ***VFINE;	/*Fine gridded velocity model (nodes from 0,0,0)*/
  float ***UL_dscm;	/*Upper physical layer 3D coarse grid cells slowness model (cells from 1,1,1)*/
  float **UL_ft;	/*Upper physical layer top fine grid position (nodes from 0,0)*/
  float **Zintf;	/*Interface fine grid position (nodes from 0,0)*/
  float Tmin, Zmax;	/*Minimum and maximum of the interface position*/
  float UL_H_min;	/*Minimum upper physical layer model depth*/
  char fi_VelMod[FLNLN], fi_fiBas[FLNLN], fi_fiTop[FLNLN], fo_VelMod[FLNLN];
  char buf[STRLN];
  FILE *fp;


  fp = file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fscanf(fp,"%ld",&k_gen);
  fscanf(fp,"%ld",&k_gen_v);
  fclose(fp);
  
  fp=file_open(COARSE_UL_GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&NNx,&NNy,&NNz);
  fscanf(fp,"%ld%ld",&uk_min,&uk_max);
  fscanf(fp,"%g",&UL_H_min);
  fscanf(fp,"%g%g",&dd,&ddz);
  fclose(fp);
  
  N=NNx*NNy;

  /*------------Read parameters----------------*/
  if (narg < 2) errquit("Insufficient parameters.\n");
  if (strcmp(argv[1],"-f")==0) /*Use parameters file with the name given as 2nd parameter*/
  {
    if (narg < 3) errquit("Parameters file name not given.\n");
    fp=file_open(argv[2],"rt");
  }
  else
  {
    if (strcmp(argv[1],"-s")==0) /*Read parameters from standart input*/
      fp=stdin;
    else
    {
      printf("Unrecognized flag %s.\n",argv[1]);
      exit(1);
    }
  }

  fprintf(stderr,"Reading parameters.\n");
  /*1st section: reference model*/
  fscanf(fp,"%s",fi_VelMod); fgets(buf,STRLN,fp);  /*1st LINE: Fine velocity model*/
  fscanf(fp,"%s",fi_fiTop);	 fgets(buf,STRLN,fp);  /*2nd LINE: Fine layer top position*/
  fscanf(fp,"%s",fi_fiBas);	 fgets(buf,STRLN,fp);  /*3rd LINE: Fine reference interface position*/
  fscanf(fp,"%s",fo_VelMod);	fgets(buf,STRLN,fp);  /*4th LINE: Output: coarse velocity model*/
  fprintf(stderr,"Parameters read.\n");


  if (strcmp(argv[1],"-f")==0) fclose(fp);
  /*--------------------------------------------------------*/

  fprintf(stderr,"Allocating memory...");
  alloc_float_3d(&VFINE,Nz,Ny,Nx);
  alloc_float_matrix(&Zintf,Ny,Nx);
  alloc_float_matrix(&UL_ft,Ny,Nx);
  fprintf(stderr,"Success.\n");

   /*------Read reference true velocity model-------------------------------------------*/
  fprintf(stderr,"Reading fine velocity model...");
  read_3D_volume(VFINE,Nz,Ny,Nx,fi_VelMod);

  /*------Read fine grids of layer top and interface*/
  fprintf(stderr,"Read surfaces...");
  read_2D_surface(UL_ft,Ny,Nx,fi_fiTop);
  read_2D_surface(Zintf,Ny,Nx,fi_fiBas);
  fprintf(stderr,"Success.\n");

  /*Averaging*/
 
  fprintf(stderr,"%ld cells with dimentions of %gx%g units in each layer.\n",N,dd,dd);
  fprintf(stderr,"layer model span: k %ld-%ld, Z %g-%g\n",uk_min,uk_max,
  		z0+h*(uk_min-1),z0+h*(uk_max-1));
  fprintf(stderr,"%ld grid layers with %g units thickness each.\n",NNz,ddz);

  fprintf(stderr,"Allocating memory for UL_dscm (NNz=%ld NNy=%ld NNx=%ld)...\n",NNz,NNy,NNx);
  UL_dscm=f3tensor(1,NNz,1,NNy,1,NNx);
  fprintf(stderr,"Memory for UL_dscm allocated successfully.\n");

  fprintf(stderr,"Averaging VFINE...\n");
  if (av3d_ins_lay(VFINE,UL_ft,Zintf,UL_dscm,Nz,Ny,Nx,k_gen,uk_min,k_gen_v,NNz,z0,h)!=1)
  return EXIT_FAILURE;
  fprintf(stderr,"Vfine averaged.\n");

  write_3D_volume_from_1(UL_dscm,NNz,NNy,NNx,fo_VelMod);
  fprintf(stderr,"Success.\n");
}
