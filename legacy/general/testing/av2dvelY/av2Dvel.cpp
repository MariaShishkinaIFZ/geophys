//Программа берет 3D скоростную модель по детельной сетке и создает Gmt файл 2D (вырезанный в середине) по грубой сетке

#include <iostream>
#include <cstdlib>

#include "gstdata.h"
#include "gsterror.h"
#include "gstdefs.h"
//#include "gstmath.h"
#include "gstio.h"
#include "subio.h"
#include "submem.h"
#include "nrutil.h"

#define GEOMFILE "model-geometry.gras"
#define DUMMY_LIM 1e9
#define DUMMY_VALUE 1e10
#define STRLN 256

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

int av2d(float **src, float **dest, long nrow, long ncol, long kav)
{
  long i, j, ii, jj;
  long M, N;
  float w;

  N=(ncol-1)/kav;
  M=(nrow-1)/kav;

  if (((N*kav)!=(ncol-1)) || ((M*kav)!=(nrow-1)))
    return -1; /*Error*/

  for (jj=1;jj<=M;jj++)
    for (ii=1;ii<=N;ii++)
    {
      dest[jj][ii]=0.0; w=0.0;
      /*Internal cell points*/
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
        for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
	  if (src[j][i]<DUMMY_LIM)
	  {
            dest[jj][ii]+=src[j][i];
	    w+=1.0;
	  }
      /*Lower cell boundary*/
      j=(jj-1)*kav;
      i=(ii-1)*kav;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      i=ii*kav-1;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      /*Upper cell boundary*/
      j=jj*kav-1;
      i=(ii-1)*kav;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      i=ii*kav-1;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      /*Left cell boundary*/
      i=(ii-1)*kav;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      /*Right cell boundary*/
      i=ii*kav-1;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      if (w>0)
        dest[jj][ii]/=w;
      else
        dest[jj][ii]=DUMMY_VALUE;
    }
    return 1;
}


int main(int narg, char ** argv)
{
  long i,j;
  long Nx,Ny,Nz, NNy, NNx, NNz;  
  long k_gen, k_gen_v;
  long z;
  
  float x0,y0,z0;
  float h;
  
  char f_3D[GST_STRLN];
  char f_2D_av[GST_STRLN]; 
  char buf[STRLN];
 
  float **V2D;            /*2D Fine gridded velocity model (nodes from 0,0)*/
  float ***VFINE;	/*Fine gridded velocity model (nodes from 0,0,0)*/
  float ** V2D_av ;     /* averagd  2D Fine gridded velocity model (nodes from 1,1)*/
 
  FILE *fp;   
  
 
  
 /*
  if (strcmp(argv[1],"-s")==0) //Read parameters from standart input
   {
   fp=stdin;
   fscanf(fp,"%s%s",f_3D,f_2D_av); 
   }
  else
  {
   if  (strcmp(argv[1],"-f")==0) //Use parameters file with the name given as 2nd parameter
     {
      fp=fopen(argv[1],"rt");  
      fscanf(fp,"%s%s",f_3D,f_2D_av);   
      fclose(fp);  
     }
   else
    {
      GSTerror("Unrecognized flag \n");
    }
  }
  
  */
  
//  fprintf(stderr,"%s %s.\n",f_3D,f_2D_av);  
  
  fp = fopen(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fscanf(fp,"%ld",&k_gen);
  fscanf(fp,"%ld",&k_gen_v);
  fclose(fp);
  
  if (k_gen != k_gen_v)
    errquit("Generalization constants must be equal for XY and Z\n");

  NNx=(Nx-1)/k_gen;
  NNy=(Ny-1)/k_gen;
  NNz=(Nz-1)/k_gen;

  if (((NNx*k_gen)!=(Nx-1)) || ((NNy*k_gen)!=(Ny-1)))
  {
    fprintf(stderr,"Invalid generalization constant.\n");
    exit(EXIT_FAILURE);
  }
  
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
  
  
  fscanf(fp,"%s",f_3D);		fgets(buf,STRLN,fp); /*1th LINE: Fine velocity model to averaging*/
  fscanf(fp,"%s",f_2D_av);	fgets(buf,STRLN,fp); /*2th LINE: Out *.gmt file*/

  if (strcmp(argv[1],"-f")==0) fclose(fp);
  /*--------------------------------------------------------*/ 
  
  
 //   if (narg<2)  GSTerror("Insufficient parameters \n");
  
 // fp=fopen(argv[1],"rt");  
 // fscanf(fp,"%s%s",f_3D,f_2D_av);   
 // fclose(fp);  
  
  
  
  
  fprintf(stderr,"Allocating memory...");
  alloc_float_3d(&VFINE,Nz,Ny,Nx);
  alloc_float_matrix(&V2D,Nz,Nx);
  V2D_av=matrix(1,Nz,1,Nx);  
  fprintf(stderr,"Success.\n");
  
  z=(Ny/2);
   
  
  read_3D_volume(VFINE,Nz,Ny,Nx,f_3D);  
  for (j=0; j <Nz; j++)
    for (i=0; i<Nx; i++)
     {
      V2D[j][i]=VFINE[j][z][i];  
      //V2D[j][i]=0.0;
     }     
   
   av2d(V2D,V2D_av,Nz,Nx,k_gen); 
   GSTgriddedData *V= new GSTgriddedData(NNz,NNx,-(z0+h*k_gen*(NNz)),-z0,x0, (x0+h*k_gen*(NNx)),area); // Was -z0,-(z0+h*k_gen*(NNz))
  
   for (j=0; j <NNz; j++)
    for (i=0; i<NNx; i++)
     {
      V->data()[NNz-1-j][i] =V2D_av[j+1][i+1];   
     }
  
    writeGMT( f_2D_av, *V);
    delete V;
   free_float_matrix(V2D,Nz,Nx);
   free_float_3d(VFINE,Nz,Ny,Nx);
   free_matrix(V2D_av,1,Nz,1,Nx);
   fprintf(stderr,"Success.\n");
}
