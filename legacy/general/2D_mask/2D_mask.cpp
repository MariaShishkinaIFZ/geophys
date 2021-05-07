 
/*
 программа вырезающая по файлу точек покрытия (2d файл по грубой сетке),
 файл глубин границ заданный на мелкой сетке.

-читаю модель и  2D  покрытие
- интерполирую покрытие
-записываю полько нужные (выерезанные маской значения)
 (верхние --1 оставляю как есть,
 Нижние дополняю знаячениеми из 2 D)

файл параметров- 
входные файлы:
-скоростная модель
-луечевое покрытие для верхнего слоя 3D
-лучевое покрытие для нижнего слоя 2D
файл с верхней поверхностью, 
                 файл с нижней поверхностью
*/


#include <iostream> 
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>
#include <limits>
#include <map>

#include <cstdio>
#include <cstring>
#include <values.h>

#include "subio.h"
#include "submem.h"
#include "nrutil.h"

#include "gmtio.h"

using namespace std;

#define GEOMFILE "model-geometry.gras"


//#define STATISTICS_FILE "gras.stat"

#define FLNLN 80 	/*Maximum length of the file name*/
#define STRLN 256
#define TSTRLN 256

#define DUMMY_LIM 1e9
#define DUMMY_VALUE 1e10

#define ACC_CONST 1e-5

#define ASSUMED_ACCURACY 1e-14
void read_3D_volume(float ***Buf, long Nz, long Ny, long Nx, char *filename)
{
  long j,k;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fread(Buf[k][j],sizeof(float),Nx,fi);
  fclose(fi);
}

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

void read_2D_surface(float **Buf, long Ny, long Nx, char *filename)
{
  long j;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (j=0;j<Ny;j++)
    fread(Buf[j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_2D_surface(float **Buf, long Ny, long Nx, char *filename)
{
  long j;
  FILE *fo;

  fo=file_open(filename,"wb");
  for (j=0;j<Ny;j++)
    fwrite(Buf[j],sizeof(float),Nx,fo);
  fclose(fo);
}

long cln(long j, long i, long Ncol)
{
  return (j-1)*Ncol+i;
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

void minmax(float const *buf, long N1, long N2, float *min, float *max)
{
 long i;

 *min=*max=buf[N1];
 for (i=N1+1;i<=N2;i++)
 {
   if (buf[i]<*min) *min=buf[i];
   if (buf[i]>*max) *max=buf[i];
 }
}

long roun(double x)
{
  long i;
  i = floor(x);
  if (x-i >= 0.5) i++;
  return i;
}




int main(int narg, char **argv)
{
  int N;  
  int NNx,NNy;
  int i,j,ii,jj;
  int iii, jjj;
  int R;
  int Nc;  

  long Nx,Ny,Nz;
  long k_gen, k_gen_v;
  long uk_min, uk_max;  /*Minimum and maximum fine grid k (Z) indices, bounding the upper physical layer*/ 

  float x,y;
  float x0,y0,z0;
  float h;
 
 // float z_node;		/*Depth of the node*/
  float dd;
  //float ddz;		/*Generalized grid depth (Z) discretisation step*/
  float ss;
 // float Tmin, Tmax;	/*Minimum and maximum of the upper layer top position*/
  //float Zmin, Zmax;	/*Minimum and maximum of the interface position*/ 

  float pz; // критерий значимости 

  //float ** Zintf, ** UL_ft; 
 // float **I_dzcm;	/*Interface depth 2D coarse grid cells slowness model (cells from 1,1)*/
 // float *T;		/*Upper layer top position (may be topography)*/
 // float *Z;		/*Reference interface position (obtained by the averaging of fine Z model)*/
  //float *** Buff_top_c; // N rays of Cells_hw_dw_top + Cells_reflected;(NNz,NNx,NNy)
  float **Buff_points_c; // по грубой сетке
  //float *** Rcave;
  float ** points_cave;

  float **Mask_borders;
  float **Borders; 

  char  fi_borders[FLNLN], fi_points_cave[FLNLN];
  char  fi_out_Mask_borders[FLNLN];
  char buf[STRLN];
  char I_parametr[STRLN];
  FILE *fp; 

  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fscanf(fp,"%ld",&k_gen);
  fscanf(fp,"%ld",&k_gen_v);
  fclose(fp);
 
  NNx=(Nx-1)/k_gen;
  NNy=(Ny-1)/k_gen;

  if (((NNx*k_gen)!=(Nx-1)) || ((NNy*k_gen)!=(Ny-1)))
   {
    fprintf(stderr,"Invalid generalization constant.\n");
    exit(EXIT_FAILURE);
   }

  N=NNy*NNx;
  dd=h*k_gen;
  //N_int=(NNy-2)*(NNx-2); /*Number of internal points*/

 
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
  
                                                                
  fscanf(fp,"%s",fi_borders);	fgets(buf,STRLN,fp); /*1th LINE: borders file*/ 
  fscanf(fp,"%s",fi_points_cave); fgets(buf,STRLN,fp); /*2th LINE: Ray caverege- reflected points */


  
  fscanf(fp,"%s",fi_out_Mask_borders);	fgets(buf,STRLN,fp);/*5th LINE: out file */

  fscanf(fp,"%s",I_parametr);	fgets(buf,STRLN,fp);/*6th LINE: interpolation parametr L ,N */
  if (strcmp(I_parametr,"L")==0)  printf("Linear interpolation.\n");
   else
    {
     if (strcmp(I_parametr,"N")==0) 
      {
        fscanf(fp,"%d",&R);	fgets(buf,STRLN,fp); /*7th LINE R for naer neibhour interpolation  */
        if (R<=0) 
         {
          printf("R<=0  R=%d.\n",R);
          exit(1);
         }
        printf("naer neibhour interpolation R=%d.\n",R);
       }
      
     else
    {
      printf("Unrecognized parametr intrpolation %s.\n",I_parametr);
      exit(1);
    }
    }
  
  if(fscanf(fp,"%g",&pz)!=1) pz=0;	fgets(buf,STRLN,fp);/*7 (8) LINE */  
  if (strcmp(argv[1],"-f")==0) fclose(fp);
  /*--------------------------------------------------------*/ 
  

  /*------I STAGE: reading reference model values and averaging them onto coarse grid-------*/
  fprintf(stderr,"Allocating memory...");
 
  //alloc_float_matrix(&Zintf,Ny,Nx);
//  alloc_float_matrix(&UL_ft,Ny,Nx);
 
// I_dzcm=matrix(1,NNy,1,NNx);
  
  fprintf(stderr,"Success.\n");

 

  /*------Read fine grids of layer top and interface*/
  fprintf(stderr,"Read surfaces..."); 
  fprintf(stderr,"Success.\n");

  


  

 alloc_float_matrix(&Buff_points_c,NNy, NNx);//значения  покрытия по грубой сетке

 alloc_float_matrix(&Mask_borders,Ny,Nx);
 alloc_float_matrix(&Borders,Ny,Nx);// маскированная  модель границ

  read_2D_surface(Buff_points_c, NNy, NNx, fi_points_cave);

 /*------Read reference true velocity model-------------------------------------------*/
  fprintf(stderr,"Reading fine Borders model...\n");
 
  read_2D_surface(Borders,  Ny,  Nx, fi_borders);
  fprintf(stderr,"Success...");

  fprintf(stderr,"First interpolate sub-interface slowness variation...");
  //alloc_float_matrix(&SI_dsv_fine,Ny,Nx);
  /*3D-W14: new interpolation procedure with nodes centered in the centers of coarse grid cells*/

  alloc_float_matrix(&points_cave,Ny, Nx);  //интерполированные значения лучевого покрытия(нижний слой)

  if (strcmp(I_parametr,"L")==0) 
   {
    for (j=0;j<Ny;j++)
    {
      jj=floor((h*(float)j)/dd-0.5)+1;
      if (jj<1) jj=1; if (jj>=NNy) jj=NNy-1;
      y=(float)j*h-((float)jj-0.5)*dd;
      for (i=0;i<Nx;i++)
      {
        ii=floor((h*(float)i)/dd-0.5)+1;
        if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;
        x=(float)i*h-((float)ii-0.5)*dd;
        points_cave[j][i]=(1.0/(dd*dd))*
	            (Buff_points_c[jj-1][ii-1]*(dd-x)*(dd-y)+
	             Buff_points_c[jj-1][ii+1-1]*x*(dd-y)+
	             Buff_points_c[jj-1+1][ii-1]*(dd-x)*y+
	             Buff_points_c[jj+1-1][ii+1-1]*x*y);    
      }
    }    
    fprintf(stderr,"Success.\n");
   }

  else 
   {
    for (j=0;j<Ny;j++)
     {
       jj=floor((h*(float)j)/dd-0.5)+1;
       if (jj<1) jj=1; if (jj>=NNy) jj=NNy-1;     
       for (i=0;i<Nx;i++)
       {
        ii=floor((h*(float)i)/dd-0.5)+1;
        if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;       
        if( Buff_points_c[jj-1][ii-1]>0)    points_cave[j][i]=Buff_points_c[jj-1][ii-1];
        else  
        {
         ss=0; Nc=0;
         for( jjj=jj-R; jjj<=jj+R; jjj++)
          for( iii=ii-R; iii<=ii+R; iii++ )
           {
             //fprintf(stderr,"iii=%d jjj=%d.\n", iii, jjj); 
            if (jjj>0 && jjj<NNy && iii>0 && iii<NNx)
             {
              ss+=Buff_points_c[jjj-1][iii-1];
              Nc++;  
              }
            }
           points_cave[j][i]=ss/Nc;  
         } 

      }
    }    
    fprintf(stderr,"Success.\n");
   }

  free_float_matrix(Buff_points_c, NNy, NNx);   

  
   for (j=0; j<Ny; j++)
    for (i=0; i<Nx; i++)
     {
       if( points_cave[j][i]>pz ) Mask_borders[j][i]=Borders[j][i];
       else Mask_borders[j][i]=DUMMY_VALUE ;
      } 
  
  write_2D_surface( Mask_borders , Ny,  Nx, fi_out_Mask_borders);

  free_float_matrix( Mask_borders, Ny, Nx);
  free_float_matrix( Borders, Ny, Nx);
  
  free_float_matrix( points_cave,Ny, Nx);  
  //free_matrix(I_dzcm,1,NNy,1,NNx);
  

 
}
