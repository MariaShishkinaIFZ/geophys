 
/*
-читаю модель и  2D 3D покрытие
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

#define FLNLN 256 	/*Maximum length of the file name*/
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

int eq_digital(float x1, float x2, float eps)
{
  if (fabs(x1-x2)<eps) return 1;
  return 0;
}

int main(int narg, char **argv)
{
  int N;  
  int NNx,NNy, NNz;
  int i,j,k,ii,jj,kk;
  int iii, jjj, kkk;
  int R;
  int Nc;  

  long Nx,Ny,Nz;
  long k_gen, k_gen_v;
  long uk_min, uk_max;  /*Minimum and maximum fine grid k (Z) indices, bounding the upper physical layer*/ 

  float x,y,z;
  float x0,y0,z0;
  float h;
 
  float z_node;		/*Depth of the node*/
  float dd;
  float ddz;		/*Generalized grid depth (Z) discretisation step*/
  float ss;
  float Tmin, Tmax;	/*Minimum and maximum of the upper layer top position*/
  float Zmin, Zmax;	/*Minimum and maximum of the interface position*/ 

  float pz; // критерий значимости 

  float ** Zintf, ** UL_ft; 
  float **I_dzcm;	/*Interface depth 2D coarse grid cells slowness model (cells from 1,1)*/
  float *T;		/*Upper layer top position (may be topography)*/
  float *Z;		/*Reference interface position (obtained by the averaging of fine Z model)*/
  float *** Buff_top_c; // N rays of Cells_hw_dw_top + Cells_reflected;(NNz,NNx,NNy)
  float **Buff_2D_bot_c;
  float *** Rcave;
  float ** bot2D_Rcave;

  float ***Mask_velmod;
  float ***Velmod; 

  char  fi_fiTop[FLNLN], fi_fiBot[FLNLN], fi_VelMod[FLNLN], fi_Rcave[FLNLN],fi_bot2D_Rcave[FLNLN],
        fi_out_coverage_fine_file_2D[FLNLN], fi_out_coverage_fine_file_3D[FLNLN];
  char  fi_out_Mask_velmod[FLNLN];
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
  if (narg < 2) errquit("Insufficient parameters.\n		USE: 3D_mask <-f|s> \n ");
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
  
                                                                
  fscanf(fp,"%s",fi_VelMod);	fgets(buf,STRLN,fp); /*1th LINE: Fine true reference velocity model*/
  fscanf(fp,"%s",fi_Rcave);	fgets(buf,STRLN,fp); /*2th LINE:Ray caverege for top layer(3D) */
  fscanf(fp,"%s",fi_bot2D_Rcave); fgets(buf,STRLN,fp); /*3th LINE: Ray caverege for botom layr*/


  fscanf(fp,"%s",fi_fiBot);	fgets(buf,STRLN,fp); /*4th LINE: Fine reference interface position*/
  fscanf(fp,"%s",fi_fiTop);	fgets(buf,STRLN,fp); /*5th LINE: Fine layer top position*/
  fscanf(fp,"%s",fi_out_Mask_velmod);	fgets(buf,STRLN,fp);/*6th LINE: out file */
  fscanf(fp,"%s",fi_out_coverage_fine_file_3D);	fgets(buf,STRLN,fp); /*7th LINE: out coverage_fine_file_2D */
  fscanf(fp,"%s",fi_out_coverage_fine_file_2D);	fgets(buf,STRLN,fp); /*8th LINE: out coverage_fine_file_3D */

  puts(fi_VelMod);
  puts(fi_Rcave);
  puts(fi_bot2D_Rcave);
  puts(fi_fiBot);
  puts(fi_fiTop);
  puts(fi_out_Mask_velmod);

  fscanf(fp,"%s",I_parametr);/*9th LINE: interpolation parametr -L ,-N */
  if (strcmp(I_parametr,"L")==0)  printf("Linear interpolation.\n");
   else
    {
     if (strcmp(I_parametr,"N")==0) 
      {
        fscanf(fp,"%d",&R);
        if (R<=0) 
         {
          printf("R<=0  R=%d.\n",R);
          //exit(1);
         }
        printf("naer neibhour interpolation R=%d.\n",R);
       }
      
     else
    {
      printf("Unrecognized parametr intrpolation %s.\n",I_parametr);
      exit(1);
    }
    }
  fgets(buf,STRLN,fp); 
  
  if(fscanf(fp,"%g",&pz)!=1) pz=0;	fgets(buf,STRLN,fp);/*10 LINE */  
  if (strcmp(argv[1],"-f")==0) fclose(fp);
  /*--------------------------------------------------------*/ 
  

  /*------I STAGE: reading reference model values and averaging them onto coarse grid-------*/
  fprintf(stderr,"Allocating memory...");
 // alloc_float_3d(&VFINE,Nz,Ny,Nx);
  //alloc_float_3d(&dSFINE,Nz,Ny,Nx);
  alloc_float_matrix(&Zintf,Ny,Nx);
  alloc_float_matrix(&UL_ft,Ny,Nx);
  //SI_dscm=matrix(1,NNy,1,NNx);
  I_dzcm=matrix(1,NNy,1,NNx);
  T=fvector(1,N);
  Z=fvector(1,N);
  fprintf(stderr,"Success.\n");

 

  /*------Read fine grids of layer top and interface*/
  fprintf(stderr,"Read surfaces...");
  read_2D_surface(UL_ft,Ny,Nx,fi_fiTop);
  read_2D_surface(Zintf,Ny,Nx,fi_fiBot);
  fprintf(stderr,"Success.\n");

  /*Averaging*/
  av2d(UL_ft,I_dzcm,Ny,Nx,k_gen);
  for (jj=1;jj<=NNy;jj++)
    for (ii=1;ii<=NNx;ii++)
      T[cln(jj,ii,NNx)]=I_dzcm[jj][ii]; /*From here now T contains topografy*/
  av2d(Zintf,I_dzcm,Ny,Nx,k_gen);
  for (jj=1;jj<=NNy;jj++)
    for (ii=1;ii<=NNx;ii++)
      Z[cln(jj,ii,NNx)]=I_dzcm[jj][ii]; /*From here now Z and I_dzcm contains TOTAL Z of interface*/
  minmax(T,1,N,&Tmin,&Tmax);
  fprintf(stderr,"Layer top position read. Tmin=%g Tmax=%g\n",
  		Tmin,Tmax);
  minmax(Z,1,N,&Zmin,&Zmax);
  fprintf(stderr,"Reference interface position read. Zmin=%g Zmax=%g\n",
  		Zmin,Zmax);
  
  if (eq_digital((Tmin-z0)/h,roun((Tmin-z0)/h),1e-7*h)) {fprintf(stderr,"%g %ld ",(Tmin-z0)/h,roun((Tmin-z0)/h)); uk_min=roun((Tmin-z0)/h)+1;}
  else uk_min=floor((Tmin-z0)/h)+1; fprintf(stderr,"uk_min=%ld\n",uk_min);
  if (eq_digital((Zmax-z0)/h,roun((Zmax-z0)/h),1e-7*h)) {fprintf(stderr,"%g %ld ",(Zmax-z0)/h,roun((Zmax-z0)/h)); uk_max=roun((Zmax-z0)/h);}
  else uk_max=floor((Zmax-z0)/h); fprintf(stderr,"uk_max=%ld\n",uk_max);
  
  
  NNz=(uk_max-uk_min)/k_gen_v; if ((uk_max>uk_min) && (NNz<1)) NNz=1;
  //fprintf(stderr,"NNz=%ld\n",NNz);
  if ((uk_min+NNz*k_gen_v)<uk_max) {NNz++; uk_max=uk_min+NNz*k_gen_v; fprintf(stderr,"New is uk_max=%ld\n",uk_max);}
  ddz=h*k_gen_v;


  if (((Tmin-z0)/h - roun((Tmin-z0)/h)) < 1e-7*h) uk_min=floor((Tmin-z0)/h)+2;
  else uk_min=floor((Tmin-z0)/h)+1;
  if (((Zmax-z0)/h - roun((Zmax-z0)/h)) < 1e-7*h) uk_max=floor((Zmax-z0)/h)+2;
  else uk_max=floor((Zmax-z0)/h)+2;


 alloc_float_3d(&Buff_top_c, NNz, NNy, NNx);//значения лучевого покрытия по грубой сетке
 alloc_float_matrix(&Buff_2D_bot_c,NNy, NNx);

 
 


 alloc_float_3d(&Mask_velmod,Nz,Ny,Nx);// маскированная скоростная модель
 alloc_float_3d(&Velmod,Nz,Ny,Nx);


// read_3D_volume(float ***Buf, long Nz, long Ny, long Nx, char *filename)
  read_3D_volume(Buff_top_c, NNz, NNy,  NNx, fi_Rcave);
  read_2D_surface(Buff_2D_bot_c, NNy, NNx, fi_bot2D_Rcave);

 /*------Read reference true velocity model-------------------------------------------*/
  fprintf(stderr,"Reading fine velocity model...\n");
 // read_3D_volume(VFINE,Nz,Ny,Nx,fi_VelMod);  
  read_3D_volume(Velmod, Nz, Ny,  Nx, fi_VelMod);
 fprintf(stderr,"Success...");

 fprintf(stderr,"First interpolate sub-interface slowness variation...");
  //alloc_float_matrix(&SI_dsv_fine,Ny,Nx);
  /*3D-W14: new interpolation procedure with nodes centered in the centers of coarse grid cells*/

  alloc_float_matrix(&bot2D_Rcave,Ny, Nx);  //интерполированные значения лучевого покрытия(нижний слой)

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
        bot2D_Rcave[j][i]=(1.0/(dd*dd))*
	            (Buff_2D_bot_c[jj-1][ii-1]*(dd-x)*(dd-y)+
	             Buff_2D_bot_c[jj-1][ii+1-1]*x*(dd-y)+
	             Buff_2D_bot_c[jj-1+1][ii-1]*(dd-x)*y+
	             Buff_2D_bot_c[jj+1-1][ii+1-1]*x*y);    
      }
    }    
    fprintf(stderr,"Success.\n");
   }

  else 
   {
    for (j=0;j<Ny;j++)
     {
       //jj=floor((h*(float)j)/dd-0.5)+1;
       jj=j/k_gen+1;
       if (jj<1) jj=1; if (jj>=NNy) jj=NNy-1;     
       for (i=0;i<Nx;i++)
       {
        //ii=floor((h*(float)i)/dd-0.5)+1;
        ii=i/k_gen+1;
        if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;
	

	if (jj==0) jj=1;
	if( Buff_2D_bot_c[jj-1][ii-1]>0)   
	{
	  bot2D_Rcave[j][i]=Buff_2D_bot_c[jj-1][ii-1];
	}
        else  
        {
         if (R>0)
	 {
	 ss=0; Nc=0;
         for( jjj=jj-R; jjj<=jj+R; jjj++)
          for( iii=ii-R; iii<=ii+R; iii++ )
           {
             //fprintf(stderr,"iii=%d jjj=%d.\n", iii, jjj); 
            if (jjj>0 && jjj<=NNy && iii>0 && iii<=NNx)
             {
              ss+=Buff_2D_bot_c[jjj-1][iii-1];
              Nc++;  
              }
            }
          bot2D_Rcave[j][i]=ss/Nc;  
         }
         else
           ss=0.0;
	}
      }
    }    
    fprintf(stderr,"Success.\n");
   }

  write_2D_surface(bot2D_Rcave,Ny,Nx,fi_out_coverage_fine_file_2D);
  
  free_float_matrix(Buff_2D_bot_c, NNy, NNx);
  alloc_float_3d(&Rcave, Nz, Ny, Nx);  //интерполированные значения лучевого покрытия (верхний слой дополненный нижним)



  fprintf(stderr,"Now interpolate UPL slowness variations...");
  if (strcmp(I_parametr,"L")==0) 
   {
    for (k=0;k<Nz;k++)
     {
      /*kk=(k-(uk_min-1))/k_gen_v; if (kk==NNz) kk--;*/
      kk=floor(h*(float)(k-uk_min+1)/ddz-0.5)+1;
      if (kk<1) kk=1; if (kk>=NNz) kk=NNz-1;
      z_node=z0+h*(float)k;
      z=h*(float)(k-uk_min+1)-((float)kk-0.5)*ddz;
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
	    if ((k>floor((UL_ft[j][i]-z0)/h)) && (k<floor((Zintf[j][i]-z0)/h))) /*Node inside reference UPL*/
	     {
               
	       ss= (1.0/(dd*dd*ddz))*
	        (Buff_top_c[kk-1][jj-1][ii-1]*(dd-x)*(dd-y)*(ddz-z)+
	        Buff_top_c[kk-1][jj-1][ii+1-1]*x*(dd-y)*(ddz-z)+
	        Buff_top_c[kk-1][jj+1-1][ii-1]*(dd-x)*y*(ddz-z)+
	        Buff_top_c[kk+1-1][jj-1][ii-1]*(dd-x)*(dd-y)*z+
	        Buff_top_c[kk+1-1][jj-1][ii+1-1]*x*(dd-y)*z+
	        Buff_top_c[kk+1-1][jj+1-1][ii-1]*(dd-x)*y*z+
	        Buff_top_c[kk-1][jj+1-1][ii+1-1]*x*y*(ddz-z)+
	        Buff_top_c[kk+1-1][jj+1-1][ii+1-1]*x*y*z);
	     }
	    else /*Node outside reference UPL*/
	     {
	        if ((k<=floor((UL_ft[j][i]-z0)/h))) ss=-1.0;
	        else /*Node below the interface*/
	         ss=bot2D_Rcave[j][i];
	     }
	      Rcave[k][j][i]=ss;
	  //  VFINE[k][j][i]=VFINE[k][j][i]/(1.0+VFINE[k][j][i]*ss);
        }
       }
     }
    fprintf(stderr,"Success.\n");
   }
  else
   {
    for (k=0;k<Nz;k++)
     {
      
       //kk=floor(h*(float)(k-uk_min+1)/ddz-0.5)+1;
       kk=(k-uk_min)/k_gen_v+1;
       if (kk>=NNz) kk=NNz-1; if (kk<1) kk=1;
       //fprintf(stderr,"%d\n",kk);
       //z_node=z0+h*(float)k;
       //z=h*(float)(k-uk_min+1)-((float)kk-0.5)*ddz;
       for (j=0;j<Ny;j++)
        {
         //jj=floor((h*(float)j)/dd-0.5)+1;
         jj=j/k_gen+1;
         if (jj>=NNy) jj=NNy-1; if (jj<1) jj=1;
         // y=(float)j*h-((float)jj-0.5)*dd;
         for (i=0;i<Nx;i++)
          {
           //ii=floor((h*(float)i)/dd-0.5)+1;
           ii=i/k_gen+1;
           if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;
            //x=(float)i*h-((float)ii-0.5)*dd;
	    if ((k>floor((UL_ft[j][i]-z0)/h)) && (k<floor((Zintf[j][i]-z0)/h))) /*Node inside reference UPL*/
	     {
               //fprintf(stderr,"%ld %ld %ld\n",kk,jj,ii);
	       if(Buff_top_c[kk-1][jj-1][ii-1]>0)
	       {
		   ss=Buff_top_c[kk-1][jj-1][ii-1];
	       }
               else
                {
                  if (R>0)
		  {
		   ss=0; Nc=0;
                   for(kkk=kk-R; kkk<=kk+R; kkk++) 
                    for( jjj=jj-R; jjj<=jj+R; jjj++)
                     for( iii=ii-R; iii<=ii+R; iii++ )
                     {
                      if (jjj>0 && jjj<NNy && iii>0 && iii<NNx &&  (kkk>0) && (kkk<=NNz))
                       {
			 //fprintf(stderr,"%d %d %d\n",kkk-1,jjj-1,iii-1);
			 ss+=Buff_top_c[kkk-1][jjj-1][iii-1];
                         Nc++;  
                        }
                     } 
                  ss=ss/Nc;  
                 }
                 else
		   ss=0.0;
		}
	     }
	    else /*Node outside reference UPL*/
	     {
	        if ((k<=floor((UL_ft[j][i]-z0)/h))) ss=-1.0;
	        else /*Node below the interface*/
	         ss=bot2D_Rcave[j][i];
	     }
	      Rcave[k][j][i]=ss;
	  //  VFINE[k][j][i]=VFINE[k][j][i]/(1.0+VFINE[k][j][i]*ss);
        }
       }
     }
    fprintf(stderr,"Success.\n");
   }



  free_float_3d( Buff_top_c, NNz, NNy, NNx); 
  free_float_matrix( bot2D_Rcave,Ny, Nx);



  for (k=0; k<Nz; k++)
   for (j=0; j<Ny; j++)
    for (i=0; i<Nx; i++)
     {
       if (Rcave[k][j][i]>pz || Rcave[k][j][i]==-1.0 ) Mask_velmod[k][j][i]=Velmod[k][j][i];
       else  Mask_velmod[k][j][i]=DUMMY_VALUE ;
       if (Rcave[k][j][i]==-1.0) Rcave[k][j][i]=0;
      } 



  write_3D_volume(Rcave,Nz,Ny,Nx,fi_out_coverage_fine_file_3D);
  free_float_3d(Rcave,Nz,Ny,Nx);
 //sprintf(buf,"top_layer_ray_coverage.3dc");
  write_3D_volume( Mask_velmod, Nz, Ny,  Nx, fi_out_Mask_velmod);

  free_float_3d( Mask_velmod, Nz, Ny, Nx);
  free_float_3d( Velmod, Nz, Ny, Nx);


  free_float_matrix(Zintf,Ny,Nx);
  free_float_matrix(UL_ft,Ny,Nx);
  free_matrix(I_dzcm,1,NNy,1,NNx);
  free_fvector(T,1,N);
   free_fvector(Z,1,N);

return 0;
 
}
