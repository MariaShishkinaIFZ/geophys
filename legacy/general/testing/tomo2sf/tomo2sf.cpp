//#include <math.h>
//#include <exception>
#include <iostream>


#include "gstdata.h"
#include "gsterror.h"
#include "gstdefs.h"
//#include "gstmath.h"
#include "gstio.h"
#include "subio.h"
//#include "lsqr.h"
#include <stdlib.h>
//#define REL_MATR_ERR 1E-4
//#define MAX_CONDITION 1E7
#define STRLN 256


int main(int narg, char ** argv)
{
  int i,j;
 // int count;
  int kx,kz;
  int Mph,Nph;
  int Mvp,Nvp;
  int Mvs,Nvs;
  
  int Nx1,Nx2;
  int Nz1,Nz2;
  
  char f_pho[GST_STRLN];
  char f_Vp[GST_STRLN]; 
  char f_Vs[GST_STRLN];  
  char f_out[GST_STRLN];
  char buf[STRLN];
  
  GSTgriddedData *pho;
  GSTgriddedData *Vp;
  GSTgriddedData *Vs;
  
  FILE *fp; 
  
  
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
  

  fscanf(fp,"%s%s%s",f_pho,f_Vp,f_Vs); 	fgets(buf,STRLN,fp); /*1th LINE: Rho.gmt Vp.gmt Vs.gmt*/
  fscanf(fp,"%d%d",&kx, &kz);  		fgets(buf,STRLN,fp); /*2nd LINE: */
  fscanf(fp,"%s",f_out); 		fgets(buf,STRLN,fp); /*3rd LINE: Text out file for spec2d*/
  
  
  if (strcmp(argv[1],"-f")==0) fclose(fp);
  /*--------------------------------------------------------*/ 
  
  
  
 //if (narg<2)  GSTerror("Insufficient parameters \n");
 //fp=fopen(argv[1],"rt");  
 //fscanf(fp,"%s%s%s",f_pho,f_Vp,f_Vs); 
 //fscanf(fp,"%d%d",&kx, &kz);  
 //fscanf(fp,"%s",f_out); 
 //fclose(fp); 
 //fprintf(stderr,"read file par\n");
 
  pho=readGMT(f_pho);
  Vp=readGMT(f_Vp);
  Vs=readGMT(f_Vs);
  
  pho->size(Mph,Nph);
  Vp->size(Mvp,Nvp);
  Vs->size(Mvs,Nvs);
   if (Mph!=Mvp || Mvp!=Mvs || Nph!=Nvp || Nvp!=Nvs)  GSTerror("Grids don't coincide \n");
  fp=fopen(f_out,"wt");
  
  fprintf(fp,"nbmodels                        = %d              # nb of different models\n",Mph*Nph);
  
  for (i=0; i<Mph; i++)
   {
     for (j=0; j<Nph; j++)
     {
      fprintf(fp,"%d 1 %10g %10g %10g 0 0  9999 9999 0 0 0 0 0 0\n",(i)*Nph+j+1,pho->value(i,j),Vp->value(i,j), Vs->value(i,j));
     }
   }
  
  fprintf(fp,"nbregions                        = %d              #nb of regions and model number for each\n",Mph*Nph);
  
  for (i=0; i<Mph; i++)
   {
     for (j=0; j<Nph; j++)
     {
       Nx1=(j)*kx+1;
       Nx2=(j+1)*kx;
       Nz1=(i)*kz+1;
       Nz2=(i+1)*kz;
       
      fprintf(fp,"%10d %10d %10d %10d %10d \n",Nx1,Nx2,Nz1,Nz2, (i)*Nph+j+1);
     }
   }
  fclose(fp);
   
}
