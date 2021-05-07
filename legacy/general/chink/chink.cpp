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


#define FLNLN 130 /*Maximum length of the file name*/
#define STRLN 256

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

int main(int narg, char **argv)
{
  int i,j,k;
  long Nx,Ny,Nz, k_gen, k_gen_v;
  float x0,y0,z0;
  float h; 
  float x,y,z;  
  float ***A;
  
  char  fi_inp[FLNLN], fi_out[FLNLN];
  char buf[STRLN];
  FILE* fp;

  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fscanf(fp,"%ld",&k_gen);
  fscanf(fp,"%ld",&k_gen_v);
  fclose(fp);
  
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
  
                                                                
  fscanf(fp,"%s",fi_inp);	fgets(buf,STRLN,fp); 
  fscanf(fp,"%s",fi_out);	fgets(buf,STRLN,fp); 
  fscanf(fp,"%g",&x); fscanf(fp,"%g",&y); fgets(buf,STRLN,fp);
  if (strcmp(argv[1],"-f")==0) fclose(fp);
  
  
  alloc_float_3d(&A, Nz, Ny, Nx);//значения лучевого покрытия по грубой сетке

  read_3D_volume(A, Nz, Ny,  Nx, fi_inp);
  fp=file_open(fi_out,"wt");   
  
  
  i=floor((x-x0)/h);
  j=floor((y-y0)/h);
  
  fprintf(stderr,"%g %g %d %d\n",x,y,i,j);
   
  for(k=0; k<Nz; k++)
  {
   z=z0+h*k; 
   fprintf(fp,"%10g  %10g \n",A[k][j][i],-z);
  }
 fclose(fp);
}
