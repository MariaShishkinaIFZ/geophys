#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "subio.h"
#include "submem.h"
#include "nr.h"
#include "nrutil.h"


int main(int narg, char **argv)
{
  float sigma_X, sigma_Y, sigma_Z, sigma_T, X, Y, Z, s_RMS;
  int N_string, i,N_rec;
  
  long idum=(-13);
  long long src_name;
  long double T;
  
  FILE *fi, *fo;
  
  if (narg<7)
  {
    fprintf(stderr,"Insufficient parameters.\n rand_cur_ev <current_events.dat> sigma_X sigma_Y sigma_Z sigma_T <randomised_current_events.dat>\n");
    exit(EXIT_FAILURE);
  }
  
  fi=file_open(argv[1],"rt");
   
  sscanf(argv[2],"%f",&sigma_X);
  sscanf(argv[3],"%f",&sigma_Y);
  sscanf(argv[4],"%f",&sigma_Z);
  sscanf(argv[5],"%f",&sigma_T);
  
  fo=file_open(argv[6],"wt");
  
  fscanf (fi,"%d",&N_string);
  fprintf(fo,"%d\n",N_string);
  
  idum=-rand();
  
  for(i=0; i<N_string; i++)
   {
     fscanf (fi,"%lld %f %f %f %Lf %f %d",&src_name,&X,&Y,&Z,&T,&s_RMS,&N_rec);
     
     X+=sigma_X*gasdev(&idum);
     Y+=sigma_Y*gasdev(&idum);
     Z+=sigma_Z*gasdev(&idum);
     T+=(long double)sigma_T*gasdev(&idum);
     
     fprintf (fo,"%lld %10.2f %10.2f %10.2f %18.2Lf %10.2f %5d\n",src_name,X,Y,Z,T,s_RMS,N_rec);
   }
  
fclose(fi);
fclose(fo);  
  
}