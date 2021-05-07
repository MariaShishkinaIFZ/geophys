#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "subio.h"
#include "submem.h"
#include "nr.h"
#include "nrutil.h"


int main(int narg, char **argv)
{
  float sigma_T;
  int N_string, i,N_rec;
  
  long idum=(-13);
  
  int id_p, prof_p, prof_s, tpk_number;
  float d, az, tp;
  char name_s[100];
  
  FILE *fi, *fo;
  
  if (narg<4)
  {
    fprintf(stderr,"Insufficient parameters.\n rand_tpk <picks.tpk>  sigma_T <randomised_picks.tpk>\n");
    exit(EXIT_FAILURE);
  }
  
  fi=file_open(argv[1],"rt");
   

  sscanf(argv[2],"%f",&sigma_T);
  
  fo=file_open(argv[3],"wt");
  
  srand(time(0));
 
  idum=-rand();
  

  while (!feof(fi))
   {
     tpk_number=fscanf(fi,"%d %s %d  %d %g %g %g",&prof_s,&name_s,&prof_p,&id_p,&d,&az,&tp);
     if (tpk_number==7) 
     {
      tp+=(float)sigma_T*gasdev(&idum);
      fprintf(fo,"%5d %s %5d  %5d %10g %10g %10g\n",prof_s,name_s,prof_p,id_p,d,az,tp);
     }
   }
  
fclose(fi);
fclose(fo);  
return 0;
  
}