#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "subio.h"
#include "submem.h"

#define STRLN 256
#define LN 160		/*Maximum length of string in file*/



void cross_correlation (float *x, float *y, float*result, float maxdelay, int n )
{

  int i,j,delay;
  double mx,my,sx,sy,sxy,denom;
   
   /* Calculate the mean of the two series x[], y[] */
   mx = 0;
   my = 0;   
   for (i=0;i<n;i++) {
      mx += x[i];
      my += y[i];
   }
   mx /= n;
   my /= n;

   /* Calculate the denominator */
   sx = 0;
   sy = 0;
   for (i=0;i<n;i++) {
      sx += (x[i] - mx) * (x[i] - mx);
      sy += (y[i] - my) * (y[i] - my);
   }
   denom = sqrt(sx*sy);

   /* Calculate the correlation series */
   for (delay=-maxdelay;delay<maxdelay;delay++) {
      sxy = 0;
      for (i=0;i<n;i++) {
         j = i + delay;
         if (j < 0 || j >= n)
            continue;
         else
            sxy += (x[i] - mx) * (y[j] - my);
         /* Or should it be (?)
         if (j < 0 || j >= n)
            sxy += (x[i] - mx) * (-my);
         else
            sxy += (x[i] - mx) * (y[j] - my);
         */
      }
      result[delay] = sxy / denom;
      //fprintf(stdout,"%g\n", r);
      /* r is the correlation coefficient at "delay" */

   }
}


int main (int narg, char **argv)
{
  int NR, i;
  char string[LN];
  float *s1_t, *s1_a, *s2_t, *s2_a, *result;

  FILE *fd1, *fd2;
  
  if (narg<3)
  {
    printf("Insufficient parameters.\nUSAGE: cross_corr <signal 1 data file> <signal 2 data file>\n");
    exit(1);
  }
  
  
  
 // parameters reading
 
 
  fd1=file_open(argv[1],"rt");
  fd2=file_open(argv[2],"rt");

  i=0;
  while(!feof(fd1))
  {
    fgets(string,LN,fd1);
    i++;
  }
  
  NR=i;
  
  fprintf(stdout,"Number of strings in input file is %d\n", NR);
  
  rewind(fd1);
   
  s1_t=(float*)fmemalloc((NR),sizeof(float));
  s1_a=(float*)fmemalloc((NR),sizeof(float));
  s2_t=(float*)fmemalloc((NR),sizeof(float));
  s2_a=(float*)fmemalloc((NR),sizeof(float));
  result=(float*)fmemalloc((2*NR),sizeof(float));
  
  for(i=0;i<=NR;i++)
  {
    fgets(string,LN,fd1);
    sscanf(string,"%g%g",&s1_t[i], &s1_a[i]); 
    
    fgets(string,LN,fd2);
    sscanf(string,"%g%g",&s2_t[i], &s2_a[i]); 
  }  
  fclose(fd1);
  fclose(fd2);
 
 
 // cross correlation
  
  cross_correlation(s1_a, s2_a, result, NR, NR);
 
 // length of resulting signal adjustment
 
 // printing result to std out
 
 for(i=0;i<=2*NR;i++) fprintf(stdout,"%g\n", result[i]);
  
 // cleaning
  fmemfree(s1_t);
  fmemfree(s1_a);
  fmemfree(s2_t);
  fmemfree(s2_a);
  fmemfree(result);
  
  return 1;
    
}
