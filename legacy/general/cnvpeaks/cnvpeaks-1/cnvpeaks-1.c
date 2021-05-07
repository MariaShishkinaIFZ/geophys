#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "subio.h"

#define STRLN 256

int gras_round(float x)
{
  return x-floor(x) >= 0.5 ? ceil(x) : floor(x);
}

int main(int narg, char **argv)
{
  long long i;
  int nprof, c, Np;
  float xs,ys,zs;
  float xp, yp, zp, tp;
  float d, az;
  char s[STRLN];

  FILE *fi, *fo;

  if (narg<6)
  {
    printf("Insufficient parameters.\nUSAGE: cnvpeaks <peaks data file> <out file> xs ys zs\n");
    exit(1);
  }

  fi=file_open(argv[1],"rt");
  fo=file_open(argv[2],"wt");
  sscanf(argv[3],"%g",&xs);
  sscanf(argv[4],"%g",&ys);
  sscanf(argv[5],"%g",&zs);

  c=0;
  fscanf(fi,"%d",&Np);
  printf("Peaks data file consists of %d peaks.\n",Np);
  fgets(s,STRLN,fi);
  fgets(s,STRLN,fi);

  for (c=0;c<Np;c++)
  {
    sscanf(s,"%d%lld%g%g%g%*d%*d%*d%g",&nprof,&i,&xp,&yp,&zp,&tp);
    d=sqrt((xs-xp)*(xs-xp)+(ys-yp)*(ys-yp)+(zs-zp)*(zs-zp));
    az=180.0/M_PI*atan2(xp-xs,yp-ys);
    fprintf(fo,"%5d %10lld %10g %10g %10g\n",nprof,i,d,az,tp);
    fgets(s,STRLN,fi);
  }
  printf("%d peaks done.\n",c);

  fclose(fo);
  fclose(fi);
  return 1;
}
