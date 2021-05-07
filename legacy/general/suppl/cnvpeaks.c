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
  int i, c;
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

  fgets(s,STRLN,fi);
  c=0;
  while (!feof(fi))
  {
    sscanf(s,"%d%g%g%g%*d%*d%*d%g",&i,&xp,&yp,&zp,&tp);
    d=sqrt((xs-xp)*(xs-xp)+(ys-yp)*(ys-yp)+(zs-zp)*(zs-zp));
    az=180.0/M_PI*atan2(xp-xs,yp-ys);
    fprintf(fo,"%5i %5i %10g %10g %10g\n",0,i,d,az,tp);
    c++;
    fgets(s,STRLN,fi);
  }
  printf("%d peaks done.\n",c);

  fclose(fo);
  fclose(fi);
  return 1;
}
