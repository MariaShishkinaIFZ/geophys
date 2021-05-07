
//convertor desined for Vladov data
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

  if (narg<4)
  {
    printf("Insufficient parameters.\nUSAGE: cnvpeaks <peaks data file> <out file> xs ys zs\n");
    exit(1);
  }

  fi=file_open(argv[1],"rt");
  fo=file_open(argv[2],"wt");

  sscanf(argv[3],"%g",&ys);


xs=0;
zs=0;

  fgets(s,STRLN,fi);
  c=0;
  while (!feof(fi))
  {
    sscanf(s,"%g%g",&yp,&tp);

    i=(int)yp;

    xp=0;
    yp=yp/1000;
    zp=0;

    tp=tp/1000;

    d=sqrt((xs-xp)*(xs-xp)+(ys-yp)*(ys-yp)+(zs-zp)*(zs-zp));
    az=180.0/M_PI*atan2(xp-xs,yp-ys);
   if ((c%2)==0) fprintf(fo,"1   %5i %10g %10g %10g\n",i,d,az,tp);
    c++;
    fgets(s,STRLN,fi);
  }
  printf("%d peaks done.\n",c);

  fclose(fo);
  fclose(fi);
  return 1;
}
