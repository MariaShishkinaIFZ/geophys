
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
  float right_boarder, left_boarder;
  float right_cut, left_cut;
  float xp, yp, zp, tp;
  float d, az;
  char s[STRLN];

  FILE *fi, *fo_dw, *fo_hw;

  if (narg<7)
  {
    printf("Insufficient parameters.\nUSAGE: cnvpeaks <peaks data file> <dw out file> <hw out file> ys a b c d \nwhere a = left boarder of diving wave and b = right boarder of diving wave \n");
    exit(1);
  }

  fi=file_open(argv[1],"rt");
  fo_dw=file_open(argv[2],"wt");
  fo_hw=file_open(argv[3],"wt");

  sscanf(argv[4],"%g",&ys);
  sscanf(argv[5],"%g",&left_boarder);
  sscanf(argv[6],"%g",&right_boarder);
  sscanf(argv[7],"%g",&left_cut);
  sscanf(argv[8],"%g",&right_cut);


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


   if (((c%1)==0) && ((yp<left_boarder) || (right_boarder<yp)) && ((left_cut<yp) && (yp<right_cut))) fprintf(fo_hw,"1   %5i %10g %10g %10g\n",i,d,az,tp);

   if (((c%1)==0) && ((left_boarder<yp) && (yp<right_boarder)) && ((left_cut<yp) && (yp<right_cut))) fprintf(fo_dw,"1   %5i %10g %10g %10g\n",i,d,az,tp);


    c++;
    fgets(s,STRLN,fi);
  }
  printf("%d peaks done.\n",c);

  fclose(fo_hw);
  fclose(fo_dw);
  fclose(fi);
  return 1;
}
