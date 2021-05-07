/*Making sclices from the 3d file*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "submem.h"
#include "subio.h"

#define GEOMFILE "model-geometry.gras"
#define STRLN 256

int gras_round(float x)
{
  
  return (x-floor(x)) >= 0.5 ? ceil(x) : floor(x);
}

int main(int narg, char **argv)
{
  long Np, i;
  float Xa, Ya, Xb, Yb, Yp, dist;
  float BW;
  float Lprof,tg_a,sin_a,b,b1;
  float x_c, y_c, x_prf;
  float *X, *Y, *Z;
  char buf[STRLN];
  char fln_in[STRLN], fln_o[STRLN];
  
  FILE *fp, *fi, *fo;
  
  if (narg < 2) errquit("Insufficient parameters\nUSAGE: sl3d -f parfilename | -s\n");
  if (strcmp(argv[1],"-f")==0)
  {
    if (narg < 3) errquit("Insufficient parameters\nUSAGE: sl3d -f parfilename | -s\n");
    fp=file_open(argv[2],"rt");
  }
  else
    if (strcmp(argv[1],"-s")!=0) errquit("Unknown option\nUSAGE: sl3d -f parfilename | -s\n");
    else fp=stdin;
    
  fscanf(fp,"%s",fln_in); fgets(buf,STRLN,fp);
  fscanf(fp,"%s",fln_o); fgets(buf,STRLN,fp);
  fscanf(fp,"%g%g",&Xa,&Ya);
  fscanf(fp,"%g%g",&Xb,&Yb);
  fscanf(fp,"%g",&BW);
  fgets(buf,STRLN,fp);
  if (fp!=stdin) fclose(fp);

  fi=file_open(fln_in,"rt");
  Np=0;
  do
  {
    if (fgets(buf,STRLN,fi)!=NULL) Np++;
    else break;
  }
  while (!feof(fi));
  
  fprintf(stderr,"%ld events in file.\n",Np);
  
  X=fmemalloc(Np,sizeof(float));
  Y=fmemalloc(Np,sizeof(float));
  Z=fmemalloc(Np,sizeof(float));
  
  rewind(fi);
  for (i=0;i<Np;i++)
    if (fscanf(fi,"%g%g%g",X+i,Y+i,Z+i)!=3)
      errquit("Error reading data file\n");
  fclose(fi);
  fprintf(stderr,"Events read.\n");
  
  fo=file_open(fln_o,"wt");
  if (fabs(Xa-Xb)<1e-13)
  {
    for (i=0;i<Np;i++)
      if (fabs(X[i]-Xa)<=BW)
	fprintf(fo,"%13g %13g %13g\n",X[i],Y[i],Z[i]);
  }
  else {
    if (fabs(Ya-Yb)<1e-13)
    {
      for (i=0;i<Np;i++)
        if (fabs(Y[i]-Ya)<=BW)
	  fprintf(fo,"%13g %13g %13g\n",X[i],Y[i],Z[i]);
    }
    else {
      tg_a=(Ya-Yb)/(Xa-Xb);
      b=Ya-tg_a*Xa;
      sin_a=((Ya-Yb)/pow((Xa-Xb)*(Xa-Xb)+(Ya-Yb)*(Ya-Yb),0.5));
      fprintf(stderr,"tg_a: %g b: %g sin_a: %g\n",tg_a,b,sin_a);
      for (i=0;i<Np;i++)
      {
	Yp=X[i]*tg_a+b;
	dist=fabs((Yp-Y[i])/sin_a);
	if (dist<=BW)
	{
	  b1=Y[i]+tg_a*X[i];
          x_c=(b1-b)/(2*tg_a);
          y_c=tg_a*x_c+b;
	  x_prf=pow((x_c-Xa)*(x_c-Xa)+(y_c-Ya)*(y_c-Ya),0.5);
	  fprintf(fo,"%13g %13g\n",x_prf,-Z[i]);
	}
      }
    }
  }
  fclose(fo);
  
  fmemfree(Z);
  fmemfree(Y);
  fmemfree(X);
  
  return 0;
}
