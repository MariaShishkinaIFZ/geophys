#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "iogras.h"
#include "submem.h"
#include "subio.h"

#define FLNLN 80 	/*Maximum length of the file name*/


int main(int narg, char **argv)
{
  int  l, i;
  int ask;
  float x0, y0, z0, h;
  float xs, ys, zs;
  float nx, ny, nz;
  float V0, gradV;
  float c0, beta, td;
  float **Int;
  char src_file[FLNLN], rec_file[FLNLN], out_file[FLNLN];
  int Np;
  float *xp, *yp, *zp ,*t;
  int *np, *ii, *jj, *kk;
  FILE *fp, *fi;
 
 if (narg < 2) errquit("Insufficient parameters.\n");
 
 
  if (strcmp(argv[1],"-f")==0) /*Use parameters file with the name given as 2nd parameter*/
  {
    if (narg < 3) errquit("Parameters file name is not given.\n");
    fp=file_open(argv[2],"rt");
    ask=0;
  }
  else
  {
    if (strcmp(argv[1],"-s")==0) /*Read parameters from standart input*/
    {
      fp=stdin;
      ask=1;
    }
    else
    {
      printf("Unrecognized flag %s.\n",argv[1]);
      exit(1);
    }
  }
  
  if (ask) fprintf(stdout,"Input punch parameter file name (%d chars max) >  \n",FLNLN);
  fscanf(fp,"%s",src_file);
	xs=1.5;
	ys=0.5;
	zs=0.1;
  if (ask) fprintf(stdout,"Input receiver parameter file name (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",rec_file);

    fi=file_open(rec_file,"rt");  //np[l],xp[l],yp[l],zp[l],i,j,k,t

      fscanf(fi,"%d",&Np);
      np=(int*)fmemalloc(Np,sizeof(int));
      xp=(float*)fmemalloc(Np,sizeof(float));
      yp=(float*)fmemalloc(Np,sizeof(float));
      zp=(float*)fmemalloc(Np,sizeof(float));
      ii=(int*)fmemalloc(Np,sizeof(int));
      jj=(int*)fmemalloc(Np,sizeof(int));
      kk=(int*)fmemalloc(Np,sizeof(int));
      t=(float*)fmemalloc(Np,sizeof(float));

      for (i=0;i<Np;i++)
        fscanf(fi,"%d%g%g%g%d%d%d%g",np+i,xp+i,yp+i,zp+i,ii+i,jj+i,kk+i,t+i);
    fclose(fi);
  printf("Reciver parameters read.\n");

  if (ask) fprintf(stdout,"Input velocity on the upper surface V0 and gradient > \n");
  fscanf(fp,"%g%g",&V0, &gradV);
  if (ask) fprintf(stdout,"Input resulting file name (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",out_file);
  if (!ask) fclose(fp);
 
//V=V0+z*gradV
//C=c0(1+z*beta)
c0=V0+zs*gradV;
beta=gradV/c0;

  fp=fopen(out_file,"wt");
  for (l=0;l<Np;l++)
  {
      //td=(yp[l]-ys)/v1*sin(teta+ffi)+t1d;
      td=2/c0/beta*asinh(sqrt((yp[l]-ys)*(yp[l]-ys)+(xp[l]-xs)*(xp[l]-xs))*beta/2);

      fprintf(fp,"%5d %10g %10g %10g %5d %5d %5d %10g %10g %10g %10g \n", np[l],xp[l],yp[l],zp[l],ii[l],jj[l],kk[l],t[l],td,(td-t[l]),((td-t[l])/td));

  }
  fclose(fp);
}
