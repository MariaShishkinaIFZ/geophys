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
  float a1, a2, a3;
  float xs, ys, zs;
  float nx, ny, nz;
  float alfa, delta, ffi, teta;
  float v1, v2;
  float A,B,C,D;
  float t1d,zd,td,zd1,zd2;
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
  
  if (ask) fprintf(stdout,"Input source point coordinates: xs ys zs > \n");
  fscanf(fp,"%g%g%g",&xs,&ys,&zs);

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
  if (ask) fprintf(stdout,"Input velocity of upper layer > \n");
  fscanf(fp,"%g",&v1);
  if (ask) fprintf(stdout,"Input velocity of lower layer > \n");
  fscanf(fp,"%g",&v2);

  if (ask) fprintf(stdout,"Input boundary mark poin x y z, - coordinates > \n");
  fscanf(fp,"%g%g%g",&a1,&a2,&a3);
  if (ask) fprintf(stdout,"Input boundary hade and strech angles > \n");
  fscanf(fp,"%g%g",&alfa, &delta);
  if (ask) fprintf(stdout,"Input resulting file name (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",out_file);
  if (!ask) fclose(fp);
 
  if (delta>=90 || delta<0 )
  {
    fprintf(stderr,"Hade should be  in interval between 0 and 90 degrees!\n");
    exit(0);
  }


  alfa =(alfa*(M_PI/180.0));
  delta=(delta*(M_PI/180.0));

  nx=sin(delta)*sin(alfa); 	//nx(x - a1)+ny(y-a2)+nz(z-a3)=0
  ny=sin(delta)*cos(alfa);
  nz=cos(delta);

  A=nx;				//Ax+By+Cz+D=0
  B=ny;
  C=nz;
  D=(-A*a1-B*a2-C*a3);

  ffi=delta;

  zd=sqrt((A*xs+B*ys+C*zs+D)*(A*xs+B*ys+C*zs+D))/sqrt(A*A+B*B+C*C);

   fprintf(stdout,"A=nx= %g >  \n",A);
   fprintf(stdout,"B=ny= %g >  \n",B);
   fprintf(stdout,"C=nz= %g >  \n",C);
   fprintf(stdout,"D= (-nx*a1-ny*a2-nz*a3) = %g >  \n",D);

   fprintf(stdout,"a1= %g >  \n",a1);
   fprintf(stdout,"a2= %g >  \n",a2);
   fprintf(stdout,"a3= %g >  \n",a3);

   fprintf(stdout,"zd= %g >  \n",zd);

  zd1=cos(ffi)*(tan(ffi)*ys+0.4);

   fprintf(stdout,"zd1= %g >  \n",zd1);

  zd2=cos(ffi)*(a3-zs);

   fprintf(stdout,"zd2= %g >  \n",zd2);

  teta=asin(v1/v2);
 
  t1d=2.0*zd2/v1*cos(teta);

  fp=fopen(out_file,"wt");
  for (l=0;l<Np;l++)
  {

      if( (sqrt((yp[l]-ys)*(yp[l]-ys)+(xp[l]-xs)*(xp[l]-xs))/v1) >= (sqrt((yp[l]-ys)*(yp[l]-ys)+(xp[l]-xs)*(xp[l]-xs))/v1*sin(teta+ffi)+t1d) )
      {
      td=sqrt((yp[l]-ys)*(yp[l]-ys)+(xp[l]-xs)*(xp[l]-xs))/v1*sin(teta+ffi)+t1d;

      fprintf(fp,"%5d %10g %10g %10g %5d %5d %5d %10g %10g %10g %10g \n", np[l],xp[l],yp[l],zp[l],ii[l],jj[l],kk[l],t[l],td,(td-t[l]),((td-t[l])/td));
      }

      else           td=sqrt((yp[l]-ys)*(yp[l]-ys)+(xp[l]-xs)*(xp[l]-xs))/v1;

  }
  fclose(fp);
}


