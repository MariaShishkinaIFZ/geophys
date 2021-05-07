#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*x0 - начало отсчёта, w - ширина клетки, t - полуширина переходной зоны*/
double ch1d(double x, double x0, double w, double t)
{
  long n = floor((x-x0)/w); /*номер клетки*/
  long a = pow(-1.0,n);
  
  if ((x-x0>=n*w+t) && (x-x0<=(n+1)*w-t)) return a;
  else
  {
    if (x-x0<n*w+t) return a*sin(M_PI*(x-x0-n*w)/(2.0*t));
    else return a*sin(M_PI*((n+1)*w-(x-x0))/(2.0*t));
  }
}

int main(int narg, char **argv)
{
  int N, i;
  double x0, dx, w, t;
  
  if (narg<6)
  {
    fprintf(stderr,"Use: chessplot x0 dx N w t\n");
    exit(0);
  }
  
  sscanf(argv[1],"%lg",&x0);
  sscanf(argv[2],"%lg",&dx);
  sscanf(argv[3],"%d",&N);
  sscanf(argv[4],"%lg",&w);
  sscanf(argv[5],"%lg",&t);
  
  for (i = 0;i<=N;i++)
    fprintf(stdout,"%lg %lg %lg %lg\n",dx*i,dx*i-x0,floor((dx*i-x0)/w),ch1d(dx*i,x0,w,t));
}
