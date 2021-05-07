#include <cmath>
#include <cstdio>
#include <cstdlib>

/*Returns 2^n, n must be integer less than 31 and non-negative*/
long l2n(long n)
{
  long s=1;
  if (n==0) return 1;
  return s << n;
}

float haar_w(long l, long m, long n, long k, long i, long j)
{
  i = (i-1)/l2n(l-1)-2*m; j=(j-1)/l2n(l-1)-2*n;
  switch (k)
  {
    case (1) : if ((j==1) && ((i==0) || (i==1))) return 0.25/(float)l2n(l);
    	       if ((j==0) && ((i==0) || (i==1))) return -0.25/(float)l2n(l);
	       return 0.0;
    case (2) : if ((i==1) && ((j==0) || (j==1))) return 0.25/(float)l2n(l);
    	       if ((i==0) && ((j==0) || (j==1))) return -0.25/(float)l2n(l);
	       return 0.0;
    case (3) : if (((i==0) && (j==0)) || ((i==1) && (j==1))) return 0.25/(float)l2n(l);
    	       if (((i==0) && (j==1)) || ((i==1) && (j==0))) return -0.25/(float)l2n(l);
	       return 0.0;
    default  : fprintf(stderr,"Error: k=%ld\n",k); exit(EXIT_FAILURE);
  } 
}

float haar_f(long l, long m, long n, long i, long j)
{
  i = (i-1)/l2n(l-1)-2*m; j=(j-1)/l2n(l-1)-2*n;
  if (((i==0) || (i==1)) && ((j==0) || (j==1))) return 1.0/l2n(l);
  else return 0.0;
}

void haar_supp(long l, long m, long n, long *imin, long *imax, long *jmin, long *jmax)
{
  *imin = l2n(l)*m+1; *imax = l2n(l)*(m+1);
  *jmin = l2n(l)*n+1; *jmax = l2n(l)*(n+1);
}

