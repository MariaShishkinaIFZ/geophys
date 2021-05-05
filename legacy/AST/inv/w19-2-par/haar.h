/*Haar wavelet functions definitions*/
/*Strasbourg, 05.2005*/

/*Returns 2^n, n must be integer less than 31 and non-negative*/
long l2n(long n);

/*Returns Haar wavelet function at the scale l, and shifts m and n, and type k*/
/*i and j are pixel indices in the image, supposed to change from 1 to N=2^J*/
/*No change is performed on the validity of the arguments, use at you own risk*/
/*Valid values are: 1<=l<=J, 0<=m,n<=2^(J-l)-1, 1<=i,j<=n*/
float haar_w(long l, long m, long n, long k, long i, long j);

/*Returns Haar scaling function at the scale l, and shifts m and n*/
/*i and j are pixel indices in the image, supposed to change from 1 to N=2^J*/
/*No change is performed on the validity of the arguments, use at you own risk*/
/*Valid values are: 1<=l<=J, 0<=m,n<=2^(J-l)-1, 1<=i,j<=n*/
float haar_f(long l, long m, long n, long i, long j);

/*Calculates the minimum and maximum values of pixel indices that bounds the Haar wavelet*/
/*support area at the scale l, and shifts m and n, and type k*/
/*These values are returned in imin, imax, jmin, jmax*/
/*No change is performed on the validity of the arguments, use at you own risk*/
/*Valid values are: 1<=l<=J, 0<=m,n<=2^(J-l)-1, 1<=i,j<=n*/
void haar_supp(long l, long m, long n, long *imin, long *imax, long *jmin, long *jmax);
