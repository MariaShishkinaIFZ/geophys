/*---Interface+sub-interface+3d upper layer velocity and density inversion----*/
/*---Version W.10 last edited 20.06.08. by S.Tikhotsky, Strasbourg------------*/
/*----------------------------------------------------------------------------*/
/*---NEW from 3D-W1: Smooth. constraints in vertical direction for upper layer*/
/*---NEW from 3D-W2: Back-interpolation on the find grid and construction of--*/
/*-------------------the next model added.------------------------------------*/
/*---NEW from 3D-W4: Takes fine grid velocity and interface model on the------*/
/*-------------------input, averages it to get reference velocity and---------*/
/*-------------------interface models. Can either apply smoothing constrain to*/
/*-------------------the solution (i.e. slowness variance at the current------*/
/*-------------------iteration) or to the target slowness model (i.e. Sref+dS)*/
/*---NEW from 3D-W5: Bugs at the stage of the wavelet paramerisation----------*/
/*-------------------construction corrected: in W5 version number of rays per-*/
/*-------------------wavelet support was actually calculated as number of-----*/
/*-------------------inversion grid primary cells crossed by rays. Also shifts*/
/*-------------------of wavelet pattern wrt layer were incorrect--------------*/
/*---NEW from 3D-W6: Structural similarity (colinearity of the gradients) of--*/
/*-------------------density and slowness fields added instead of linear corr.*/
/*-------------------Linearized expression of vector product used (3 eq.) for-*/
/*-------------------each model point-----------------------------------------*/
/*---------07.04.06--Smoothing shifted from velocity do density coeffs.-------*/
/*---06.2008, S.Tikhotsky, EOST, Strasbourg-----------------------------------*/
/*---NEW IN 3D-W10:  Computing of the resolution matrix-----------------------*/
/*---------20.06.08--Wavelets numbering system was incorrect. Fixed and tested*/
/*-------------------Scaling function is now always included as the 1st member*/
/*-------------------of the series (may be not wise, but at least work corr.).*/
/*-------------------Other wavelet functions with l=J are NOT included, search*/
/*-------------------starts from l=J-1.---------------------------------------*/
/*---------23.06.08--Two critical bugs fixed. See BUG230608 label in the text.*/
/*-----------3D-W12--Hit counts normalized by 2^n-----------------------------*/
/*-------------------"Treppe" choosing rule based on normalized hit conts &---*/
/*-------------------angular coverage-----------------------------------------*/
/*---------05.07.09. Petelino-------------------------------------------------*/
/*-------------------W14: New interpolation procedure without pre-averaging---*/
/*-------------------vortices and pre-averaging over neighbouring cells-------*/
/*-------------------Some bugs were fixed meanwhile (from 3D-W12)-------------*/
/*-------------------Different weights for smoothing in X and Y directions are*/
/*-------------------introduced for better perofrmance in quasi-2D models-----*/
/*--3D-W14: New approach to prepare UL_dscm. Only those cells that were really*/
/*----inverted are used for wavelet summation - therefore no "tails", then any*/
/*-------cell that have reconstructed neighbours is set to mean value of those*/
/*-------------------neighbours - for future interplation onto fine grid------*/
/*--3D-W18: More intence resolution martix analysis, incl. non-diagonal terms-*/
/*--according to suggestions of T.B.Yanovskaya--------------------------------*/
/*--3D-W19 skipped------------------------------------------------------------*/
/*--3D-W20: Rearranging wavelet basis based on the actual resolution estimated*/
/*--from the resolution matrix. First inversion and resolution matrix are made*/
/*--based on the ray and angular coverage as in W18---------------------------*/
/*--p2_Nr; p2_Np; p2_Theta are now float (as in W19) allouing the fractional--*/
/*--powers of two normalization of empirical resolution measures--------------*/
/*----------------------------------------------------------------------------*/
/*---BASIC FEATURES FROM GRAS-W4:---------------------------------------------*/
/*---Haar wavelets decompositon used for the reparametrization----------------*/
/*---Different adaptive wavelet basises for slowness and interface------------*/
/*---NEW from W.2: same basis for slowness and density, based on ray coverage-*/
/*---Reflection and refraction traveltimes used-------------------------------*/
/*---LSQR solution------------------------------------------------------------*/
/*---Correlation between ALL slowness and density wavelet coefficients--------*/
/*---A priory information about the interface position at some points---------*/
/*---NEW from W.1: Smoothing constraints on the slowness and interface--------*/
/*---OPTIMIZED from W.3: Equations matrix computation optimized taking the----*/
/*-----------------------advantage of limited Haar wavelets support-----------*/
/*----------------------------------------------------------------------------*/

#define NOMINMAX
#include <iostream> 
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>
#include <limits>
#include <map>

#include "datatypes.h"

#include "cilk_util.h"
#include "utility.h"

#include "triple.h"
#include "csc.h"
#include "bicsb.h"
#include "spvec.h"

#include <cstdio>
#include <cstring>
#include <values.h>

#include "subio.h"
#include "submem.h"
#include "nrutil.h"
#include "lsqr.h"
#include "haar.h"

#include "gmtio.h"

using namespace std;

#define GEOMFILE "model-geometry.gras"
#define STATISTICS_FILE "gras.stat"

#define RESOLUTION_MATRIX_FILE "R_matrix.txt"
#define RESOLUTION_MATRIX_GMT_FILE "R_matrix.gmt"
#define RESOLUTION_MARKUP_FILE "R_matrix.xy"
#define WAVELETS_STRUCTURE_FILE "wavelets_used.txt"
#define COV_Z_FILE "cov_z.txt"
#define FISHER_SI_FILE "fish_si.txt"
#define FISHER_UPL_FILE "fish_upl.txt"
#define RESOLUTION_SI_FILE "resolution_SI.txt"
#define RESOLUTION_Z_FILE "resolution_Z.txt"
#define RESOLUTION_UPL_FILE "resolution_UPL.txt"

#define FLNLN 256 	/*Maximum length of the file name*/
#define STRLN 256
#define TSTRLN 256

#define G 6.672

#define DUMMY_LIM 1e9
#define DUMMY_VALUE 1e10

#define ACC_CONST 1e-5

#define ASSUMED_ACCURACY 0

/*Wavelet scale and shifts information*/
/*W11: k=0 stays for the scaling function, l=J,m,n=0*/
typedef struct
	{
	  short l, m, n, k;
	}
	Ww_type;

typedef struct
	{
	  long ii_min, ii_max, jj_min, jj_max; /*Position in the original gen. grid*/
	  float h_top, h_bot; /*Layer minimum top and bottom depth*/
	  long NC;	/*Actual number of cells in the layer to be inverted for*/
	  long *CV;	/*Layer cells indices*/
	  long nx, ny; /*Layer dimensions*/
	  long WI_J;   /*Wavelet image side length power of 2*/
	  long WI_N;   /*Wavelet image side length*/
	  long w_N;    /*Number of wavelet SUPPORTS in the expansion*/
	  long w_ncf;  /*Number of wavelet COEFFICIENTS = 3*w_N+1*/
	  long W_Si, W_Sj; /*Shifts of the model area over the wavelets image in X and Y direction*/
	  Ww_type *WV; /*Wavelets used in the expansion*/
	}
	Grid_Layer_type;

/*Defines the structure of block sparse matrix, used to optimize matrix-vector multiplication*/
typedef struct
	{
	  long NBR;		/*Number of block rows*/
	  long NBC;		/*Number of block columns*/
	  long *Vbr;		/*First indices of each block row (NBR+1 vector)*/
	  long *Vbc;		/*First indices of each block column (MBR+1 vector)*/
	  short **ABS;		/*Matrix blocks indicators (0 - zero block, 1 - non-zero block)*/
	  float **A;		/*Matrix data*/
	}
	block_matrix_structure;

typedef struct
	{	  
	 float l,x,y,z;
	}
	ray_components;	
 
typedef std::map<long, ray_components> rays_map;

typedef std::map<long, float> reflct_points_map;

double delta_kron(long i, long j)
{
  return (i==j) ? 1.0 : 0.0;
}

int eq_digital(float x1, float x2, float eps)
{
  if (fabs(x1-x2)<eps) return 1;
  return 0;
}

long roun(double x)
{
  long i;
  i = floor(x);
  if (x-i >= 0.5) i++;
  return i;
}

int imax(int i1, int i2)
{
  return i1>i2 ? i1 : i2;
}

/*Return the index of cell from its position in [1..Nrow][1..Ncol] matrix*/
/*Do NOT check for the possible runs outside the matrix dimentions*/
long cln(long j, long i, long Ncol)
{
  return (j-1)*Ncol+i;
}

/*Looks for the position (i.e. i & j indices) in [1..Nrow][1..Ncol] matrix*/
/*from the global cell index*/
void pos_from_cln(long *j, long *i, long c, long Ncol)
{
  *j = (c-1)/Ncol+1;
  *i = c - (*j-1)*Ncol;
}

void alloc_bms(block_matrix_structure *V)
{
  long i;
  V->Vbr=(long*)fmemalloc(V->NBR+1,sizeof(long))-1;
  V->Vbc=(long*)fmemalloc(V->NBC+1,sizeof(long))-1;
  V->ABS=(short**)fmemalloc(V->NBR,sizeof(short*))-1;
  for (i=1;i<=V->NBR;i++) V->ABS[i]=(short*)fmemalloc(V->NBC,sizeof(short))-1;
}

void free_bms(block_matrix_structure *V)
{
  long i;
  for (i=1;i<=V->NBR;i++) fmemfree((char*)(V->ABS[i]+1));
  fmemfree((char*)(V->ABS+1));
  fmemfree((char*)(V->Vbc+1));
  fmemfree((char*)(V->Vbr+1));
}

void aprod_bicsb( long  mode,
	          dvec  *x, 
                  dvec  *y,
                  void  *prod )
/*
*     ------------------------------------------------------------------
*     This is the matrix-vector product routine required by LSQR
*     for a sparse matrix in BiCSB format.
*     BiCSB matrix class is passed through prod
*     ------------------------------------------------------------------
*/
{
  BiCsb<VALUETYPE, INDEXTYPE> *data;
  
  data = (BiCsb<VALUETYPE, INDEXTYPE> *) prod;  
/*  
*     Compute  Y = Y + A*X
*/
  if( mode == 0 )
     bicsb_gaxpy(*data, x->elements, y->elements);
/*  
*     Compute  X = X + A^T*Y
*/  
  if( mode == 1 )
     bicsb_gaxpy_trans(*data,y->elements,x->elements);
  return;
}

int ww_compare(const Ww_type *w1, const Ww_type *w2)
{
  if ((w1->l==w2->l) && (w1->m==w2->m) && (w1->n==w2->n))
    return 1;
  else
    return 0;
}

void ww_set(const Ww_type *src, Ww_type *dest)
{
  dest->l=src->l; dest->m=src->m; dest->n=src->n;
}

/*Finds the common elements in V1 and V2 and place them in the corresponding cells at the head*/
/*of both vectors. Returns the number of common elements found.*/
/*Both vectors are supposed to be initialy sorted in the decreasing order of scale l*/
int ww_common(Ww_type *V1, int n1, Ww_type *V2, int n2)
{
  int i,j,k,found,Ncommon,nr;
  Ww_type wx;

  Ncommon=0;
  if (n1<=n2)
  {
     nr=n1;
     for (i=1;i<=nr;i++)
     {
       if (ww_compare(V1+i,V2+i)==0)
       {
         found=0;
	 for (j=i+1;j<n2;j++)
	   if (ww_compare(V1+i,V2+j)==1)
	   {
	     for (k=j-1;k>=i;k--) ww_set(V2+k,V2+k+1);
	     ww_set(V1+i,V2+i);
	     found=1; Ncommon++; break;
	   }
	 if (found==0)
	 {
	   ww_set(V1+i,&wx);
	   for (k=i+1;k<=n1;k++) ww_set(V1+k,V1+k-1);
	   ww_set(&wx,V1+n1);
	   i--; nr--;
	 }
       }
       else
         Ncommon++;
     }
  }
  else
  {
     nr=n2;
     for (i=1;i<=nr;i++)
     {
       if (ww_compare(V2+i,V1+i)==0)
       {
         found=0;
	 for (j=i+1;j<n1;j++)
	   if (ww_compare(V2+i,V1+j)==1)
	   {
	     for (k=j-1;k>=i;k--) ww_set(V1+k,V1+k+1);
	     ww_set(V2+i,V1+i);
	     found=1; Ncommon++; break;
	   }
	 if (found==0)
	 {
	   ww_set(V2+i,&wx);
	   for (k=i+1;k<=n2;k++) ww_set(V2+k,V2+k-1);
	   ww_set(&wx,V2+n2);
	   i--; nr--;
	 }
       }
       else
         Ncommon++;
     }
  }
  return Ncommon;
}

float ev_haar(long wq, Ww_type *V, int i, int j)
{
  if (V[wq].k==0) return haar_f(V[wq].l,V[wq].m,V[wq].n,i,j);
  else	return haar_w(V[wq].l,V[wq].m,V[wq].n,V[wq].k,i,j);
}

void find_haar_indices(long wq, Ww_type *V, long *l, long *m, long *n, long *k)
{
  *l=V[wq].l;
  *m=V[wq].m;
  *n=V[wq].n;
  *k=V[wq].k;
}

void find_haar_supp(long wq, Ww_type *V, long *imin, long *imax, long *jmin, long *jmax)
{
  haar_supp(V[wq].l,V[wq].m,V[wq].n,imin,imax,jmin,jmax);
}

void minmax(float const *buf, long N1, long N2, float *min, float *max)
{
 long i;

 *min=*max=buf[N1];
 for (i=N1+1;i<=N2;i++)
 {
   if (buf[i]<*min) *min=buf[i];
   if (buf[i]>*max) *max=buf[i];
 }
}

long lmax(long n1, long n2)
{
  return n1>n2 ? n1 : n2;
}

long lmin(long n1, long n2)
{
  return n1<n2 ? n1 : n2;
}

long find_long(long *buf, long n, long val)
{
  long i, il, ir, i1;

  il=1; ir=n; i=n/2;
  if (buf[il]==val) return il;
  if (buf[ir]==val) return ir;
  do
  {
    if (buf[i]==val) return i;
    if (buf[i]>val)
    {
      if ((i==il) || (i-1==il)) return -1;
      i1=i;
      i=(i+lmax(1,il))/2;
      ir=i1;
    }
    else
    {
      if ((i==ir) || (i+1==ir)) return -1;
      i1=i;
      i=(i+lmin(n,ir))/2;
      il=i1;
    }
  }
  while (1);
}

/*Return index of cell with the throught_index in cells array of layer nl*/
/*return -1 if cell with the throught_index is not presented in the layer nl*/
long find_cell_in_layer(Grid_Layer_type *L, long NNl, long nl, long through_index)
{
  if ((nl<1) || (nl>NNl)) return -1;
  return find_long(L[nl].CV,L[nl].NC,through_index);
}

/*Return 1 if cell with the index ind in the cells array of the layer nl is internal*/
/*0 otherwise*/
/*ATTENTION: ind is the array index and NOT the through index of cell in (Nrow)x(Ncol) matrix*/
/*Ncol is the number of columns in the general inversion matrix*/
long cell_is_internal(Grid_Layer_type *L, long NNl, long Ncol, long nl, long ind)
{
  long i,j;
  pos_from_cln(&j,&i,L[nl].CV[ind],Ncol);
  return
  (find_cell_in_layer(L,NNl,nl,cln(j-1,i,Ncol))>=0) &&
  (find_cell_in_layer(L,NNl,nl,cln(j+1,i,Ncol))>=0) &&
  (find_cell_in_layer(L,NNl,nl,cln(j,i-1,Ncol))>=0) &&
  (find_cell_in_layer(L,NNl,nl,cln(j,i+1,Ncol))>=0) &&
  (find_cell_in_layer(L,NNl,nl-1,cln(j,i,Ncol))>=0) &&
  (find_cell_in_layer(L,NNl,nl+1,cln(j,i,Ncol))>=0);
}

void read_3D_volume(float ***Buf, long Nz, long Ny, long Nx, const char *filename)
{
  long j,k;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fread(Buf[k][j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_3D_volume(float ***Buf, long Nz, long Ny, long Nx, const char *dir, const char *filename)
{
  long j,k;
  FILE *fo;

  fo=dir_file_open(dir,filename,"wb");
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      fwrite(Buf[k][j],sizeof(float),Nx,fo);
  fclose(fo);
}

void read_2D_surface(float **Buf, long Ny, long Nx, const char *filename)
{
  long j;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (j=0;j<Ny;j++)
    fread(Buf[j],sizeof(float),Nx,fi);
  fclose(fi);
}

void write_2D_surface(float **Buf, long Ny, long Nx, const char *dir, const char *filename)
{
  long j;
  FILE *fo;

  fo=dir_file_open(dir,filename,"wb");
  for (j=0;j<Ny;j++)
    fwrite(Buf[j],sizeof(float),Nx,fo);
  fclose(fo);
}

int av2d(float **src, float **dest, long nrow, long ncol, long kav)
{
  long i, j, ii, jj;
  long M, N;
  float w;

  N=(ncol-1)/kav;
  M=(nrow-1)/kav;

  if (((N*kav)!=(ncol-1)) || ((M*kav)!=(nrow-1)))
    return -1; /*Error*/

  for (jj=1;jj<=M;jj++)
    for (ii=1;ii<=N;ii++)
    {
      dest[jj][ii]=0.0; w=0.0;
      /*Internal cell points*/
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
        for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
	  if (src[j][i]<DUMMY_LIM)
	  {
            dest[jj][ii]+=src[j][i];
	    w+=1.0;
	  }
      /*Lower cell boundary*/
      j=(jj-1)*kav;
      i=(ii-1)*kav;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      i=ii*kav-1;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      /*Upper cell boundary*/
      j=jj*kav-1;
      i=(ii-1)*kav;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      i=ii*kav-1;
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[j][i];
	w+=0.25;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      /*Left cell boundary*/
      i=(ii-1)*kav;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      /*Right cell boundary*/
      i=ii*kav-1;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if (src[j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.5*src[j][i];
	w+=0.5;}
      if (w>0)
        dest[jj][ii]/=w;
      else
        dest[jj][ii]=DUMMY_VALUE;
    }
    return 1;
}

int av3d_ins_lay(float ***src,float **top,float **bot, float ***dest,
                 long nlay, long nrow, long ncol, long kav,
		 long kmin, long kz, long NL,
		 float z0, float h)
{
  long i, j, k, k1, k2, ii, jj, kk;
  long M, N;
  float w, z;

  N=(ncol-1)/kav;
  M=(nrow-1)/kav;

  if (((N*kav)!=(ncol-1)) || ((M*kav)!=(nrow-1)))
    return -1; /*Error*/

  for (kk=1;kk<=NL;kk++)
  {
   k1=kmin+(kk-1)*kz;
   k2=kmin+kk*kz-1;
   
   //fprintf(stderr,"KK=%d NL=%d k1=%d k2=%d\n",kk,NL,k1,k2);
   
   for (jj=1;jj<=M;jj++)
    for (ii=1;ii<=N;ii++)
    {
      dest[kk][jj][ii]=0.0; w=0.0;
      /*Internal layers*/
      for (k=k1+1;k<=k2-1;k++)
      {
       z=z0+k*h;
       /*Internal cell points*/
       for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
        for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
	  if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h))
	  {
            dest[kk][jj][ii]+=src[k][j][i];
	    w+=1.0;
	  }
       /*Front cell boundary*/
       j=(jj-1)*kav;
       i=(ii-1)*kav;
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       i=ii*kav-1;
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
       /*Rear cell boundary*/
       j=jj*kav-1;
       i=(ii-1)*kav;
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       i=ii*kav-1;
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.25*src[k][j][i];
	 w+=0.25;}
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
       /*Left cell boundary*/
       i=(ii-1)*kav;
       for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
       /*Right cell boundary*/
       i=ii*kav-1;
       for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
         dest[kk][jj][ii]+=0.5*src[k][j][i];
	 w+=0.5;}
      }

      //fprintf(stderr,"Internal nodes OK\n",kk,NL,k1,k2);

      k=k1; /*Cell top*/
      z=z0+k*h;
      /*Internal cell points*/
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
	 if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h))
	 {
           dest[kk][jj][ii]+=0.5*src[k][j][i];
	   w+=0.5;
	 }
      /*Front cell boundary*/
      j=(jj-1)*kav;
      i=(ii-1)*kav;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav-1;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Rear cell boundary*/
      j=jj*kav-1;
      i=(ii-1)*kav;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav-1;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Left cell boundary*/
      i=(ii-1)*kav;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Right cell boundary*/
      i=ii*kav-1;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}

      //fprintf(stderr,"Top OK\n");

      k=k2; /*Cell bottom*/
      z=z0+k*h;
      /*Internal cell points*/
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
       for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
	 if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h))
	 {
           dest[kk][jj][ii]+=0.5*src[k][j][i];
	   w+=0.5;
	 }
      /*Front cell boundary*/
      j=(jj-1)*kav;
      i=(ii-1)*kav;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav-1;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Rear cell boundary*/
      j=jj*kav-1;
      i=(ii-1)*kav;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      i=ii*kav;
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.125*src[k][j][i];
	w+=0.125;}
      for (i=(ii-1)*kav+1;i<=ii*kav-2;i++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Left cell boundary*/
      i=(ii-1)*kav;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}
      /*Right cell boundary*/
      i=ii*kav-1;
      for (j=(jj-1)*kav+1;j<=jj*kav-2;j++)
      if ((src[k][j][i]<DUMMY_LIM) && (z >= top[j][i]) && (z < bot[j][i]-ACC_CONST*h)) {
        dest[kk][jj][ii]+=0.25*src[k][j][i];
	w+=0.25;}

      //fprintf(stderr,"Bottom OK\n");

      if (w>0)
        dest[kk][jj][ii]/=w;
      else
        dest[kk][jj][ii]=DUMMY_VALUE;
    }
   } /*for kk*/
   return 1;
}

int av3d_sub_intf(float ***src,float **intf, float **dest,
                  long nlay, long nrow, long ncol, long kav,
		  float z0, float h)
{
  long i, j, k, ii, jj;
  long M, N;
  float w;

  N=(ncol-1)/kav;
  M=(nrow-1)/kav;

  if (((N*kav)!=(ncol-1)) || ((M*kav)!=(nrow-1)))
    return -1; /*Error*/

  for (jj=1;jj<=M;jj++)
    for (ii=1;ii<=N;ii++)
    {
      dest[jj][ii]=0.0; w=0.0;
      /*Internal cell points*/
      for (j=(jj-1)*kav+2;j<=jj*kav-1;j++)
        for (i=(ii-1)*kav+2;i<=ii*kav-1;i++)
	{
	  k = (intf[j][i]-z0)/h;
	  if (k+1>=nlay) return 0;
	  if (src[k+1][j][i]<DUMMY_LIM)
	  {
            dest[jj][ii]+=src[k+1][j][i];
	    w+=1.0;
	  }
	}
      /*Lower cell boundary*/
      j=(jj-1)*kav+1;
      i=(ii-1)*kav+1;
      k = (intf[j][i]-z0)/h;
      if (k+1>=nlay) return 0;
      if (src[k+1][j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[k+1][j][i];
	w+=0.25;}
      i=ii*kav;
      k = (intf[j][i]-z0)/h;
      if (k+1>=nlay) return 0;
      if (src[k+1][j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[k+1][j][i];
	w+=0.25;}
      for (i=(ii-1)*kav+2;i<=ii*kav-1;i++)
      {
        k = (intf[j][i]-z0)/h;
        if (k+1>=nlay) return 0;
	if (src[k+1][j][i]<DUMMY_LIM) {
          dest[jj][ii]+=0.5*src[k+1][j][i];
	  w+=0.5;}
      }
      /*Upper cell boundary*/
      j=jj*kav;
      i=(ii-1)*kav+1;
      k = (intf[j][i]-z0)/h;
      if (k+1>=nlay) return 0;
      if (src[k+1][j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[k+1][j][i];
	w+=0.25;}
      i=ii*kav;
      k = (intf[j][i]-z0)/h;
      if (k+1>=nlay) return 0;
      if (src[k+1][j][i]<DUMMY_LIM) {
        dest[jj][ii]+=0.25*src[k+1][j][i];
	w+=0.25;}
      for (i=(ii-1)*kav+2;i<=ii*kav-1;i++)
      {
        k = (intf[j][i]-z0)/h;
        if (k+1>=nlay) return 0;
	if (src[k+1][j][i]<DUMMY_LIM) {
          dest[jj][ii]+=0.5*src[k+1][j][i];
	  w+=0.5;}
      }
      /*Left cell boundary*/
      i=(ii-1)*kav+1;
      for (j=(jj-1)*kav+2;j<=jj*kav-1;j++)
      {
        k = (intf[j][i]-z0)/h;
        if (k+1>=nlay) return 0;
	if (src[k+1][j][i]<DUMMY_LIM) {
          dest[jj][ii]+=0.5*src[k+1][j][i];
	  w+=0.5;}
      }
      /*Right cell boundary*/
      i=ii*kav;
      for (j=(jj-1)*kav+2;j<=jj*kav-1;j++)
      {
        k = (intf[j][i]-z0)/h;
        if (k+1>=nlay) return 0;
	if (src[k+1][j][i]<DUMMY_LIM) {
          dest[jj][ii]+=0.5*src[k+1][j][i];
	  w+=0.5;}
      }
      if (w>0)
        dest[jj][ii]/=w;
      else
        dest[jj][ii]=DUMMY_VALUE;
    }
    return 1;
}

int smooth_interface(float **V, int Nrow, int Ncol, float x0, float y0, float h, float R, float p)
{
  int i,j,k,l,dd;
  float x,y,x1,y1,s,w,c;
  float **Buf;

  alloc_float_matrix(&Buf,Nrow,Ncol);

  dd=(int)floor(R/h);

  for (j=0;j<Nrow;j++)
  {
    y=y0+h*j;
    for (i=0;i<Ncol;i++)
    {
      x=x0+h*i;
      s=0.0; w=0.0;
      for (l=-dd;l<=dd;l++)
      {
        if ((j+l >= 0) && (j+l < Nrow))
	{
	  y1=y0+h*(j+l);
	  for (k=-dd;k<=dd;k++)
	  {
	    if ((i+k >= 0) && (i+k < Ncol))
	    {
	      x1=x0+h*(i+k);
	      if (V[j+l][i+k]<DUMMY_LIM)
	      {
	        c=exp(-pow((double)sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y))/R,(double)p));
	        s+=V[j+l][i+k]*c;
	        w+=c;
	      }
	    }
 	  }
	}
      }
      if (w!=0) Buf[j][i]=s/w;
      else Buf[j][i]=DUMMY_VALUE;
    }
  }

  for (j=0;j<Nrow;j++)
    for (i=0;i<Ncol;i++)
      V[j][i]=Buf[j][i];

  free_float_matrix(Buf,Nrow,Ncol);

  return 1;
}

/*Returns 1 if n and T lay above the "treppe" determined by n1,n2,T1,T2*/
int treppe_crit(float n1, float n2, float T1, float T2, float n, float T)
{
  int reject = 0;

  if (n<=n1) reject=1;
  else {
    if (n<=n2) {if (T<T2) reject=1;}
    else {if (T<T1) reject=1;}
  }
  if (reject==1) return 0; else return 1;
}

int resolution_test(char value_type, float threshold, float value)
{
  if (value>DUMMY_LIM) return 0;
  if (value_type == 's')
  {
    if (value < threshold) return 1; else return 0;
  }
  else
  {
    if (value > threshold) return 1; else return 0;
  }
}

int cilk_main(int narg, char **argv)
{
  long matrix_count;
  long Nx, Ny, Nz;	/*Seismic velocities grid dimentions (equal to number of nodes)*/
  long NNx, NNy, NNz;	/*Generalized grid dimentions (equal to number of cells)*/
  long k_gen;      	/*Generalization coefficient for inversion procedure*/
  long k_gen_v;		/*Generalization coefficient in depth (Z direction) for inversion procedure*/
  long N;		/*Number of EACH GENERALIZED GRID LAYER cells*/
  long N_int;		/*Number of EACH GENERALIZED GRID LAYER internal cells*/
  long NNall;		/*Number of UPPER PHYSICAL LAYER generalized grid cells (all grid layers)*/
  long Nl;		/*Number of reflected rays*/
  long Nr;		/*Number of head wave rays*/
  long Ncrp;        /*Number of basement crosspoints PAIRS*/
  long i,j,k,m;     /*Indices (j -> y, i -> x, k -> z, m spans data)*/
  long ii, jj, kk;  /*Indices for generalized model (jj -> y, ii -> x, kk -> z)*/
  long bigcelln;    /*Cell number in generalized model*/
  long r, c;        /*True indices in the LS equations matrix*/
  long nray;        /*Ray number (must be the same for same rays in ray path's and crosspoints)*/
  long Nw;          /*Number of a priory interface position points*/
  long nrec;        /*Receiver number (not used)*/
  long Num_E, Num_U;	/*Number of equations and number of unknowns*/
  long N0;		/*Starting row number - suppl. index for matrix composition*/
  long c0;		/*Starting column number - suppl. index for matrix composition*/
  long Num_SMZ_E;	/*Number of equations used for smoothing in Z direction*/
  long Num_SMZ_cells;   /*Number of cells used for smoothing in Z direction*/
  /*---------------------Wavelets staff-------------------------------*/
  long W_Ni, W_J;       /*Wavelets image dimention W_Ni=2^W_J*/
  long W_Si, W_Sj;      /*Shifts of the model area over the wavelets image in X and Y direction*/
  long wi, wj;          /*Indices over wavelets image*/
  long wl, wm, wn, wk;  /*Wavelets scale, shifts in X and Y direction and type*/
  long W_nl_1;          /*Number of wavelet SUPPORT areas used in the inversion for the sub-int. velocity*/
  long w_s1_N;		/*Number of wavelet coefficients used in the inversion for the sub-int. velocity*/
  int wq;		/*Wavelet sequence index*/
  int W_Nr;		/*Number of rays per wavelet support*/
  float p2_Nr;		/*Number of rays per wavelet support is nomalised by 2^(p2_Nr*l), l - wavelet order*/
  float W_Nr_Min_1;	/*Minimum number 1 (first step) of rays per wavelet support*/
  float W_Nr_Min_2;	/*Minimum number 2 (second step) of rays per wavelet support*/
  long wimin, wimax, wjmin, wjmax; /*Wavelet support bounding values*/
  Ww_type *WV_1;        /*Wavelets indices used in the sub-int. velocity inversion*/
  Ww_type *WV_Z;	/*Wavelets indices used in the interface depth inversion*/
  long wz_N;		/*Number of wavelet SUPPORT areas used in the interface depth inversion*/
  long w_z1_N;		/*Number of wavelet coefficients  in the interface depth inversion*/
  long W_Np_h, W_Np_r;	/*Number of refaraction and reflection points per wavelet support*/
  float p2_Np;		/*Number of points per wavelet support is nomalised by 2^(p2_Np*l), l - wavelet order*/
  float W_Np_Min;	/*Minimum number of refraction+reflection points per wavelet support*/
  int scalefunc_include;/*Include (1) or not (0) the scaling function in all wavelet expansions*/
  int **wv_map;		/*Wavelet support "map" to mark cells affected by the data inside the support*/
  int wv_cells_filled;	/*Number of filled cells in the wavelet support map*/
  float wv_map_sum;	/*Sum of all elements of the wavelet support map*/
  float wi_mean, wj_mean;/*Mean-weigted coordinates of the wavelet support map*/
  float wv_scatter;	/*Estimate of the scattering of the non-empty cells in the wavelet support map*/
  /*-----------------------Here comes the specific 3D model staff-------------------*/
  Grid_Layer_type *GL;  /*Grid layers information*/
  long uk_min, uk_max;  /*Minimum and maximum fine grid k (Z) indices, bounding the upper physical layer*/
  float UL_H_min;	/*Minimum upper physical layer model depth*/
  /*--------------------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  int raypos;           /*Indicator of ray's position with respect to basement*/
  int ncell;            /*Cell number in the ray path*/
  int use_down_hr, use_up_hr; /*Flags to use (1) or not (0) downgoing and upgoing headwave rays*/
  int use_down_rr, use_up_rr; /*Flags to use (1) or not (0) downgoing and upgoing reflected rays*/
 // int ray_included;     /*Flag used in the calculation of rays per wavelet support procedure*/
  int* ray_included;
  float x0, y0, z0;	/*Origin of the grid*/
  float h;              /*Seismic velocities grid discretisation step*/
  float dd;             /*Generalized grid plane (XY) discretisation step*/
  float ddz;		/*Generalized grid depth (Z) discretisation step*/
 // float **L;		/*Rays length's in the cells below the interface*/
 // float ***UHL;		/*Head wave rays ength's in the cells above the interface*/
 // float ***URL;		/*Reflected wave rays length's in the cells above the interface*/
  //float ***RVX;         /*Ray vector mean X component in cell for UPL*/
  //float ***RVY;		/*Ray vector mean Y component in cell for UPL*/
  //float ***RVZ;		/*Ray vector mean Z component in cell for UPL*/
  float ll;             /*Ray path length*/
  //float **WL;           /*Values of dt/dz derivatives for reflected waves*/
  //float **W;            /*Values of dt/dz derivatives for head waves basement-crossing points*/
  float *X;             /*Unknowns*/
  float *Err;           /*Estimates of the unknowns errors*/
  float *B;        /*Combined influence matrix and RHS*/
  float *dtl, *stl;	/*Delay times vector for reflected waves and corresponding standart deviations*/
  float *dt,*st;        /*Delay times vector for head waves and corresponding standart deviations*/
  float *T;		/*Upper layer top position (may be topography)*/
  float *Z;		/*Reference interface position (obtained by the averaging of fine Z model)*/
  float Tmin, Tmax;	/*Minimum and maximum of the upper layer top position*/
  float Zmin, Zmax;	/*Minimum and maximum of the interface position*/
  float *wx, *wy, *wz,
         *ws;           /*A priory information on the interface position*/
  float ww;             /*General weighting factor for interface information*/
  float lambda_x, lambda_y;         /*Roughness regularization parameters for x and y directions*/
  float lambda_z;	/*Roughness regularization parameter for vertical direction*/
  float beta;           /*Global regularization parameter*/
  float x,y,z,xp,yp,zp;
  float *mod_dt;
  //float **Xr, **Yr;     /*Ray vector mean X and Y components in cell*/
  float *wXr, *wYr;     /*Ray vector mean X and Y components in wavelet support*/
  float *wRVX;          /*Ray vector mean X component in wavelet support for UPL*/
  float *wRVY;          /*Ray vector mean Y component in wavelet support for UPL*/
  float *wRVZ;          /*Ray vector mean Z component in wavelet support for UPL*/
  float wXr_mean, wYr_mean, wZr_mean; /*Used for computin Fisher statistics*/
  float Theta_1;        /*Minimum value 1 (first step) of angular deviation for ray vectors (2D)*/
  float Theta_2;        /*Minimum value 2 (second step) of angular deviation for ray vectors (2D)*/
  float Theta_UPL_1; 	/*Minimum value 1 (first step) of angular deviation for ray vectors for UPL (3D)*/
  float Theta_UPL_2; 	/*Minimum value 2 (second step) of angular deviation for ray vectors for UPL (3D)*/
  float p2_Theta;		/*11th LINE: angular deviation is nomalised by 2^(p2_Theta*l), l - wavelet order*/
  float ss, ess;        /*Summation data and error buffers*/
  float *ss_l, *ss_r;
  float rL2, rChi2;	/*Residual vector L2 and Chi^2 norms*/
  time_t tm0, tmE;	/*Time veriable*/
  char time_string[TSTRLN]; /*String for time information*/
  /*---------Back-interpolation stuff------------------------------------------*/
  long jjj,iii,kkk,kkkk;	/*Additional counters*/
  float z_node;		/*Depth of the node*/
  float ***VFINE;	/*Fine gridded velocity model (nodes from 0,0,0)*/
  float ***dSFINE;  /*dS fine gridded velocity model (nodes from 0,0,0)*/
  float **SI_dsv_fine; /*Interpolated sub-interface slowness variation (nodes from 0,0)*/
  float ***UL_dscm;	/*Upper physical layer 3D coarse grid cells slowness model (cells from 1,1,1)*/
  float **UL_ft;	/*Upper physical layer top fine grid position (nodes from 0,0)*/
  float **Zintf;	/*Interface fine grid position (nodes from 0,0)*/
  float **Zintf_new;    /*Corrected interface fine grid position (nodes from 0,0)*/
  float **SI_dscm;	/*Sub-interface 2D coarse grid cells slowness model (cells from 1,1)*/
  float **I_dzcm;	/*Interface depth 2D coarse grid cells slowness model (cells from 1,1)*/
  float Rsmooth;	/*Radius for the additional smoothing of interface*/
  float Psmooth;        /*Power of exponent for the addiional smoothing of interface*/
  float Rsmooth_ds;	/*Radius for the additional smoothing of sub-interface slowness variation*/
  float Psmooth_ds;     /*Power of exponent for the addiional smoothing of sub-interface slowness variation*/
  char fi_VelMod[FLNLN], fi_fiBas[FLNLN], fi_fiTop[FLNLN];
  /*---------------------------------------------------------------------------*/
  int Sm_target=0;	/*Smooth model flag (1 - smooth target model, 0 - smooth model variance)*/
  int estimate_resolution; /*Estimate (1) or not (0) the resolution matrix - time consuming!!!*/
  int write_si_pattern;	/*Write to file or not the subinterface wavelet coverage pattern*/

  float **RM;		/*Resolution matrix*/
  float R_min, R_max;	/*Minimum and maximum elements of resolution matrix*/
  float *REdiff;	/*R-E diference by-row norm*/
  float *ND_R_norm;	/*Non-diagonal R elements by-row norm*/
  float SumL;		/*Sum of all ray length's in wavelet support*/
  int element_found;	/*Logical flag needed while setting up RHS for R computation*/

  float ***Buff_top; // N rays of Cells_hw_dw_top + Cells_reflected;
  float **Buff_2D_bot;
  float **Buff_reflect_points;	// N points   
  
  /*----------W20----------*/
  float ****UL_RV;	// Upper layer per-wavelet resolution value
  float ***SI_RV;	// Sub-interface per-wavelet resolution value
  float ***Z_RV;	// Interface per-wavelet resolution value
  
  int second_inversion = 0; //Flag of the second inversion pass with the re-arranged basis
  char resolution_value = 'p'; //Type of the resolution value used for basis rearrangement 'r' - r_ii, 's' - s_i, 'n' - nd_i, 'p' - phi_i
  float resolution_threshold; //Threshold of the corresponding resolution value
  /*----------W20----------*/
  
  char output_path[STRLN];	// Path to write all deflaut files
  
  rays_map ** Cells_hw_top;
  rays_map * Cells_2D_bott;
  rays_map ** Cells_reflected;

  reflct_points_map *W_Cells;
  reflct_points_map *WL_Cells;

  reflct_points_map:: const_iterator it_p;
  rays_map::const_iterator it;

  /*------------------------LSQR stuff-----------------------------------------*/
  lsqr_input  *in_lsqr;
  lsqr_output *out_lsqr;
  lsqr_work   *work_lsqr;
  lsqr_func   *func_lsqr;
  void        *prod_lsqr;
  block_matrix_structure CB;
  /*---------------------------------------------------------------------------*/
  char  fiRays[FLNLN],
        fiWells[FLNLN], fiCrsp[FLNLN], foZ[FLNLN], foS[FLNLN],
	fiRflRays[FLNLN], fiRflp[FLNLN];
  char buf[STRLN];

  GMT_native_header gmthead;

  FILE *fp, *fi, *fo, *fstat, *ferr, *fo_gmt, *fo_fish, *fo_resol, *fo_pattern;

#ifdef	CILKPARALLEL 
	int gl_nworkers = cilk::current_worker_count();
#else
	int gl_nworkers = 0;
#endif

  INDEXTYPE forcelogbeta = 0;

  time(&tm0);

  fp = file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fscanf(fp,"%ld",&k_gen);
  fscanf(fp,"%ld",&k_gen_v);
  fclose(fp);

  NNx=(Nx-1)/k_gen;
  NNy=(Ny-1)/k_gen;

  if (((NNx*k_gen)!=(Nx-1)) || ((NNy*k_gen)!=(Ny-1)))
  {
    fprintf(stderr,"Invalid generalization constant.\n");
    exit(EXIT_FAILURE);
  }

  N=NNy*NNx;
  dd=h*k_gen;
  N_int=(NNy-2)*(NNx-2); /*Number of internal points*/

  /*------------Read parameters----------------*/
  if (narg < 2) errquit("Insufficient parameters.\n");
  if (strcmp(argv[1],"-f")==0) /*Use parameters file with the name given as 2nd parameter*/
  {
    if (narg < 3) errquit("Parameters file name not given.\n");
    fp=file_open(argv[2],"rt");
  }
  else
  {
    if (strcmp(argv[1],"-s")==0) /*Read parameters from standart input*/
      fp=stdin;
    else
    {
      printf("Unrecognized flag %s.\n",argv[1]);
      exit(1);
    }
  }

  fscanf(fp,"%s",fiRflRays); fgets(buf,STRLN,fp);  /*1st LINE: Reflected rays traces*/
  fscanf(fp,"%s",fiRflp);    fgets(buf,STRLN,fp);  /*2nd LINE: Reflection points and derivatives file*/
  fscanf(fp,"%s",fiRays);    fgets(buf,STRLN,fp);  /*3rd LINE: Rays path's file*/
  fscanf(fp,"%s",fiCrsp);    fgets(buf,STRLN,fp);  /*4th LINE: Basement crosspoints and derivatives file*/
  fscanf(fp,"%s",foS);       fgets(buf,STRLN,fp);  /*5th LINE: Output slowness variations file*/
  fscanf(fp,"%s",foZ);       fgets(buf,STRLN,fp);  /*6th LINE: Output interface position file*/
  if (fscanf(fp,"%g%g%g",&lambda_x,&lambda_y,&lambda_z)!=3) errquit("Error 3");
  fgets(buf,STRLN,fp); /*7th LINE: Rougness regularization parameters for horizontal and vertical dir.*/
  fscanf(fp,"%g",&beta);     fgets(buf,STRLN,fp);  /*8th LINE: Global regularization parameter*/
  fscanf(fp,"%s%g",fiWells,&ww);   fgets(buf,STRLN,fp);  /*9th LINE: Wells (reflection) data file and
  							   general weighting factor*/
  fscanf(fp,"%g",&W_Nr_Min_1); fscanf(fp,"%g",&W_Nr_Min_2);
  fgets(buf,STRLN,fp);  	/*10th LINE: Minimum numbers 1 and 2 of rays per wavelet support
                                                     to include it into inversion*/
  fscanf(fp,"%g",&p2_Nr); fgets(buf,STRLN,fp);	/*11th LINE: Number of rays per wavelet support is nomalised by
  						      2^(p2_Np*l), l - wavelet order*/
  fscanf(fp,"%g",&W_Np_Min); fgets(buf,STRLN,fp);   /*12th LINE: Minimum number of refraction+reflection
  						      points per wavelet support*/
  fscanf(fp,"%g",&p2_Np); fgets(buf,STRLN,fp);	/*13th LINE: Number of points per wavelet support is nomalised by
  						      2^(p2_Np*l), l - wavelet order*/
  fscanf(fp,"%d%d",&use_down_hr,&use_up_hr);
  fgets(buf,STRLN,fp);				    /*14th LINE: Flags to use (1) or not (0) downgoing
                                                                 and upgoing headwave rays*/
  fscanf(fp,"%d%d",&use_down_rr,&use_up_rr);
  fgets(buf,STRLN,fp);				    /*15st LINE: Flags to use (1) or not (0) downgoing
                                                                 and upgoing reflected rays*/
  fscanf(fp,"%s",fi_VelMod);	fgets(buf,STRLN,fp); /*16th LINE: Fine true reference velocity model*/
  fscanf(fp,"%s",fi_fiBas);	fgets(buf,STRLN,fp); /*17th LINE: Fine reference interface position*/
  fscanf(fp,"%s",fi_fiTop);	fgets(buf,STRLN,fp); /*18th LINE: Fine layer top position*/
  fscanf(fp,"%d",&Sm_target);	fgets(buf,STRLN,fp); /*19th LINE: Smooth model flag (1 -target, 0 - variance)*/
  fscanf(fp,"%g",&Theta_1); fscanf(fp,"%g",&Theta_2);
  fgets(buf,STRLN,fp);		 /*20th LINE: Minimum angular deviation values 1 and 2 (radian)*/
  fscanf(fp,"%g",&Theta_UPL_1); fscanf(fp,"%g",&Theta_UPL_2);
  fgets(buf,STRLN,fp);		/*21st LINE: Minimum angular deviation values 1 and 2 (radian) for UPL*/
  fscanf(fp,"%g",&p2_Theta); fgets(buf,STRLN,fp);	/*22nd LINE: angular deviation is nomalised by
  						      2^(p2_Theta*l), l - wavelet order*/
  fscanf(fp,"%g%g",&Rsmooth,&Psmooth); fgets(buf,STRLN,fp); /*23rd LINE: Radius and power for additional
  								smoothing of interface*/
  fscanf(fp,"%g%g",&Rsmooth_ds,&Psmooth_ds); fgets(buf,STRLN,fp); /*24th LINE: Radius and power for additional
  								smoothing of sub-interface slowness variation*/
  fscanf(fp,"%d",&scalefunc_include); fgets(buf,STRLN,fp); /*25th LINE: Include (1) or not (0) the scaling function,
  								i.e. constant in all expansions*/
  fscanf(fp,"%d",&estimate_resolution); fgets(buf,STRLN,fp); /*26th LINE: Estimate (1) or not (0) the resolution
  								matrix: time consuming!*/
  fscanf(fp,"%d",&write_si_pattern); fgets(buf,STRLN,fp);	/*27th LINE Write to file or not the subinterface wavelet coverage pattern*/
  if (fscanf(fp,"%s",output_path)!=1) output_path[0]='\0'; fgets(buf,STRLN,fp);	/*28th LINE: Path for output files*/
  fscanf(fp,"%c",&resolution_value);	/*29th LINE: Type of the resolution value used for basis rearrangement 'r' - r_ii, 's' - s_i, 'n' - nd_i, 'p' - phi_i*/
  fscanf(fp,"%g",&resolution_threshold);	/*30th LINE: Threshold of the corresponding resolution value*/
    
  fprintf(stderr,"Parameters read.\n");

  if (strcmp(argv[1],"-f")==0) fclose(fp);
  /*--------------------------------------------------------*/
  
  second_inversion = 0; //Starting first pass of inversion

  fprintf(stderr,"Scalefunction inclusion %d\n",scalefunc_include);
  fprintf(stderr,"Resolution matrix estimation %d\n",estimate_resolution);

  //fstat=fopen(STATISTICS_FILE,"wt");
  fstat=dir_file_open(output_path,STATISTICS_FILE,"wt");
  //ctime_r(&tm0,time_string);
  //fprintf(fstat,"Program call at: "); fputs(time_string,fstat);

  /*------I STAGE: reading reference model values and averaging them onto coarse grid-------*/
  fprintf(stderr,"Allocating memory...");
  alloc_float_3d(&VFINE,Nz,Ny,Nx);
  alloc_float_3d(&dSFINE,Nz,Ny,Nx);
  alloc_float_matrix(&Zintf,Ny,Nx);
  alloc_float_matrix(&UL_ft,Ny,Nx);
  SI_dscm=matrix(1,NNy,1,NNx);
  I_dzcm=matrix(1,NNy,1,NNx);
  T=fvector(1,N);
  Z=fvector(1,N);
  fprintf(stderr,"Success.\n");

  /*------Read reference true velocity model-------------------------------------------*/
  fprintf(stderr,"Reading fine velocity model...");
  read_3D_volume(VFINE,Nz,Ny,Nx,fi_VelMod);

  /*------Read fine grids of layer top and interface*/
  fprintf(stderr,"Read surfaces...");
  read_2D_surface(UL_ft,Ny,Nx,fi_fiTop);
  read_2D_surface(Zintf,Ny,Nx,fi_fiBas);
  fprintf(stderr,"Success.\n");

  /*Averaging*/
  av2d(UL_ft,I_dzcm,Ny,Nx,k_gen);
  for (jj=1;jj<=NNy;jj++)
    for (ii=1;ii<=NNx;ii++)
      T[cln(jj,ii,NNx)]=I_dzcm[jj][ii]; /*From here now T contains topografy*/
  av2d(Zintf,I_dzcm,Ny,Nx,k_gen);
  for (jj=1;jj<=NNy;jj++)
    for (ii=1;ii<=NNx;ii++)
      Z[cln(jj,ii,NNx)]=I_dzcm[jj][ii]; /*From here now Z and I_dzcm contains TOTAL Z of interface*/
  minmax(T,1,N,&Tmin,&Tmax);
  fprintf(stderr,"Layer top position read. Tmin=%g Tmax=%g\n",
  		Tmin,Tmax);
  minmax(Z,1,N,&Zmin,&Zmax);
  fprintf(stderr,"Reference interface position read. Zmin=%g Zmax=%g\n",
  		Zmin,Zmax);

  if (eq_digital((Tmin-z0)/h,roun((Tmin-z0)/h),1e-7*h)) uk_min=roun((Tmin-z0)/h)+1;
  else uk_min=floor((Tmin-z0)/h)+1;
  if (eq_digital((Zmax-z0)/h,roun((Zmax-z0)/h),1e-7*h)) uk_max=roun((Zmax-z0)/h)+1;
  else uk_max=floor((Zmax-z0)/h);
  
  /*
  if (((Tmin-z0)/h - roun((Tmin-z0)/h)) < 1e-7*h) uk_min=floor((Tmin-z0)/h)+2;
  else uk_min=floor((Tmin-z0)/h)+1;
  if (((Zmax-z0)/h - roun((Zmax-z0)/h)) < 1e-7*h) uk_max=floor((Zmax-z0)/h)+2;
  else uk_max=floor((Zmax-z0)/h)+1;
  */
  
  NNz=(uk_max-uk_min)/k_gen_v;
  if ((uk_min+NNz*k_gen_v)<uk_max) {NNz++; uk_max=uk_min+NNz*k_gen_v;}
  NNall=N*NNz;
  ddz=h*k_gen_v;
  UL_H_min=z0+(uk_min-1)*h;

  if (uk_min<0)
  {
    fprintf(stderr,"Layer top is above the model uk_min=%ld\n",uk_min);
    exit(EXIT_FAILURE);
  }
  if (uk_max>=Nz)
  {
    fprintf(stderr,"Layer bottom is below the model uk_max=%ld\n",uk_max);
    exit(EXIT_FAILURE);
  }
  
  
  /*
  fprintf(stderr,"Reading UPL density model...");
  UL_ro=f3tensor(1,NNz,1,NNy,1,NNx);
  if (strcmp(fi_UL_ro_mod,"EMPTY")==0)
  {
    fprintf(stderr,"\nZero reference density model...\n");
    for (kk=1;kk<=NNz;kk++)
      for (jj=1;jj<=NNy;jj++)
        for (ii=1;ii<=NNx;ii++)
	  UL_ro[kk][jj][ii]=0.0;
  }
  else
  {
    fi=file_open(fi_UL_ro_mod,"rb");
    for (kk=1;kk<=NNz;kk++)
      for (jj=1;jj<=NNy;jj++)
        fread(UL_ro[kk][jj],sizeof(float),NNx,fi);
    fclose(fi);
  }
  fprintf(stderr,"Success.\n");
  */

  fprintf(stderr,"%ld cells with dimentions of %gx%g units in each layer.\n",N,dd,dd);
  fprintf(stderr,"upper layer model span: k %ld-%ld, Z %g-%g\n",uk_min,uk_max,
  		z0+h*(uk_min-1),z0+h*(uk_max-1));
  fprintf(stderr,"%ld grid layers with %g units thickness each.\n",NNz,ddz);
  fprintf(stderr,"%ld cells used for the upper physical layer approximation.\n",NNall);

  UL_dscm=f3tensor(1,NNz,1,NNy,1,NNx);

  fprintf(stderr,"check point\n");
  
  if (Sm_target==1)
  {
    av3d_ins_lay(VFINE,UL_ft,Zintf,UL_dscm,Nz,Ny,Nx,k_gen,uk_min,k_gen_v,NNz,z0,h);
    fprintf(stderr,"check point\n");
    if (av3d_sub_intf(VFINE,Zintf,SI_dscm,Nz,Ny,Nx,k_gen,z0,h)!=1)
      errquit("Error averaging sub-interface velocity model.\n");
    for (jj=1;jj<=NNy;jj++)
      for (ii=1;ii<=NNx;ii++)
      {
        if (SI_dscm[jj][ii]<=0) errquit("Zero or negative sub-intf. velocity encountered");
	if (SI_dscm[jj][ii] < DUMMY_LIM) SI_dscm[jj][ii]=1.0/SI_dscm[jj][ii];
	for (kk=1;kk<=NNz;kk++)
	{
	  if (UL_dscm[kk][jj][ii]<=0) errquit("Zero or negative sub-intf. velocity encountered");
	  //NEEDED FOR DEBUG PURPOSES
	  if ((UL_dscm[kk][jj][ii]>=3.0) && (UL_dscm[kk][jj][ii]<DUMMY_LIM))
	  {
	    fprintf(stderr,"%ld %ld %ld %g\n",kk,jj,ii,UL_dscm[kk][jj][ii]);
	  }
	  //-------------------------
	  if (UL_dscm[kk][jj][ii] < DUMMY_LIM) UL_dscm[kk][jj][ii]=1.0/UL_dscm[kk][jj][ii];
	}
      }
  }
  else
  {
    for (jj=1;jj<=NNy;jj++)
      for (ii=1;ii<=NNx;ii++)
      {
        I_dzcm[jj][ii]=SI_dscm[jj][ii]=0.0;
	for (kk=1;kk<=NNz;kk++)
	  UL_dscm[kk][jj][ii]=0.0;
      }
  }
  
  fprintf(stderr,"check point\n");

  /*------Writing out debugging info------*/
  /*
  fo=fopen("iav.2d","wb");
  for (jj=1;jj<=NNy;jj++)
    fwrite(I_dzcm[jj]+1,sizeof(float),NNx,fo);
  fclose(fo);
  fo=fopen("dsav.2d","wb");
  for (jj=1;jj<=NNy;jj++)
    fwrite(SI_dscm[jj]+1,sizeof(float),NNx,fo);
  fclose(fo);
  for (kk=1;kk<=NNz;kk++)
  {
    sprintf(buf,"ulav%0d.2d",kk);
    fo=fopen(buf,"wb");
    for (jj=1;jj<=NNy;jj++)
      fwrite(UL_dscm[kk][jj]+1,sizeof(float),NNx,fo);
    fclose(fo);
  }
  */
  /*DEBUG INFO END------------------------*/

  free_float_3d(VFINE,Nz,Ny,Nx);
  free_float_matrix(Zintf,Ny,Nx);
  free_float_matrix(UL_ft,Ny,Nx);
  /*----------------------------------------------------------------------------------------*/

  /*-------------------Read reference bulk density----------------*/
  /*
  fi=file_open(fiRO,"rb");
  fread(ro+1,sizeof(float),N,fi);
  fclose(fi);
  fprintf(stderr,"Reference bulk densities read.\n");
  */
  /*--------------------------------------------------------------*/

  /*-----Determination of wavelets image parameters and memory allocation-----*/
  W_J = floor(log(imax(NNx,NNy))/log(2));
  if ((imax(NNx,NNy)!=l2n(W_J))) W_J=W_J+1;
  W_Ni = l2n(W_J);
  W_Si = (W_Ni-NNx)/2;
  W_Sj = (W_Ni-NNy)/2;
  /*W11: all coefficients have the corresponding indices record*/
  /*reserve memory for the maximum possible number N=2^(2J)*/
  WV_1 = (Ww_type*)fmemalloc((l2n(2*W_J)),sizeof(Ww_type))-1;
  WV_Z = (Ww_type*)fmemalloc((l2n(2*W_J)),sizeof(Ww_type))-1;
  fprintf(stderr,"Model XY area:  %ldx%ld cells\n",NNx,NNy);
  fprintf(stderr,"Wavelets image: %ldx%ld pixels\n",W_Ni,W_Ni);
  fprintf(stderr,"Area wrt Image shifts: X: %ld, Y: %ld\n",W_Si,W_Sj);
  /*--------------------------------------------------------------------------*/

  /*First open reflected waves path's file to read number of rays only*/
  /*needed to allocate memory for rays mean vector arrays*/
  fi=file_open(fiRflRays,"rt");
  fscanf(fi,"%ld",&Nl);
  fprintf(stderr,"Number of reflected wave ray path's: %ld\n",Nl);
  fclose(fi);
  /*-------------------------------------------------------------------*/

  /*--------------Reading head wave rays information-------------------*/
  fi=file_open(fiRays,"rt");
  fscanf(fi,"%ld",&Nr);
  fprintf(stderr,"Number of head wave ray path's: %ld\n",Nr);

  dt=fvector(1,Nr);
  st=fvector(1,Nr);


/*
  L=matrix(1,Nr,1,N);

  for (m=1;m<=Nr;m++)
    for (i=1;i<=N;i++) L[m][i]=0.0;


  UHL=(float***)fmemalloc(Nr,sizeof(float**))-1;
  for (i=1;i<=Nr;i++) UHL[i]=matrix(1,NNz,1,N);

  RVX=(float***)fmemalloc(Nl+Nr,sizeof(float**))-1;
  for (i=1;i<=Nl+Nr;i++) RVX[i]=matrix(1,NNz,1,N);
  RVY=(float***)fmemalloc(Nl+Nr,sizeof(float**))-1;
  for (i=1;i<=Nl+Nr;i++) RVY[i]=matrix(1,NNz,1,N);
  RVZ=(float***)fmemalloc(Nl+Nr,sizeof(float**))-1;
  for (i=1;i<=Nl+Nr;i++) RVZ[i]=matrix(1,NNz,1,N);

  for (m=1;m<=Nr;m++)
    for (kk=1;kk<=NNz;kk++)
      for (i=1;i<=N;i++)
        UHL[m][kk][i]=0.0;

  for (m=1;m<=Nl+Nr;m++)
    for (kk=1;kk<=NNz;kk++)
      for (i=1;i<=N;i++)
      {
        RVX[m][kk][i]=0.0;
        RVY[m][kk][i]=0.0;
        RVZ[m][kk][i]=0.0;
      }

  Xr=matrix(1,Nr,1,N);
  Yr=matrix(1,Nr,1,N);
  for (m=1;m<=Nr;m++)
    for (i=1;i<=N;i++)
      Xr[m][i]=Yr[m][i]=0.0;


*/

  Cells_hw_top=new rays_map*[NNz];
  if (Cells_hw_top==NULL) fprintf(stderr,"Insufficient memory Cells-rows\n");
  Cells_hw_top[0]=new rays_map[NNz*N];
  if (Cells_hw_top[0]==NULL) fprintf(stderr,"Insufficient memory Cells-columns\n");
  for(i=1;i<NNz;i++)
  {
    Cells_hw_top[i]=Cells_hw_top[0]+i*(N);
  }
  fprintf(stderr,"Memory for head_wawes rays top layer  allocated\n");


  Cells_2D_bott=new rays_map[N];

  fprintf(stderr,"Memory for head_wawes rays bottom layer  allocated\n");

  Cells_reflected=new rays_map*[NNz];
  if (Cells_reflected==NULL) fprintf(stderr,"Insufficient memory Cells-rows\n");
  Cells_reflected[0]=new rays_map[NNz*N];
  if (Cells_reflected[0]==NULL) fprintf(stderr,"Insufficient memory Cells-columns\n");
  for(i=1;i<NNz;i++)
  {
    Cells_reflected[i]=Cells_reflected[0]+i*(N);
  }
  fprintf(stderr,"Memory for head_wawes rays bottom layer  allocated\n");


  fprintf(stderr,"Memory allocated. Reading rays.\n");
  for (m=1;m<=Nr;m++)
  {
    fscanf(fi,"%ld",&nray);
    if (nray!=m)
    {
      fprintf(stderr,"Incorrect value of ray number: %ld (record %ld in file %s)",nray,m,fiRays);
      exit(EXIT_FAILURE);
    }
    fscanf(fi,"%ld%ld%ld%ld%g%g%g%g%g",&nrec,&i,&j,&k,&x,&y,&z,dt+nray,st+nray);

    xp=x; yp=y; zp=z;
    fscanf(fi,"%d%ld%ld%ld%g%g%g%g%*g%d",&ncell,&i,&j,&k,&x,&y,&z,&ll,&raypos);
    while (ncell!=-1)
    {
      ii=(i-1)/k_gen+1;
      jj=(j-1)/k_gen+1;
      kk=(k-uk_min)/k_gen_v+1;
      if ((ii<1) || (jj<1) || (k<1) || (ii>NNx) || (jj>NNy) || (k>Nz))
      {
        fprintf(stderr,"Running out of model at ray %ld in cell %d. m=%ld\n",nray,ncell,m);
        exit(EXIT_FAILURE);
      }
      if (((raypos==0) || (raypos==4)) && ((kk<1) || (kk>NNz)))
      {
        fprintf(stderr,"WARNING: ray %ld with raypos %d out of layer model in cell %d\n",
			nray,raypos,ncell);
      }
      /*If ray below interface and not at the intersection point add its length*/
      bigcelln=cln(jj,ii,NNx);
      if (ll>0)
      {
       switch (raypos)
       {
         case 2:
                it=Cells_2D_bott[bigcelln-1].find(nray);
                if( it!=Cells_2D_bott[bigcelln-1].end())
                 {
                   Cells_2D_bott[bigcelln-1][nray].l+=ll;
	           Cells_2D_bott[bigcelln-1][nray].x+=(x-xp);
	           Cells_2D_bott[bigcelln-1][nray].y+=(y-yp);
                   Cells_2D_bott[bigcelln-1][nray].z= 0.0;
	        }
	        else
	        {
	           Cells_2D_bott[bigcelln-1][nray].l=ll;
	           Cells_2D_bott[bigcelln-1][nray].x=(x-xp);
	           Cells_2D_bott[bigcelln-1][nray].y=(y-yp);
                   Cells_2D_bott[bigcelln-1][nray].z=0.0;
	        }
                /*
                L[nray][bigcelln]+=ll;//нижний слой
		if (ll>0) {Xr[nray][bigcelln]+=(x-xp); Yr[nray][bigcelln]+=(y-yp);}*/
                break;
         case 4: if (use_down_hr) if ((kk>=1) && (kk<=NNz)) //верхний солой
		{
                   it=Cells_hw_top[kk-1][bigcelln-1].find(nray);
                   if( it!=Cells_hw_top[kk-1][bigcelln-1].end())
                    {
                     Cells_hw_top[kk-1][bigcelln-1][nray].l+=ll;
	             Cells_hw_top[kk-1][bigcelln-1][nray].x+=(x-xp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].y+=(y-yp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].z+=(z-zp);
	             }
	           else
	            {
	             Cells_hw_top[kk-1][bigcelln-1][nray].l=ll;
	             Cells_hw_top[kk-1][bigcelln-1][nray].x=(x-xp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].y=(y-yp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].z=(z-zp);
	            }
                   /*
		  UHL[nray][kk][bigcelln]+=ll;
		  if (ll>0) {RVX[Nl+nray][kk][bigcelln]+=(x-xp);
			     RVY[Nl+nray][kk][bigcelln]+=(y-yp);
			     RVZ[Nl+nray][kk][bigcelln]+=(z-zp);}*/
		}
                break;
	case 0: if (use_up_hr) if ((kk>=1) && (kk<=NNz)) //верхний слой
		{
                  it=Cells_hw_top[kk-1][bigcelln-1].find(nray);
                   if( it!=Cells_hw_top[kk-1][bigcelln-1].end())
                    {
                     Cells_hw_top[kk-1][bigcelln-1][nray].l+=ll;
	             Cells_hw_top[kk-1][bigcelln-1][nray].x+=(x-xp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].y+=(y-yp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].z+=(z-zp);
	             }
	           else
	            {
	             Cells_hw_top[kk-1][bigcelln-1][nray].l=ll;
	             Cells_hw_top[kk-1][bigcelln-1][nray].x=(x-xp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].y=(y-yp);
	             Cells_hw_top[kk-1][bigcelln-1][nray].z=(z-zp);
	            }
		 /* UHL[nray][kk][bigcelln]+=ll;
		  if (ll>0) {RVX[Nl+nray][kk][bigcelln]+=(x-xp);
			     RVY[Nl+nray][kk][bigcelln]+=(y-yp);
			     RVZ[Nl+nray][kk][bigcelln]+=(z-zp);}*/
		}
	 	break;
       }
      }


      xp=x; yp=y; zp=z;
      fscanf(fi,"%d%ld%ld%ld%g%g%g%g%*g%d",&ncell,&i,&j,&k,&x,&y,&z,&ll,&raypos);
    }
  }
  fclose(fi);
  for (m=1;m<=Nr;m++)
    if (st[m]==0)
    {
      fprintf(stderr,"Zero standart deviation value found for ray %ld.\n",m);
      exit(EXIT_FAILURE);
    }
  fprintf(stderr,"Rays read.\n");
  /*--------------------------------------------------------------*/

  /*----------------Reading basement crosspoint information-------*/
  fprintf(stderr,"Reading basement crosspoint information\n");

 W_Cells=new reflct_points_map [N];
 /*
 W=matrix(1,Nr,1,N);
  for (m=1;m<=Nr;m++)
    for (k=1;k<=N;k++)
      W[m][k]=0.0;*/

  fi=file_open(fiCrsp,"rt");
  fscanf(fi,"%ld",&Ncrp);
  fprintf(stderr,"Number of basement crosspoints pairs: %ld\n",Ncrp);
  for (m=1;m<=Ncrp;m++)
  {
    fscanf(fi,"%ld",&nray);
    if ((nray<1) || (nray>Nr))
    {
      printf("Incorrect value of ray number: %ld (record %ld in file %s)",nray,m,fiRays);
      exit(EXIT_FAILURE);
    }
    /*First crosspoint*/
    fscanf(fi,"%*d%*d%g%g%g%*g%g",&x,&y,&z,&ll);
    //------------------
    //ll/=2.0;
    //------------------
    ii = floor((x-x0)/dd)+1;
    jj = floor((y-y0)/dd)+1;
    if ((ii<1) || (jj<1) || (ii>NNx) || (jj>NNy))
    {
      printf("Crosspoint out of model at ray %ld.\n",nray);
      exit(EXIT_FAILURE);
    }
     W_Cells[cln(jj,ii,NNx)-1][nray]=ll;
    //W[nray][cln(jj,ii,NNx)]=ll;

    /*Second crosspoint*/
    //---------------------------------------------
    fscanf(fi,"%*d%*d%g%g%g%*g%g",&x,&y,&z,&ll);
    //THIS IS A TEST DEBUG VERSION
    //TRY TO KEEP VALUE OF DT/DZ FOR THE SECOND CROSSPOINT
    //SAME AS FOR THE FIRSTH ONE
    //TO RETURN TO "CORRECT" VERSION UNCOMMENT THE PREC. LINE
    //AND COMMENT THE NEXT LINE
    //fscanf(fi,"%*d%*d%g%g%g%*g%*g",&x,&y,&z);
    //---------------------------------------------
    ii = floor((x-x0)/dd)+1;
    jj = floor((y-y0)/dd)+1;
    if ((ii<1) || (jj<1) || (ii>NNx) || (jj>NNy))
    {
      printf("Crosspoint out of model at ray %ld.\n",nray);
      exit(EXIT_FAILURE);
    }

    //W[nray][cln(jj,ii,NNx)]=ll;
    //W[nray][cln(jj,ii,NNx)]=0.0;
  }
  fclose(fi);
  fprintf(stderr,"Crosspoints information read.\n");
  /*--------------------------------------------------------------*/

  /*--------------Reading reflected wave rays information-------------------*/
  fi=file_open(fiRflRays,"rt");
  fscanf(fi,"%ld",&Nl);
  fprintf(stderr,"Number of reflected wave ray path's: %ld\n",Nl);
  dtl=fvector(1,Nl);
  stl=fvector(1,Nl);
/*
  URL=(float***)fmemalloc(Nl,sizeof(float**))-1;
  for (i=1;i<=Nl;i++) URL[i]=matrix(1,NNz,1,N);

  for (m=1;m<=Nl;m++)
    for (kk=1;kk<=NNz;kk++)
      for (i=1;i<=N;i++)
        URL[m][kk][i]=0.0;*/

  fprintf(stderr,"Memory allocated. Reading rays.\n");
  for (m=1;m<=Nl;m++)
  {
    fscanf(fi,"%ld",&nray);
    if ((nray<1) || (nray>Nl))
    {
      fprintf(stderr,"Incorrect value of ray number: %ld (record %ld in file %s)",nray,m,fiRflRays);
      exit(EXIT_FAILURE);
    }
    fscanf(fi,"%ld%ld%ld%ld%g%g%g%g%g",&nrec,&i,&j,&k,&x,&y,&z,dtl+nray,stl+nray);

    xp=x; yp=y; zp=z;
    fscanf(fi,"%d%ld%ld%ld%g%g%g%g%*g%d",&ncell,&i,&j,&k,&x,&y,&z,&ll,&raypos);
    while (ncell!=-1)
    {
      ii=(i-1)/k_gen+1;
      jj=(j-1)/k_gen+1;
      kk=(k-uk_min)/k_gen_v+1;
      if ((ii<1) || (jj<1) || (kk<1) || (ii>NNx) || (jj>NNy) || (kk>NNz))
      {
        fprintf(stderr,"Running out of model at reflected ray %ld in cell %d. m=%ld. i=%ld j=%ld k=%ld\n"
        	      ,nray,ncell,m,i,j,k);
        exit(EXIT_FAILURE);
      }
      bigcelln=cln(jj,ii,NNx);
      switch (raypos)
      {
       	case 0: if (use_down_rr && ll> 0)
       	        {
                   it=Cells_reflected[kk-1][bigcelln-1].find(nray);
                   if( it!=Cells_reflected[kk-1][bigcelln-1].end())
                    {
                     Cells_reflected[kk-1][bigcelln-1][nray].l+=ll;
	             Cells_reflected[kk-1][bigcelln-1][nray].x+=(x-xp);
	             Cells_reflected[kk-1][bigcelln-1][nray].y+=(y-yp);
	             Cells_reflected[kk-1][bigcelln-1][nray].z+=(z-zp);
	             }
	           else
	            {
	             Cells_reflected[kk-1][bigcelln-1][nray].l=ll;
	             Cells_reflected[kk-1][bigcelln-1][nray].x=(x-xp);
	             Cells_reflected[kk-1][bigcelln-1][nray].y=(y-yp);
	             Cells_reflected[kk-1][bigcelln-1][nray].z=(z-zp);
                    }

       		 /* URL[nray][kk][bigcelln]+=ll;
       		  if (ll>0) {RVX[nray][kk][bigcelln]+=(x-xp);
			     RVY[nray][kk][bigcelln]+=(y-yp);
			     RVZ[nray][kk][bigcelln]+=(z-zp);}*/
       		}
       		break;
	case 1: if (use_up_rr && ll>0)
		{
		  it=Cells_reflected[kk-1][bigcelln-1].find(nray);
                  if( it!=Cells_reflected[kk-1][bigcelln-1].end())
                    {
                     Cells_reflected[kk-1][bigcelln-1][nray].l+=ll;
	             Cells_reflected[kk-1][bigcelln-1][nray].x+=(x-xp);
	             Cells_reflected[kk-1][bigcelln-1][nray].y+=(y-yp);
	             Cells_reflected[kk-1][bigcelln-1][nray].z+=(z-zp);
	             }
	           else
	            {
	             Cells_reflected[kk-1][bigcelln-1][nray].l=ll;
	             Cells_reflected[kk-1][bigcelln-1][nray].x=(x-xp);
	             Cells_reflected[kk-1][bigcelln-1][nray].y=(y-yp);
	             Cells_reflected[kk-1][bigcelln-1][nray].z=(z-zp);
                    }
       		  /*URL[nray][kk][bigcelln]+=ll;
       		  if (ll>0) {RVX[nray][kk][bigcelln]+=(x-xp);
			     RVY[nray][kk][bigcelln]+=(y-yp);
			     RVZ[nray][kk][bigcelln]+=(z-zp);}*/
       		}
       		break;
      }
      xp=x; yp=y; zp=z;
      fscanf(fi,"%d%ld%ld%ld%g%g%g%g%*g%d",&ncell,&i,&j,&k,&x,&y,&z,&ll,&raypos);
    }
  }
  fclose(fi);
  for (m=1;m<=Nl;m++)
    if (stl[m]==0)
    {
      fprintf(stderr,"Zero standart deviation value found for reflected ray %ld.\n",m);
      exit(EXIT_FAILURE);
    }
  fprintf(stderr,"Reflected rays traces read.\n");
  /*--------------------------------------------------------------*/

  alloc_float_3d(&Buff_top, NNz, NNy, NNx);
 

  for(kk=0; kk<NNz; kk++)
   for(jj=0; jj<NNy; jj++)
    for(ii=0; ii<NNx; ii++)
     {      
      Buff_top[kk][jj][ii]=0.0;
     } 
 
  for(kk=0; kk<NNz; kk++)
   for(jj=0; jj<NNy; jj++)
    for(ii=0; ii<NNx; ii++)
     {
       bigcelln=cln(jj+1,ii+1,NNx);
       Buff_top[kk][jj][ii]=Cells_hw_top[kk][bigcelln-1].size()+Cells_reflected[kk][bigcelln-1].size();     
     } 
 
  sprintf(buf,"top_layer_ray_coverage.3dc");
  write_3D_volume(Buff_top, NNz, NNy,  NNx, output_path, buf);
  free_float_3d(Buff_top, NNz, NNy, NNx);

  alloc_float_matrix(&Buff_2D_bot,NNy, NNx); 
 
  for(jj=0; jj<NNy; jj++)
   for(ii=0; ii<NNx; ii++)
    {    
     Buff_2D_bot[jj][ii]=0;
    } 

 
  for(jj=0; jj<NNy; jj++)
   for(ii=0; ii<NNx; ii++)
    {
     bigcelln=cln(jj+1,ii+1,NNx);
     Buff_2D_bot[jj][ii]=Cells_2D_bott[bigcelln-1].size();
    } 

  sprintf(buf,"Bot_layer_ray_coverage.2dc");
  write_2D_surface(Buff_2D_bot, NNy, NNx, output_path, buf);
  free_float_matrix(Buff_2D_bot, NNy, NNx);

   

  /*----------------Reading reflection points information-------*/
  //fo=fopen("rfl_pnt_info.XYll","wt");
  fo=dir_file_open(output_path,"rfl_pnt_info.XYll","wt");
  fi=fopen(fiRflp,"rt");
  fscanf(fi,"%ld",&i);
  if (i!=Nl)
  {
    fprintf(stderr,"Number of reflection points (%ld) is not equal to number of reflected rays (%ld).\n",
    	    i,Nl);
    exit(EXIT_FAILURE);
  }

  if (Nl!=0)
  {

   WL_Cells = new reflct_points_map [N];
  // fprintf(stderr,"WL_Cells[0].size() %ld\n",WL_Cells[0].size());

   /*WL=matrix(1,Nl,1,N);
   for (m=1;m<=Nl;m++)
     for (k=1;k<=N;k++)
       WL[m][k]=0.0;*/

   fprintf(stderr,"Number of reflection points: %ld\n",Nl);
   for (m=1;m<=Nl;m++)
   {
    fscanf(fi,"%ld",&nray);
    if ((nray<1) || (nray>Nl))
    {
      printf("Incorrect value of ray number: %ld (record %ld in file %s)",nray,m,fiRflp);
      exit(EXIT_FAILURE);
    }
    fscanf(fi,"%g%g",dtl+nray,stl+nray);
    fscanf(fi,"%*d%*d%g%g%g%*g%*g%g",&x,&y,&z,&ll);
    ii = floor((x-x0)/dd)+1;
    jj = floor((y-y0)/dd)+1;
    if ((ii<1) || (jj<1) || (ii>NNx) || (jj>NNy))
    {
      printf("Reflection point out of model at ray %ld.\n",nray);
      exit(EXIT_FAILURE);
    }
     WL_Cells[cln(jj,ii,NNx)-1][nray]=ll;
     fprintf(fo,"%g %g %g\n",x0+((float)ii-0.5)*dd,y0+((float)jj-0.5)*dd,dtl[nray]/ll);
   // WL[nray][cln(jj,ii,NNx)]=ll;
   }
   fclose(fi);
   fclose(fo);
   fprintf(stderr,"Reflection points information read.\n");
  }


  alloc_float_matrix(&Buff_reflect_points, NNy,NNx);

  for(jj=0; jj<NNy; jj++)
    for(ii=0; ii<NNx; ii++)
     {    
       Buff_reflect_points[jj][ii]=0;
     } 
 
  for(jj=0; jj<NNy; jj++)
   for(ii=0; ii<NNx; ii++)
    {
     bigcelln=cln(jj+1,ii+1,NNx);
     //fprintf(stderr,"%d %d %d %d\n",jj,ii,bigcelln-1,N);
     Buff_reflect_points[jj][ii]=W_Cells[bigcelln-1].size();
    } 
  if (Nl!=0)
   {
    for(jj=0; jj<NNy; jj++)
     for(ii=0; ii<NNx; ii++)
      {
       bigcelln=cln(jj+1,ii+1,NNx);
       //fprintf(stderr,"%d %d %d %d\n",jj,ii,bigcelln-1,N);
       Buff_reflect_points[jj][ii]+=WL_Cells[bigcelln-1].size();
      } 
   }
  

  sprintf(buf,"Reflected_points_ray_coverage.2d");
  write_2D_surface(Buff_reflect_points, NNy, NNx, output_path, buf);
  free_float_matrix(Buff_reflect_points, NNy,NNx); 

 

  /*--------------------------------------------------------------*/

  wRVX=fvector(1,Nl+Nr);
  wRVY=fvector(1,Nl+Nr);
  wRVZ=fvector(1,Nl+Nr);
  ray_included=ivector(1,Nl+Nr);
//free_ivector(ray_includede,1,Nl+Nr);

  /*FILE FOR DEBUG INFO*/
  //fo_fish=file_open(FISHER_UPL_FILE,"wt");
  fo_fish=dir_file_open(output_path,FISHER_UPL_FILE,"wt");

  /*---------------Arranging 3D upper physical layer model--------*/
  GL=(Grid_Layer_type*)fmemalloc(NNz,sizeof(Grid_Layer_type))-1;
  for (kk=1;kk<=NNz;kk++)
  {

    GL[kk].h_top=UL_H_min+(float)(kk-1)*ddz;
    GL[kk].h_bot=UL_H_min+(float)(kk)*ddz;
    GL[kk].jj_min=NNy; GL[kk].jj_max=1;
    GL[kk].ii_min=NNx; GL[kk].ii_max=1;
    GL[kk].NC=0;
    GL[kk].CV=(long*)fmemalloc(N,sizeof(long))-1;
    for (jj=1;jj<=NNy;jj++)
      for (ii=1;ii<=NNx;ii++)
      {
        m=cln(jj,ii,NNx);
	if ((T[m]<GL[kk].h_bot) && (Z[m]>GL[kk].h_top))
	{
	  if (ii<GL[kk].ii_min) GL[kk].ii_min=ii;
	  if (ii>GL[kk].ii_max) GL[kk].ii_max=ii;
	  if (jj<GL[kk].jj_min) GL[kk].jj_min=jj;
	  if (jj>GL[kk].jj_max) GL[kk].jj_max=jj;
	  GL[kk].NC++;
	  GL[kk].CV[GL[kk].NC]=m;
	}
      }
    if ((GL[kk].jj_max<GL[kk].jj_min) || (GL[kk].ii_max<GL[kk].ii_min))
    {
      fprintf(stderr,"Layer %ld (depth range %g-%g) is empty\n",kk,GL[kk].h_top,GL[kk].h_bot);
      GL[kk].ny=GL[kk].nx=0;
      GL[kk].WI_J=0; GL[kk].WI_N=0;
    }
    else
    {
      GL[kk].ny=GL[kk].jj_max-GL[kk].jj_min+1;
      GL[kk].nx=GL[kk].ii_max-GL[kk].ii_min+1;

      GL[kk].WI_J = floor(log(imax(GL[kk].nx,GL[kk].ny))/log(2));
      if ((imax(GL[kk].nx,GL[kk].ny)!=l2n(GL[kk].WI_J))) GL[kk].WI_J++;
      GL[kk].WI_N = l2n(GL[kk].WI_J);
      GL[kk].W_Si = (GL[kk].WI_N-GL[kk].nx)/2;
      GL[kk].W_Sj = (GL[kk].WI_N-GL[kk].ny)/2;
      /*W11: all coefficients have the corresponding indices record*/
      /*reserve memory for the maximum possible number N=2^(2J)*/
      GL[kk].WV = (Ww_type*)fmemalloc((l2n(2*GL[kk].WI_J)),sizeof(Ww_type))-1;

      fprintf(stderr,"Layer %ld model XY area:  %ldx%ld cells\n",kk,GL[kk].nx,GL[kk].ny);
      fprintf(stderr,"Area position: X: %ld - %ld, Y: %ld - %ld\n",GL[kk].ii_min,GL[kk].ii_max,
                                                               GL[kk].jj_min,GL[kk].jj_max);
      fprintf(stderr,"Wavelets image: %ldx%ld pixels\n",GL[kk].WI_N,GL[kk].WI_N);
      fprintf(stderr,"Area wrt Image shifts: X: %ld, Y: %ld\n",GL[kk].W_Si,GL[kk].W_Sj);
    }

    /*------------------Now arranging wawelets basis-----------------*/


    fprintf(stderr,"Now calculating ray coverage for wavelets, layer %ld...",kk);
    if (GL[kk].WI_N!=0)
    {
     /*W11: no a priory ideas on whether to inlude scaling function and wavelet with l=J*/
     /*
     GL[kk].w_N=1;
     GL[kk].WV[GL[kk].w_N].l=GL[kk].WI_J;
     GL[kk].WV[GL[kk].w_N].m=0;
     GL[kk].WV[GL[kk].w_N].n=0;
     */
     GL[kk].w_N=GL[kk].w_ncf=0; //Nothing may be included

     for (wl=GL[kk].WI_J;wl>=1;wl--)
      for (wm=0;wm<=l2n(GL[kk].WI_J-wl)-1;wm++)
        for (wn=0;wn<=l2n(GL[kk].WI_J-wl)-1;wn++)
        {
          W_Nr=0;
          SumL=0.0;
	  haar_supp(wl,wm,wn,&wimin,&wimax,&wjmin,&wjmax);

	  //fprintf(stderr,"%d %d %d\n",wl,wm,wn);
	  /*Calculating number of head wave rays per wavelet support*/
	  for (nray=1;nray<=Nr;nray++)
	  {
	    //fprintf(stderr,"Q%d\n",nray);
	    wRVX[Nl+nray]=0.0;
	    wRVY[Nl+nray]=0.0;
	    wRVZ[Nl+nray]=0.0;
	    ray_included[Nl+nray]=0;
          }



             wj=wjmin;
	    do
	    {
	      jj=(wj-GL[kk].W_Sj);
	      wi=wimin;
	      do
	      {
	        ii=(wi-GL[kk].W_Si);
	        if ((ii>=1) && (ii<=GL[kk].nx) && (jj>=1) && (jj<=GL[kk].ny))
	        {
	          //BUG230608-1a the next line was:  bigcelln=cln(jj,ii,NNx);
	          bigcelln=cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx);
                  for (it=Cells_hw_top[kk-1][bigcelln-1].begin(); it!=Cells_hw_top[kk-1][bigcelln-1].end(); it++) // for rays in cell
                        {
			  wRVX[it->first]+=it->second.x;
                          wRVY[it->first]+=it->second.y;
                          wRVZ[it->first]+=it->second.z;
                          SumL+=it->second.l;
			  ray_included[it->first]=1;
			}

	       /*   if (UHL[nray][kk][bigcelln] > 1e-7*h)
		  {
		    if (!ray_included) {W_Nr++; ray_included=1;}
		    wRVX[Nl+nray]+=RVX[Nl+nray][kk][bigcelln];//головные волны верхний слой hw_top
                    wRVY[Nl+nray]+=RVY[Nl+nray][kk][bigcelln];
                    wRVZ[Nl+nray]+=RVZ[Nl+nray][kk][bigcelln];
                    SumL+=UHL[nray][kk][bigcelln];
		  }*/


		}
		wi++;
	      }
	      while ((wi<=wimax));
	      wj++;
	    }
	    while ((wj<=wjmax));
           for (nray=1;nray<=Nr;nray++)
	    {
	    x=sqrt(pow((double)wRVX[Nl+nray],(double)2.0)+pow((double)wRVY[Nl+nray],(double)2.0)+pow((double)wRVZ[Nl+nray],(double)2.0));
            if (x>0) {wRVX[Nl+nray]/=x; wRVY[Nl+nray]/=x; wRVZ[Nl+nray]/=x;}
	    }
	  /*-------------------------------------------------------------*/

	  //fprintf(stderr,"Q\n");
	  /*Calculating number of reflected wave rays per wavelet support*/
	  for (nray=1;nray<=Nl;nray++)
	  {
	    wRVX[nray]=0.0;
	    wRVY[nray]=0.0;
	    wRVZ[nray]=0.0;
	    ray_included[nray]=0;
          }

            wj=wjmin;
	    do
	    {
	      jj=(wj-GL[kk].W_Sj);
	      wi=wimin;
	      do
	      {
	        ii=(wi-GL[kk].W_Si);
	        if ((ii>=1) && (ii<=GL[kk].nx) && (jj>=1) && (jj<=GL[kk].ny))
	        {
	          //BUG230608-1b the next line was:  bigcelln=cln(jj,ii,NNx);
	          bigcelln=cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx);

                   for (it=Cells_reflected[kk-1][bigcelln-1].begin(); it!=Cells_reflected[kk-1][bigcelln-1].end(); it++) // for rays in cell
                        {
			  wRVX[it->first]+=it->second.x;
                          wRVY[it->first]+=it->second.y;
                          wRVZ[it->first]+=it->second.z;
                          SumL+=it->second.l;
			  ray_included[it->first]=1;
			}

		   /* if (!ray_included) {W_Nr++; ray_included=1;}
		    wRVX[nray]+=RVX[nray][kk][bigcelln]; //отраженные волны
                    wRVY[nray]+=RVY[nray][kk][bigcelln];
                    wRVZ[nray]+=RVZ[nray][kk][bigcelln];
                    SumL+=URL[nray][kk][bigcelln];*/

		}
		wi++;
	      }
	      while ((wi<=wimax));
	      wj++;
	    }
	    while ((wj<=wjmax));

           for (nray=1;nray<=Nl;nray++)
          {
	    x=sqrt(pow((double)wRVX[nray],(double)2.0)+pow((double)wRVY[nray],(double)2.0)+pow((double)wRVZ[nray],(double)2.0));
            if (x>0) {wRVX[nray]/=x; wRVY[nray]/=x; wRVZ[nray]/=x;}
	  }
	  /*-------------------------------------------------------------*/

	  //fprintf(stderr,"Computing statistics...\n");
	  /*Computing Fisher statistics for ray vectors found in wavelet support area*/
	  /*W11: semi-fisher statistics: using the abslute values of direction cosines*/
	  /*thus distinguishing between the co-linear and orthogonal ray vectors*/
	  /*Theta (z) must be between 0 and pi/4*/

            wXr_mean=wYr_mean=wZr_mean=0.0;
            for (nray=1;nray<=Nl+Nr;nray++)
            {
	      /*
	      wXr_mean+=fabs(wRVX[nray]);
	      wYr_mean+=fabs(wRVY[nray]);
	      wZr_mean+=fabs(wRVZ[nray]);
	      */
	      /*Full Fisher statistics check*/
	         if(ray_included[nray]>0) W_Nr++;
	         wXr_mean+=(wRVX[nray]);
	         wYr_mean+=(wRVY[nray]);
	         wZr_mean+=(wRVZ[nray]);
	       }
	      if (W_Nr>1)
          {
            x=sqrt(pow((double)wXr_mean,(double)2.0)+pow((double)wYr_mean,(double)2.0)+pow((double)wZr_mean,(double)2.0));
            y=(W_Nr-1)/(W_Nr-x); /*precision parameter*/
            z=acos(x/(float)W_Nr); /*angular standart deviation*/
          }
          else
          {
            x=0.0; y=0.0; z=0.0;
          }
          //if (W_Nr>=W_Nr_Min) fprintf(fo,"%d %d %d %d %g %g\n",wl,wm,wn,W_Nr,y,z);
          /*-------------------------------------------------------------------------*/


	  /*W11: all coefficients have corresponding indices*/
	  if (treppe_crit(W_Nr_Min_1,W_Nr_Min_2,Theta_UPL_1,Theta_UPL_2,
                          (float)W_Nr/(float)l2n(p2_Nr*wl),(float)z/(float)l2n(p2_Theta*wl))) {
	    GL[kk].w_N++;
	    /*The largest wavelet support also corresponds to the scaling function*/
	    if ((wl==GL[kk].WI_J) && (scalefunc_include))
	    {
	      GL[kk].w_ncf++;
	      GL[kk].WV[GL[kk].w_ncf].l=wl;
	      GL[kk].WV[GL[kk].w_ncf].m=0;
	      GL[kk].WV[GL[kk].w_ncf].n=0;
	      GL[kk].WV[GL[kk].w_ncf].k=0; /*This stays for the scaling function*/
	      fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",GL[kk].w_ncf,wl,wm,wn,0,W_Nr,SumL,y,z);
	    }

	    GL[kk].WV[GL[kk].w_ncf+1].l=GL[kk].WV[GL[kk].w_ncf+2].l=GL[kk].WV[GL[kk].w_ncf+3].l=wl;
	    GL[kk].WV[GL[kk].w_ncf+1].m=GL[kk].WV[GL[kk].w_ncf+2].m=GL[kk].WV[GL[kk].w_ncf+3].m=wm;
	    GL[kk].WV[GL[kk].w_ncf+1].n=GL[kk].WV[GL[kk].w_ncf+2].n=GL[kk].WV[GL[kk].w_ncf+3].n=wn;
	    GL[kk].WV[GL[kk].w_ncf+1].k=1; GL[kk].WV[GL[kk].w_ncf+2].k=2; GL[kk].WV[GL[kk].w_ncf+3].k=3;
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",GL[kk].w_ncf+1,wl,wm,wn,1,W_Nr,SumL,y,z);
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",GL[kk].w_ncf+2,wl,wm,wn,2,W_Nr,SumL,y,z);
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",GL[kk].w_ncf+3,wl,wm,wn,3,W_Nr,SumL,y,z);
            GL[kk].w_ncf+=3;

	  }
        }
     fprintf(stderr,"Total %ld supports, %ld wavelet coefficients will be used in the inversion\n",
           GL[kk].w_N,GL[kk].w_ncf);
     fprintf(stderr,"Success.\n");
    }
    else
      NNz-=1;

    /*-------------Writing out debugging information-------------*/
    /*
    sprintf(buf,"w%0d.2d",kk);
    fo=fopen(buf,"wb");
    for (jj=1;jj<=NNy;jj++)
      for (ii=1;ii<=NNx;ii++)
      {
        wj=GL[kk].W_Sj+(jj-GL[kk].jj_min+1); wi=GL[kk].W_Si+(ii-GL[kk].ii_min+1);
        ss=0.0;
        if ((jj>=GL[kk].jj_min) && (jj<=GL[kk].jj_max) && (ii>=GL[kk].ii_min) && (ii<=GL[kk].ii_max))
	  for (wq=1;wq<=GL[kk].w_N;wq++)
            ss+=haar_w(GL[kk].WV[wq].l,GL[kk].WV[wq].m,GL[kk].WV[wq].n,3,wi,wj);
        fwrite(&ss,sizeof(float),1,fo);
      }
    fclose(fo);
    */
    /*-----------------------------------------------------------*/
    /*---------------------------------------------------------------*/
  }
  /*--------------------------------------------------------------*/

/*
  for (i=1;i<=Nl+Nr;i++) free_matrix(RVX[i],1,NNz,1,N);
  fmemfree((char*)(RVX+1));
  for (i=1;i<=Nl+Nr;i++) free_matrix(RVY[i],1,NNz,1,N);
  fmemfree((char*)(RVY+1));
  for (i=1;i<=Nl+Nr;i++) free_matrix(RVZ[i],1,NNz,1,N);
  fmemfree((char*)(RVZ+1));*/

  free_fvector(wRVZ,1,Nl+Nr);
  free_fvector(wRVY,1,Nl+Nr);
  free_fvector(wRVX,1,Nl+Nr);
  free_ivector(ray_included,1,Nl+Nr);
  fclose(fo_fish); /*DEBUG INFO FILE CLOSED*/

  /*--------------------------------------------------------------------------------------*/
  /*Counting number of internal points in the upper physical layer*/
  /*Only internal points will be used for structural similarity equations*/
  /*
  NN_int=0;
  for (kk=2;kk<=NNz-1;kk++)
  {
    for (m=1;m<=GL[kk].NC;m++)
      if (cell_is_internal(GL,NNz,NNx,kk,m)) NN_int++;
  }
  fprintf(stderr,"%d internal points in upper physical layer\n",NN_int);
  */
  /*--------------------------------------------------------------------------------------*/

  /*
  fo = fopen("L.2d","wb");
  for (m=1;m<=Nr;m++)
    fwrite(L[m]+1,sizeof(float),N,fo);
  fclose(fo);

  fo = fopen("W.2d","wb");
  for (m=1;m<=Nr;m++)
    fwrite(W[m]+1,sizeof(float),N,fo);
  fclose(fo);

  if (Nl!=0)
  {
    fo = fopen("WL.2d","wb");
    for (m=1;m<=Nr;m++)
      fwrite(WL[m]+1,sizeof(float),N,fo);
    fclose(fo);
  }
  */

  /*---------------------Reading gravity anomalies data-----------*/
  /*
  fi=file_open(fiGrav,"rt");
  fscanf(fi,"%d",&Ng);
  fprintf(stderr,"Number of gravity data points: %d\n",Ng);
  dg=vector(1,Ng);
  sg=vector(1,Ng);
  gx=vector(1,Ng);
  gy=vector(1,Ng);
  gz=vector(1,Ng);
  for (i=1;i<=Ng;i++)
  {
    fscanf(fi,"%g%g%g%g%g",gx+i,gy+i,gz+i,dg+i,sg+i);
    if (sg[i]==0)
    {
      fprintf(stderr,"Zero standart deviation value found for gravity data point %d.\n",i);
      exit(EXIT_FAILURE);
    }
  }
  fclose(fi);
  fprintf(stderr,"Gravity anomalies read.\n");
  */
  /*--------------------------------------------------------------*/

  /*---------Calculating ray coverage for wavelets----------------*/
  fprintf(stderr,"Now calculating ray coverage for wavelets...");

  wXr=fvector(1,Nr);
  wYr=fvector(1,Nr);
  ray_included=ivector(1,Nr);

  /*FILE FOR DEBUG INFO*/
  //fo_fish=file_open(FISHER_SI_FILE,"wt");
  fo_fish=dir_file_open(output_path,FISHER_SI_FILE,"wt");

  w_s1_N=W_nl_1=0; /*W11: Nothing may be included*/
  for (wl=W_J;wl>=1;wl--)
  {
    if (write_si_pattern)
    {
      sprintf(buf,"si_pattern-l%ld.xy",wl);
      //fo_pattern=file_open(buf,"wt");
     fo_pattern=dir_file_open(output_path,buf,"wt");
    }
    for (wm=0;wm<=l2n(W_J-wl)-1;wm++)
      for (wn=0;wn<=l2n(W_J-wl)-1;wn++)
      {
        W_Nr=0;
        SumL=0.0;
	haar_supp(wl,wm,wn,&wimin,&wimax,&wjmin,&wjmax);

	/*Calculating number of head wave rays per wavelet support*/
	for (nray=1;nray<=Nr;nray++)
	{
	  wXr[nray]=wYr[nray]=0.0;
          ray_included[nray]=0;
        }
	  wj=wjmin;
	  do
	  {
	    jj=(wj-W_Sj);
	    wi=wimin;
	    do
	    {
	      ii=(wi-W_Si);
	      if ((ii>=1) && (ii<=NNx) && (jj>=1) && (jj<=NNy))
	      {

                  for (it=Cells_2D_bott[cln(jj,ii,NNx)-1].begin(); it!=Cells_2D_bott[cln(jj,ii,NNx)-1].end(); it++) // for rays in cell
                    {
                         wXr[it->first]+=it->second.x;
                         wYr[it->first]+=it->second.y;
                         SumL+=it->second.l;
			 ray_included[it->first]=1;
		     }
                  /*
		  if (!ray_included) {W_Nr++; ray_included=1;}
                  wXr[nray]+=Xr[nray][cln(jj,ii,NNx)];
                  wYr[nray]+=Yr[nray][cln(jj,ii,NNx)];
                  SumL+=L[nray][cln(jj,ii,NNx)]; // нижнее 2д*/

	      }
	      wi++;
	    }
	    while ((wi<=wimax));
	    wj++;
	  }
	  while ((wj<=wjmax));

          for (nray=1;nray<=Nr;nray++)
	{
          /*Normalisation to find direction cosines of mean ray vectors*/
          x=sqrt(pow((double)wXr[nray],(double)2.0)+pow((double)wYr[nray],(double)2.0));
          if (x>0) {wXr[nray]/=x; wYr[nray]/=x;}
          if (ray_included[nray]>0)  W_Nr++;
	}


	/*-------------------------------------------------------------*/

        /*Computing Fisher statistics for ray vectors found in wavelet support area*/
        /*W11: semi-fisher statistics: using the abslute values of direction cosines*/
	/*thus distinguishing between the co-linear and orthogonal ray vectors*/
	/*Theta (z) must be between 0 and pi/4*/
        if (W_Nr>1)
        {
          wXr_mean=wYr_mean=0.0;
          for (nray=1;nray<=Nr;nray++)
          {
	    wXr_mean+=fabs(wXr[nray]);
	    wYr_mean+=fabs(wYr[nray]);
	  }
          x=sqrt(pow((double)wXr_mean,(double)2.0)+pow((double)wYr_mean,(double)2.0));
          y=(W_Nr-1)/(W_Nr-x); /*precision parameter*/
          z=acos(x/W_Nr); /*angular standart deviation*/
        }
        else
        {
          x=0.0; y=0.0; z=0.0;
        }
        /*-------------------------------------------------------------------------*/

	if (treppe_crit(W_Nr_Min_1,W_Nr_Min_2,Theta_1,Theta_2,(float)W_Nr/(float)l2n(p2_Nr*wl),
							      (float)z/(float)l2n(p2_Theta*wl))) {
	  W_nl_1++;
	  if ((wl==W_J) && (scalefunc_include)) /*The largest wavelet support also corresponds to the scaling function*/
	  {
	    w_s1_N++;
	    WV_1[w_s1_N].l=wl;
	    WV_1[w_s1_N].m=wm;
	    WV_1[w_s1_N].n=wn;
	    WV_1[w_s1_N].k=0; /*This stays for the scaling function*/
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",w_s1_N,wl,wm,wn,0,W_Nr,SumL,y,z);
	    }
	    WV_1[w_s1_N+1].l=WV_1[w_s1_N+2].l=WV_1[w_s1_N+3].l=wl;
	    WV_1[w_s1_N+1].m=WV_1[w_s1_N+2].m=WV_1[w_s1_N+3].m=wm;
	    WV_1[w_s1_N+1].n=WV_1[w_s1_N+2].n=WV_1[w_s1_N+3].n=wn;
	    WV_1[w_s1_N+1].k=1; WV_1[w_s1_N+2].k=2; WV_1[w_s1_N+3].k=3;
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",w_s1_N+1,wl,wm,wn,1,W_Nr,SumL,y,z);
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",w_s1_N+2,wl,wm,wn,2,W_Nr,SumL,y,z);
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %d %g %g %g\n",w_s1_N+3,wl,wm,wn,3,W_Nr,SumL,y,z);
            w_s1_N+=3;
	    
	    if (write_si_pattern) /*Save the pattern of the wavelet coverage in the sub-interface layer*/
	    {
	      fprintf(fo_pattern,"> -Wthick\n");
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimax-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimax-W_Si-1)*dd,y0+(wjmax-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+(wjmax-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"> -Wthin\n");
	      fprintf(fo_pattern,"%g %g\n",x0+((wimin-1+wimax)/2-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+((wimin-1+wimax)/2-W_Si-1)*dd,y0+(wjmax-W_Sj)*dd);
	      fprintf(fo_pattern,"> -Wthin\n");
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+((wjmin-1+wjmax)/2-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimax-W_Si-1)*dd,y0+((wjmin-1+wjmax)/2-W_Sj)*dd);
	    }
	}
      }
    if (write_si_pattern) fclose(fo_pattern);
  }
  fprintf(stderr,"Total %ld supports, %ld wavelet coefficients will be used in the inversion\n",
          W_nl_1,w_s1_N);

  free_ivector(ray_included,1,Nr);
  free_fvector(wYr,1,Nr);
  free_fvector(wXr,1,Nr);
 
  fclose(fo_fish);

  fprintf(stderr,"Success.\n");
  /*--------------------------------------------------------------*/

  
  //fo_fish=file_open(COV_Z_FILE,"wt");
  fo_fish=dir_file_open(output_path,COV_Z_FILE,"wt");
  
  wz_N=w_z1_N=0; /*Nothing may be included*/
  fprintf(stderr,"Now looking for the wavelets basis for the interface expantion...");
  for (wl=W_J;wl>=1;wl--)
    for (wm=0;wm<=l2n(W_J-wl)-1;wm++)
      for (wn=0;wn<=l2n(W_J-wl)-1;wn++)
      {
        W_Np_h=W_Np_r=0;
	haar_supp(wl,wm,wn,&wimin,&wimax,&wjmin,&wjmax);

	/*Allocating memory for wavelet support map*/
	wv_map=(int**)fmemalloc(wjmax-wjmin+1,sizeof(int*))-wjmin;
	for (wj=wjmin;wj<=wjmax;wj++) {
	  wv_map[wj]=(int*)fmemalloc(wimax-wimin+1,sizeof(int));
	  wv_map[wj]-=wimin;
	  for (wi=wimin;wi<=wimax;wi++) wv_map[wj][wi]=0; }
	/*-----------------------------------------*/

	/*Calculating number of crossing points per wavelet support*/


	  wj=wjmin;
	  do
	  {
	    jj=(wj-W_Sj);
	    wi=wimin;
	    do
	    {
	      ii=(wi-W_Si);
	      if ((ii>=1) && (ii<=NNx) && (jj>=1) && (jj<=NNy))
	      {
           for(it_p=W_Cells[cln(jj,ii,NNx)-1].begin(); it_p!=W_Cells[cln(jj,ii,NNx)-1].end(); it_p++) // for rays in cell
            {
			  W_Np_h++;
              wv_map[wj][wi]++;
			}
	       /* if (W[nray][cln(jj,ii,NNx)] != 0)
		{
		  if (!ray_included) {
		  W_Np_h++;
		  ray_included=1;}
		  wv_map[wj][wi]++;
		}*/

	      }
	      wi++;
	    }
	    while ((wi<=wimax));
	    wj++;
	  }
	  while ((wj<=wjmax));


	/*------------------------------------------------------------*/

	/*Calculating number of reflection points per wavelet support*/
       if (Nl>0)
       {
	  wj=wjmin;
	  do
	  {
	    jj=(wj-W_Sj);
	    wi=wimin;
	    do
	    {
	      ii=(wi-W_Si);
	      if ((ii>=1) && (ii<=NNx) && (jj>=1) && (jj<=NNy))
	      {
                for(it_p=WL_Cells[cln(jj,ii,NNx)-1].begin(); it_p!=WL_Cells[cln(jj,ii,NNx)-1].end(); it_p++) // for rays in cell
                        {
			  W_Np_r++;
		          wv_map[wj][wi]++;
			}
                /*
	        if (WL[nray][cln(jj,ii,NNx)] != 0)
		{
		 if (!ray_included) {
		  W_Np_r++;
		  ray_included=1;}
		  wv_map[wj][wi]++;
		}*/
	      }
	      wi++;
	    }
	    while ((wi<=wimax));
	    wj++;
	  }
	  while ((wj<=wjmax));
        }
	/*------------------------------------------------------------*/

	wv_cells_filled=0;
	wj_mean=wi_mean=0.0;
	for (wj=wjmin;wj<=wjmax;wj++)
	  for (wi=wimin;wi<=wimax;wi++)
	    if (wv_map[wj][wi]!=0) {
	      wj_mean+=wv_map[wj][wi]*wj;
	      wi_mean+=wv_map[wj][wi]*wi;
	      wv_cells_filled++;
	      wv_map_sum+=wv_map[wj][wi];
	    }
	wj_mean/=wv_map_sum;
	wi_mean/=wv_map_sum;
	wv_scatter=0.0;
	/*
	for (wj=wjmin;wj<=wjmax;wj++)
	  for (wi=wimin;wi<=wimax;wi++)
	  {
	    if (wv_map[wj][wi]!=0)
	     // wv_scatter+=wv_map[wj][wi]*(pow(wj-wj_mean,2.0)+pow(wi-wi_mean,2.0));
	     //wv_scatter+=wv_map[wj][wi]*((wj-wj_mean)*(wj-wj_mean)+(wi-wi_mean)*(wi-wi_mean));
	  }
	*/
	//wv_scatter=pow((double)wv_scatter/wv_map_sum,(double)0.5);

	if ((float)(W_Np_h+W_Np_r)/(float)l2n(p2_Np*wl)>=W_Np_Min) {
	  wz_N++;
	  if ((wl==W_J) && (scalefunc_include)) /*The largest wavelet support also corresponds to the scaling function*/
	  {
	    w_z1_N++;
	    WV_Z[w_z1_N].l=wl;
	    WV_Z[w_z1_N].m=wm;
	    WV_Z[w_z1_N].n=wn;
	    WV_Z[w_z1_N].k=0; /*This stays for the scaling function*/
	    /*
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g %g\n"
	                   ,w_z1_N,wl,wm,wn,0,W_Np_r+W_Np_h
	                   ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl),wv_scatter);
	    */
	    fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g\n"
	                   ,w_z1_N,wl,wm,wn,0,W_Np_r+W_Np_h
	                   ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl));
	  }
	  WV_Z[w_z1_N+1].l=WV_Z[w_z1_N+2].l=WV_Z[w_z1_N+3].l=wl;
	  WV_Z[w_z1_N+1].m=WV_Z[w_z1_N+2].m=WV_Z[w_z1_N+3].m=wm;
	  WV_Z[w_z1_N+1].n=WV_Z[w_z1_N+2].n=WV_Z[w_z1_N+3].n=wn;
	  WV_Z[w_z1_N+1].k=1; WV_Z[w_z1_N+2].k=2; WV_Z[w_z1_N+3].k=3;
	  /*
	  fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g %g\n",w_z1_N+1,wl,wm,wn,1,W_Np_r+W_Np_h
	  		 ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl),wv_scatter);
	  fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g %g\n",w_z1_N+2,wl,wm,wn,2,W_Np_r+W_Np_h
	  		 ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl),wv_scatter);
	  fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g %g\n",w_z1_N+3,wl,wm,wn,3,W_Np_r+W_Np_h
	  		 ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl),wv_scatter);
	  */
	  fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g\n",w_z1_N+1,wl,wm,wn,1,W_Np_r+W_Np_h
	  		 ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl));
	  fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g\n",w_z1_N+2,wl,wm,wn,2,W_Np_r+W_Np_h
	  		 ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl));
	  fprintf(fo_fish,"%ld %ld %ld %ld %d %ld %d %g\n",w_z1_N+3,wl,wm,wn,3,W_Np_r+W_Np_h
	  		 ,wv_cells_filled,(float)wv_cells_filled/(float)l2n(2*wl));
          w_z1_N+=3;
	}
	/*Freeing memory for wavelet support map*/
	for (wj=wjmin;wj<=wjmax;wj++)
	  fmemfree((char*)(wv_map[wj]+wimin));
	fmemfree((char*)(wv_map+wjmin));
	/*-----------------------------------------*/
      }
  fprintf(stderr,"Total %ld supports, %ld interface depth wavelet coefficients will be used\n",
          wz_N,w_z1_N);
  fprintf(stderr,"Success.\n");

  fclose(fo_fish); /*DEBUG INFO FILE CLOSED*/
  
  //fo=fopen("WL_Cells.dat","wt");
  fo=dir_file_open(output_path,"WL_Cells.dat","wt");
  if (Nl>0)
   for (i=0;i<N;i++)
   {
    pos_from_cln(&jj,&ii,i+1,NNx);
    if (!WL_Cells[cln(jj,ii,NNx)-1].empty())
    {
      fprintf(fo,"%g %g",x0+((float)ii-0.5)*dd,y0+((float)jj-0.5)*dd);
      for(it_p=WL_Cells[cln(jj,ii,NNx)-1].begin(),ss=x=0.0; it_p!=WL_Cells[cln(jj,ii,NNx)-1].end(); it_p++) // for rays in cell
      {
        ss+=it_p->second;
	x+=dtl[it_p->first];
      }
      fprintf(fo," %g\n",x/ss);
    }
   }
  fclose(fo);
  
  /*----------------------------------------------------------------------*/

  fi=file_open(fiWells,"rt");
  fscanf(fi,"%ld",&Nw);
  wz=fvector(1,Nw);
  wy=fvector(1,Nw);
  wx=fvector(1,Nw);
  ws=fvector(1,Nw);
  fprintf(stderr,"%ld interface position points will be read...\n",Nw);
  for (i=1;i<=Nw;i++)
    fscanf(fi,"%g%g%g%g",wx+i,wy+i,wz+i,ws+i);
  fclose(fi);
  fprintf(stderr,"Success.\n");

  MATRIX_LABEL:
  
  /*Number of unknowns*/
  Num_U=w_s1_N+w_z1_N;
  /*Number of equations: Number of reflected rays data + Number of head wave rays data*/
  /*+ Nw interface a priory information equations*/
  /*+ 2*N ds smoothing equations*/
  /*+ 2*N dz smoothing equations*/
  Num_E=Nl+Nr+Nw+2*N+2*N;

  /*Now calculate matrix space needed for the upper physical layer model*/
  /*Number of unknowns is increased by GL[kk].w_ncf for each grid layer*/
  for (kk=1;kk<=NNz;kk++) Num_U+=GL[kk].w_ncf;
  /*Number of equations is increased by 2*GL[kk].NC for each grid layer:*/
  /*2*GL[kk].NC smoothing equations*/
  for (kk=1;kk<=NNz;kk++) Num_E+=2*GL[kk].NC;

  /*Additional equations for smooth. constraints in vertical direction*/
  Num_SMZ_E=0;
  if ((NNz>2) && (lambda_z>0)) for (kk=2;kk<=NNz-1;kk++) Num_SMZ_E+=GL[kk].NC;
  Num_E+=Num_SMZ_E;

  fprintf(stderr,"Matrix dimentions: %ldx%ld\n",Num_E,Num_U);
//  C=matrix(1,Num_E,1,Num_U);
  B=fvector(1,Num_E);
  X=fvector(1,Num_U);
  alloc_lsqr_mem(&in_lsqr,&out_lsqr,&work_lsqr,&func_lsqr,Num_E,Num_U);
  fprintf(stderr,"Memory for inversion successfully allocated.\n");
  fprintf(stderr,"Equations system matrix %ld rows x %ld columns.\n",Num_E,Num_U);
  /*Matrix will be sparse, so first set all elements of C and B to 0*/
  for (r=1;r<=Num_E;r++)
  {
    B[r]=0.0;
    /*
    for (c=1;c<=Num_U;c++)
      C[r][c]=0.0;
    */
  }

  fprintf(stderr,"Now forming the equations system matrix...\n");

  CB.NBR=7+NNz;
  if ((NNz>2) && (lambda_z>0)) CB.NBR+=NNz-2;
  CB.NBC=2+NNz;
  alloc_bms(&CB);

  std::vector<INDEXTYPE>  row; //numbers of not zero rowes of matrix
  std::vector<INDEXTYPE>  column;//numbers of not zero columnns of matrix
  std::vector<VALUETYPE> value;//numbers of not zero velues of matrix
  matrix_count=0;//number of next element in vector

  N0=0; /*1 block row*/
  CB.Vbr[1]=N0+1;
  fprintf(stderr,"Now forming BLOCK II (WQ1)...\n");
  c0=w_s1_N;
  ss_l=fvector(1,Nl);
  ss_r=fvector(1,Nr);
//free_vector(ss_l,1,Nl);
//free_vector(ss_r,1,Nr);

  for (c=1;c<=w_z1_N;c++)
  {
    for (r=1;r<=Nl;r++) ss_l[r]=0;
    
    find_haar_supp(c,WV_Z,&wimin,&wimax,&wjmin,&wjmax);      
    for (jj=wjmin-W_Sj;jj<=wjmax-W_Sj;jj++)
      for (ii=wimin-W_Si;ii<=wimax-W_Si;ii++)
         if ((ii>=1) && (ii<=NNx) && (jj>=1) && (jj<=NNy))
          for (it_p=WL_Cells[cln(jj,ii,NNx)-1].begin(); it_p!=WL_Cells[cln(jj,ii,NNx)-1].end(); it_p++) // for rays in cell            
	    ss_l[it_p->first]+= (it_p->second)*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
	  
     for (r=1;r<=Nl;r++)           
       if(fabs(ss_l[r]/stl[r])>ASSUMED_ACCURACY)
       {
              value.push_back(ss_l[r]/stl[r]);
	      row.push_back(N0+r);
              column.push_back(c0+c);   
	      matrix_count++;    
            // C[N0+r][c0+c]=ss_l[r]/stl[r];
        }     
  }
  fprintf(stderr,"Success.\n");


  fprintf(stderr,"Now forming upper physical layer blocks of the first block row...");
  c0=w_s1_N+w_z1_N;
  for (kk=1;kk<=NNz;kk++)
  {
    fprintf(stderr,"Grid layer %ld: %ld unknowns\n",kk,GL[kk].w_ncf);
    for (c=1;c<=GL[kk].w_ncf;c++)
    {
      for (r=1;r<=Nl;r++) //отраженные волны
        ss_l[r]=0.0;
      
      find_haar_supp(c,GL[kk].WV,&wimin,&wimax,&wjmin,&wjmax);
      for (jj=wjmin-GL[kk].W_Sj;jj<=wjmax-GL[kk].W_Sj;jj++)
	for (ii=wimin-GL[kk].W_Si;ii<=wimax-GL[kk].W_Si;ii++)
	  if ((ii>=1) && (ii<=GL[kk].nx) && (jj>=1) && (jj<=GL[kk].ny))
           for (it=Cells_reflected[kk-1][cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx)-1].begin(); 
                it!=Cells_reflected[kk-1][cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx)-1].end(); it++)
                {
                  ss_l[it->first]+= (it->second.l)*ev_haar(c,GL[kk].WV,GL[kk].W_Si+ii,GL[kk].W_Sj+jj);
		}
	         //ss+=URL[r][kk][cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx)]*
	        //  ev_haar(c,GL[kk].WV,GL[kk].W_Si+ii,GL[kk].W_Sj+jj);
		
       for (r=1;r<=Nl;r++)
         if(fabs(ss_l[r]/stl[r])>ASSUMED_ACCURACY)
          {
              value.push_back(ss_l[r]/stl[r]);
	      row.push_back(N0+r);
              column.push_back(c0+c);   
	      matrix_count++;    
         //C[N0+r][c0+c]=ss_l[r]/stl[r];
          }
      }
    c0+=GL[kk].w_ncf;
  }
  fprintf(stderr,"Success.\n");

  N0=Nl; /*2 block row*/
  CB.Vbr[2]=N0+1;
  fprintf(stderr,"Now forming BLOCK III (WQ2)...\n");
  c0=0;

  for (c=1;c<=w_s1_N;c++)
  {
    for (r=1;r<=Nr;r++) // нижнее полупространство
      ss_r[r]=0;
    
    find_haar_supp(c,WV_1,&wimin,&wimax,&wjmin,&wjmax);     
    for (jj=wjmin-W_Sj;jj<=wjmax-W_Sj;jj++)
      for (ii=wimin-W_Si;ii<=wimax-W_Si;ii++)
        if ((ii>=1) && (ii<=NNx) && (jj>=1) && (jj<=NNy))
         for (it=Cells_2D_bott[cln(jj,ii,NNx)-1].begin(); it!=Cells_2D_bott[cln(jj,ii,NNx)-1].end(); it++)
            ss_r[it->first]+=(it->second.l)*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
	    //ss+=L[r][cln(jj,ii,NNx)]*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
	    
     for (r=1;r<=Nr;r++) 
       if(fabs(ss_r[r])>ASSUMED_ACCURACY)
       {
              value.push_back(ss_r[r]/st[r]);
	      row.push_back(N0+r);
              column.push_back(c0+c);   
	      matrix_count++;    
              // C[N0+r][c0+c]=ss_r[r]/st[r];
        }
    }
  fprintf(stderr,"Success.\n");

  fprintf(stderr,"Now forming BLOCK IV (WQ3)...\n");
  c0=w_s1_N;
  for (c=1;c<=w_z1_N;c++)
  {
    for (r=1;r<=Nr;r++)
      ss_r[r]=0;
    
    find_haar_supp(c,WV_Z,&wimin,&wimax,&wjmin,&wjmax);  
    for (jj=wjmin-W_Sj;jj<=wjmax-W_Sj;jj++)
      for (ii=wimin-W_Si;ii<=wimax-W_Si;ii++)
       if ((ii>=1) && (ii<=NNx) && (jj>=1) && (jj<=NNy))
       {
	 for (it_p=W_Cells[cln(jj,ii,NNx)-1].begin(); it_p!=W_Cells[cln(jj,ii,NNx)-1].end(); it_p++)
	 {
            ss_r[it_p->first]+=(it_p->second)*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
	 }
       }
	    //ss+=W[r][cln(jj,ii,NNx)]*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
	    
     for(r=1;r<=Nr;r++)
       if(fabs(ss_r[r])>ASSUMED_ACCURACY)
            {
              value.push_back(ss_r[r]/st[r]);
	      row.push_back(N0+r);
              column.push_back(c0+c);   
	      matrix_count++;     
              //C[N0+r][c0+c]=ss_r[r]/st[r];
            }
    }
  fprintf(stderr,"Success.\n");

  fprintf(stderr,"Now forming upper physical layer blocks of the second block row...");
  c0=w_s1_N+w_z1_N;
  for (kk=1;kk<=NNz;kk++)
  {
    fprintf(stderr,"Grid layer %ld: 2*%ld unknowns\n",kk,GL[kk].w_ncf);

    for (c=1;c<=GL[kk].w_ncf;c++)
    {
      for (r=1;r<=Nr;r++)//головные волны
        ss_r[r]=0.0;
      
      find_haar_supp(c,GL[kk].WV,&wimin,&wimax,&wjmin,&wjmax);
      for (jj=wjmin-GL[kk].W_Sj;jj<=wjmax-GL[kk].W_Sj;jj++)
	for (ii=wimin-GL[kk].W_Si;ii<=wimax-GL[kk].W_Si;ii++)
	  if ((ii>=1) && (ii<=GL[kk].nx) && (jj>=1) && (jj<=GL[kk].ny))
             for (it=Cells_hw_top[kk-1][cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx)-1].begin(); 
                  it!=Cells_hw_top[kk-1][cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx)-1].end(); it++)
                    ss_r[it->first]+=(it->second.l)*ev_haar(c,GL[kk].WV,GL[kk].W_Si+ii,GL[kk].W_Sj+jj);
		    
	           // ss+=UHL[r][kk][cln(jj+GL[kk].jj_min-1,ii+GL[kk].ii_min-1,NNx)]*
	          //    ev_haar(c,GL[kk].WV,GL[kk].W_Si+ii,GL[kk].W_Sj+jj);
        for (r=1;r<=Nr;r++)
          if(fabs(ss_r[r])>ASSUMED_ACCURACY)
          {
           value.push_back(ss_r[r]/st[r]);
	   row.push_back(N0+r);
           column.push_back(c0+c);   
	   matrix_count++;     
	  //C[N0+r][c0+c]=ss_r[r]/st[r];
          }
    }
    c0+=GL[kk].w_ncf;
  }
  fprintf(stderr,"Success.\n");

  free_fvector(ss_l,1,Nl);
  free_fvector(ss_r,1,Nr);

  N0+=Nr; /*3 block row*/
  CB.Vbr[3]=N0+1;
  fprintf(stderr,"Now forming BLOCK XIV (dZ a priory)...\n");// не меняю
  c0=w_s1_N;
  for (r=1;r<=Nw;r++)
  {
    i = r;
    ii = floor((wx[i]-x0)/dd)+1;
    jj = floor((wy[i]-y0)/dd)+1;
    if ((ii<1) || (jj<1) || (ii>NNx) || (jj>NNy))
    {
      printf("Interface a priory point %ld out of model.\n",i);
      exit(EXIT_FAILURE);
    }
    for (c=1;c<=w_z1_N;c++)
      if(fabs(ww)>ASSUMED_ACCURACY)
         {
           value.push_back(ww/ws[i]*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj));
	   row.push_back(N0+r);
           column.push_back(c0+c);   
	   matrix_count++;     
	  // C[N0+r][c0+c]=ww/ws[i]*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
         }
     
  }
  fprintf(stderr,"Success.\n");

  N0+=Nw; /*4 block row*/
  CB.Vbr[4]=N0+1;
  fprintf(stderr,"Now forming BLOCK VII (Sm X)\n");
  std::vector<double> *sm_vec = new vector<double>(N,0.0); //??????????? N
  c0=0;
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((ii>=2) && (ii<=NNx-1)) /*Over the internal points: second derivative square minimisation*/
    {
      for (c=1;c<=w_s1_N;c++)
        (*sm_vec)[c]=lambda_x*-2.0*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);          
	 // C[N0+r][c0+c]=lambda_x*-2.0*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);         
       
      for (c=1;c<=w_s1_N;c++)
       {
        (*sm_vec)[c]+=lambda_x*(ev_haar(c,WV_1,W_Si+ii-1,W_Sj+jj)+ev_haar(c,WV_1,W_Si+ii+1,W_Sj+jj));
        // C[N0+r][c0+c]+=lambda_x*(ev_haar(c,WV_1,W_Si+ii-1,W_Sj+jj)+ev_haar(c,WV_1,W_Si+ii+1,W_Sj+jj));
        if((*sm_vec)[c]>ASSUMED_ACCURACY)
        {
           value.push_back((*sm_vec)[c]);
	   row.push_back(N0+r);
           column.push_back(c0+c);   
	   matrix_count++;   	 
         }
       }
    }
    else
    {
      if (ii==1) /*At the edges: gradient square minimisation*/
      {
        for (c=1;c<=w_s1_N;c++)
           (*sm_vec)[c]=lambda_x*-1.0*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
         // C[N0+r][c0+c]=lambda_x*-1.0*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_s1_N;c++)
         {
           (*sm_vec)[c]+=lambda_x*(ev_haar(c,WV_1,W_Si+ii+1,W_Sj+jj));
           if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            } // C[N0+r][c0+c]
         }
      }
      else
      {
        for (c=1;c<=w_s1_N;c++)
          (*sm_vec)[c]=lambda_x*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_s1_N;c++)
         {
           (*sm_vec)[c]+=lambda_x*-1.0*(ev_haar(c,WV_1,W_Si+ii-1,W_Sj+jj));
           if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
              value.push_back((*sm_vec)[c]);
	      row.push_back(N0+r);
              column.push_back(c0+c);   
	      matrix_count++;   	 
            }
          }
          //C[N0+r][c0+c]+=lambda_x*-1.0*(ev_haar(c,WV_1,W_Si+ii-1,W_Sj+jj));
      }
    }
  }
  fprintf(stderr,"Success.\n");

  N0+=N; /*5 block row*/
  CB.Vbr[5]=N0+1;
  fprintf(stderr,"Now forming BLOCK IX (Sm Y)\n");
  c0=0;
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((jj>=2) && (jj<=NNy-1)) /*Over the internal points: second derivative square minimisation*/
    {
      for (c=1;c<=w_s1_N;c++)
        (*sm_vec)[c]=lambda_y*-2.0*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
      for (c=1;c<=w_s1_N;c++)
       {
        (*sm_vec)[c]+=lambda_y*(ev_haar(c,WV_1,W_Si+ii,W_Sj+jj-1)+ev_haar(c,WV_1,W_Si+ii,W_Sj+jj+1));
	if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }
       } 
    }
    else
    {
      if (jj==1) /*At the edges: gradient square minimisation*/
      {
        for (c=1;c<=w_s1_N;c++)
          (*sm_vec)[c]=lambda_y*-1.0*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_s1_N;c++)
         {
          (*sm_vec)[c]+=lambda_y*(ev_haar(c,WV_1,W_Si+ii,W_Sj+jj+1));
	  if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }
         }
      }
      else
      {
        for (c=1;c<=w_s1_N;c++)
          (*sm_vec)[c]=lambda_y*ev_haar(c,WV_1,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_s1_N;c++)
         {
          (*sm_vec)[c]+=lambda_y*-1.0*(ev_haar(c,WV_1,W_Si+ii,W_Sj+jj-1));
          if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }
         }
      }
    }
  }
  fprintf(stderr,"Success.\n");

  N0+=N; /*6 block row*/
  CB.Vbr[6]=N0+1;
  fprintf(stderr,"Now forming BLOCK XII (Sm X dZ)\n");
  c0=w_s1_N;
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((ii>=2) && (ii<=NNx-1)) /*Over the internal points: second derivative square minimisation*/
    {
      for (c=1;c<=w_z1_N;c++)
        (*sm_vec)[c]=lambda_x*-2.0*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
      for (c=1;c<=w_z1_N;c++)
       {
        (*sm_vec)[c]+=lambda_x*(ev_haar(c,WV_Z,W_Si+ii-1,W_Sj+jj)+ev_haar(c,WV_Z,W_Si+ii+1,W_Sj+jj));
        if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            } 
       }
    }
    else
    {
      if (ii==1)
      {
        for (c=1;c<=w_z1_N;c++)
          (*sm_vec)[c]=lambda_x*-1.0*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_z1_N;c++)
         {
          (*sm_vec)[c]+=lambda_x*(ev_haar(c,WV_Z,W_Si+ii+1,W_Sj+jj));
          if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }
         }
      }
      else
      {
        for (c=1;c<=w_z1_N;c++)
          (*sm_vec)[c]=lambda_x*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_z1_N;c++)
         { 
          (*sm_vec)[c]+=lambda_x*-1.0*(ev_haar(c,WV_Z,W_Si+ii-1,W_Sj+jj));
          if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }  
         }
      }
    }
  }
  fprintf(stderr,"Success.\n");

  N0+=N; /*7 block row*/
  CB.Vbr[7]=N0+1;
  fprintf(stderr,"Now forming BLOCK XIV (Sm Y dZ)\n");
  c0=w_s1_N;
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((jj>=2) && (jj<=NNy-1))
    {
      for (c=1;c<=w_z1_N;c++)
        (*sm_vec)[c]=lambda_y*-2.0*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
      for (c=1;c<=w_z1_N;c++)
       {
         (*sm_vec)[c]+=lambda_y*(ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj-1)+ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj+1));
         if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }
       }
    }
    else
    {
      if (jj==1)
      {
        for (c=1;c<=w_z1_N;c++)
          (*sm_vec)[c]=lambda_y*-1.0*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_z1_N;c++)
         { 
           (*sm_vec)[c]+=lambda_y*(ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj+1));
           if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }
         } 
      }
      else
      {
        for (c=1;c<=w_z1_N;c++)
          (*sm_vec)[c]=lambda_y*ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj);
        for (c=1;c<=w_z1_N;c++)
         {
          (*sm_vec)[c]+=lambda_y*-1.0*(ev_haar(c,WV_Z,W_Si+ii,W_Sj+jj-1));
          if((*sm_vec)[c]>ASSUMED_ACCURACY)
            {
             value.push_back((*sm_vec)[c]);
	     row.push_back(N0+r);
             column.push_back(c0+c);   
	     matrix_count++;   	 
            }
         }
      }
    }
  }
  fprintf(stderr,"Success.\n");

  N0+=N; /*Starting row number*/
  c0=w_s1_N+w_z1_N; /*Starting column number*/

  fprintf(stderr,"Now forming smoothing blocks corresponding to upper physical layer.\n");
  for (kk=1;kk<=NNz;kk++)
  {
    fprintf(stderr,"Layer %ld. Starting row %ld Starting column %ld\n",kk,N0,c0);

    /*Smoothing block. Sm X*/
    CB.Vbr[7+kk]=N0+1;
    for (r=1;r<=GL[kk].NC;r++)
    {
      pos_from_cln(&jj,&ii,GL[kk].CV[r],NNx);
      if ((ii>=GL[kk].ii_min+1) && (ii<=GL[kk].ii_max-1))
      {
        if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	    (UL_dscm[kk][jj][ii-1] < DUMMY_LIM) &&
	    (UL_dscm[kk][jj][ii+1] < DUMMY_LIM))
	{
	  for (c=1;c<=GL[kk].w_ncf;c++)
            (*sm_vec)[c]=lambda_x*-2.0*
	              ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+1));
          for (c=1;c<=GL[kk].w_ncf;c++)
           {
             (*sm_vec)[c]+=lambda_x*
	              (ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min),GL[kk].W_Sj+(jj-GL[kk].jj_min+1))+
		       ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+2),GL[kk].W_Sj+(jj-GL[kk].jj_min+1)));
             if((*sm_vec)[c]>ASSUMED_ACCURACY)
               {
                value.push_back((*sm_vec)[c]);
	        row.push_back(N0+r);
                column.push_back(c0+c);   
	        matrix_count++;   	 
               }
            }
	}
      }
      else
      {
        if (ii==GL[kk].ii_min)
	{
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj][ii+1] < DUMMY_LIM))
	  {
	    for (c=1;c<=GL[kk].w_ncf;c++)
              (*sm_vec)[c]=lambda_x*-1.0*
	                ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+1));
            for (c=1;c<=GL[kk].w_ncf;c++)
	     { 
              (*sm_vec)[c]+=lambda_x*
	                (ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+2),GL[kk].W_Sj+(jj-GL[kk].jj_min+1)));
              if((*sm_vec)[c]>ASSUMED_ACCURACY)
               {
                value.push_back((*sm_vec)[c]);
	        row.push_back(N0+r);
                column.push_back(c0+c);   
	        matrix_count++;   	 
               }
             }
	  }
	}
	else
	{
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj][ii-1] < DUMMY_LIM))
	  {
	    for (c=1;c<=GL[kk].w_ncf;c++)
              (*sm_vec)[c]=lambda_x*
	                ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+1));
            for (c=1;c<=GL[kk].w_ncf;c++)
             {
               (*sm_vec)[c]+=lambda_x*-1.0*
	                (ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min),GL[kk].W_Sj+(jj-GL[kk].jj_min+1)));
	       if((*sm_vec)[c]>ASSUMED_ACCURACY)
                {
                 value.push_back((*sm_vec)[c]);
	         row.push_back(N0+r);
                 column.push_back(c0+c);   
	         matrix_count++;   	 
                }		  
             }
          }
	}
      }
    }
    /*Smoothing block. Sm Y*/
    N0+=GL[kk].NC;
    for (r=1;r<=GL[kk].NC;r++)
    {
      pos_from_cln(&jj,&ii,GL[kk].CV[r],NNx);
      if ((jj>=GL[kk].jj_min+1) && (jj<=GL[kk].jj_max-1))
      {
        if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	     (UL_dscm[kk][jj-1][ii] < DUMMY_LIM) &&
	     (UL_dscm[kk][jj+1][ii] < DUMMY_LIM))
	{
	  for (c=1;c<=GL[kk].w_ncf;c++)
            (*sm_vec)[c]=lambda_y*-2.0*
	              ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+1));
          for (c=1;c<=GL[kk].w_ncf;c++)
	   {   
            (*sm_vec)[c]+=lambda_y*
	              (ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min))+
		       ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+2)));
	    if((*sm_vec)[c]>ASSUMED_ACCURACY)
             {
               value.push_back((*sm_vec)[c]);
	       row.push_back(N0+r);
               column.push_back(c0+c);   
	       matrix_count++;   	 
             }
            }
	}
      }
      else
      {
        if (jj==GL[kk].jj_min)
	{
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj+1][ii] < DUMMY_LIM))
	  {
	    for (c=1;c<=GL[kk].w_ncf;c++)
              (*sm_vec)[c]=lambda_y*-1.0*
	                ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+1));
            for (c=1;c<=GL[kk].w_ncf;c++)
             { 
              (*sm_vec)[c]+=lambda_y*
	                (ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+2)));
              if((*sm_vec)[c]>ASSUMED_ACCURACY)
               {
                value.push_back((*sm_vec)[c]);
	        row.push_back(N0+r);
                column.push_back(c0+c);   
	        matrix_count++;   	 
               }
             }
	  }
	}
	else
	{
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj-1][ii] < DUMMY_LIM))
	  {
	    for (c=1;c<=GL[kk].w_ncf;c++)
              (*sm_vec)[c]=lambda_y*
	                ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min+1));
            for (c=1;c<=GL[kk].w_ncf;c++)
             { 
              (*sm_vec)[c]+=lambda_y*-1.0*
	                (ev_haar(c,GL[kk].WV,GL[kk].W_Si+(ii-GL[kk].ii_min+1),GL[kk].W_Sj+(jj-GL[kk].jj_min)));
              if((*sm_vec)[c]>ASSUMED_ACCURACY)
               {
                value.push_back((*sm_vec)[c]);
	        row.push_back(N0+r);
                column.push_back(c0+c);   
	        matrix_count++;   	 
               }
             }
	  }
	}
      }
    }
    /*Shifting to the next block*/
    N0+=GL[kk].NC;
    if (kk<NNz)
      c0+=GL[kk].w_ncf;
  }
  fprintf(stderr,"Success.\n");

  /*Additional equations for smooth. constraints in vertical direction*/
  if ((NNz>2) && (lambda_z>0))
  {
    fprintf(stderr,"Additional equations for smoothing in vertical direction...");
    c0=w_s1_N+w_z1_N+GL[1].w_ncf; /*Starting columnn number for upper physical layer*/
    Num_SMZ_cells=0;
    for (kk=2;kk<=NNz-1;kk++) /*For internal layers only*/
    {
      for (r=1;r<=GL[kk].NC;r++)
      {
        pos_from_cln(&jj,&ii,GL[kk].CV[r],NNx);
	if ((find_long(GL[kk-1].CV,GL[kk-1].NC,GL[kk].CV[r])>0) &&
	    (find_long(GL[kk+1].CV,GL[kk+1].NC,GL[kk].CV[r])>0)) /*Cells over and below exist...*/
	  if ((UL_dscm[kk][jj][ii]<DUMMY_LIM) &&
	      (UL_dscm[kk-1][jj][ii]<DUMMY_LIM) &&
	      (UL_dscm[kk+1][jj][ii]<DUMMY_LIM)) /*...and non-empty*/
	  {
	    Num_SMZ_cells++;
            pos_from_cln(&jj,&ii,GL[kk].CV[r],NNx);
	    for (c=1;c<=GL[kk].w_ncf;c++)
	    {
              ss=lambda_z*-2.0*
	          ev_haar(c,GL[kk].WV,
		            GL[kk].W_Si+(ii-GL[kk].ii_min+1),
			    GL[kk].W_Sj+(jj-GL[kk].jj_min+1));
	      if (fabs(ss)>ASSUMED_ACCURACY)
	      {
	        value.push_back(ss);
	        row.push_back(N0+r);
                column.push_back(c0+c);   
	        matrix_count++;
	      }		
	    }
	    for (c=1;c<=GL[kk-1].w_ncf;c++)
	    {
              ss=lambda_z*
	          ev_haar(c,GL[kk-1].WV,
		            GL[kk-1].W_Si+(ii-GL[kk-1].ii_min+1),
			    GL[kk-1].W_Sj+(jj-GL[kk-1].jj_min+1));
	       if (fabs(ss)>ASSUMED_ACCURACY)
	       {
	         value.push_back(ss);
	         row.push_back(N0+r);
                 column.push_back(c0-GL[kk-1].w_ncf+c);   
	         matrix_count++;
	       }		
	    }
	    for (c=1;c<=GL[kk+1].w_ncf;c++)
	    {
              ss=lambda_z*
	          ev_haar(c,GL[kk+1].WV,
		            GL[kk+1].W_Si+(ii-GL[kk+1].ii_min+1),
			    GL[kk+1].W_Sj+(jj-GL[kk+1].jj_min+1));
	      if (fabs(ss)>ASSUMED_ACCURACY)
	      {
	         value.push_back(ss);
	         row.push_back(N0+r);
                 column.push_back(c0+GL[kk].w_ncf+c);   
	         matrix_count++;
	      }
	    }
	  }
      }

      N0+=GL[kk].NC; /*Shift down to the next grid layer equations*/
      c0+=GL[kk].w_ncf; /*Zero cell of the medium layer*/
    }
    fprintf(stderr,"Success. Total %ld cells used for smoothing in Z direction.\n",Num_SMZ_cells);
  }

  /*For historical reasons matrix indices were from (1,1)*/
  /*We shift them to be from (0,0)*/
  for (i=0;i<matrix_count;i++)
  {
    --row[i];
    --column[i];
  }

  //fo=file_open("matrix.dat","wt");
  fo=dir_file_open(output_path,"matrix.dat","wt");
  for (i=0;i<matrix_count;i++)
  {
     fprintf(fo,"%10lg %10u %10u\n",value[i],row[i],column[i]);
  } 
  fclose(fo);

  fprintf(stderr,"Equations system matrix done.\n");
  delete sm_vec;

  fprintf(stderr,"Equations system matrix done.\n");

  fprintf(stderr,"\nMatrix: %ld nonzeros (%lg %%)\n\n",matrix_count,100.0*(double)matrix_count/((double)Num_U*Num_E));

  fprintf(stderr,"Now forming right-hand side...\n");
  N0=0; /*1 block row*/
  fprintf(stderr,"Now forming BLOCK I (dt reflected)\n");
  for (r=1;r<=Nl;r++)
      B[r]=dtl[r]/stl[r];
  fprintf(stderr,"Success.\n");

  N0=Nl; /*2 block row*/
  fprintf(stderr,"Now forming BLOCK II (dt head)\n");
  for (r=1;r<=Nr;r++)
      B[N0+r]=dt[r]/st[r];
  fprintf(stderr,"Success.\n");

  N0+=Nr; /*3 block row*/
  fprintf(stderr,"Now forming BLOCK III (Z-Zref)\n");
  for (r=1;r<=Nw;r++)
  {
    i = r;
    ii = floor((wx[i]-x0)/dd)+1;
    jj = floor((wy[i]-y0)/dd)+1;
    j=cln(jj,ii,NNx);
    B[N0+r]=ww*(wz[i]-Z[j])/ws[i];
  }
  fprintf(stderr,"Success.\n");

  N0+=Nw; /*4 block row*/
  fprintf(stderr,"Now forming BLOCK IV (- Sm X ds)\n");
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((ii>=2) && (ii<=NNx-1))
      B[N0+r]=-1.0*lambda_x*(-2.0*SI_dscm[jj][ii]+SI_dscm[jj][ii-1]+SI_dscm[jj][ii+1]);
    else
      if (ii==1)
        B[N0+r]=-1.0*lambda_x*(-1.0*SI_dscm[jj][ii]+SI_dscm[jj][ii+1]);
      else
        B[N0+r]=-1.0*lambda_x*(+SI_dscm[jj][ii]-1.0*SI_dscm[jj][ii-1]);
  }
  fprintf(stderr,"Success.\n");

  N0+=N; /*5 block row*/
  fprintf(stderr,"Now forming BLOCK V (- Sm Y ds)\n");
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((jj>=2) && (jj<=NNy-1))
      B[N0+r]=-1.0*lambda_y*(-2.0*SI_dscm[jj][ii]+SI_dscm[jj-1][ii]+SI_dscm[jj+1][ii]);
    else
      if (jj==1)
        B[N0+r]=-1.0*lambda_y*(-1.0*SI_dscm[jj][ii]+SI_dscm[jj+1][ii]);
      else
        B[N0+r]=-1.0*lambda_y*(+SI_dscm[jj][ii]-1.0*SI_dscm[jj-1][ii]);
  }
  fprintf(stderr,"Success.\n");

  N0+=N; /*6 block row*/
  fprintf(stderr,"Now forming BLOCK VIII (- Sm X dZ0)\n");
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((ii>=2) && (ii<=NNx-1))
      B[N0+r]=-1.0*lambda_x*(-2.0*I_dzcm[jj][ii]+I_dzcm[jj][ii-1]+I_dzcm[jj][ii+1]);
    else
      if (ii==1)
        B[N0+r]=-1.0*lambda_x*(-1.0*I_dzcm[jj][ii]+I_dzcm[jj][ii+1]);
      else
        B[N0+r]=-1.0*lambda_x*(I_dzcm[jj][ii]-1.0*I_dzcm[jj][ii-1]);
  }
  fprintf(stderr,"Success.\n");

  N0+=N; /*7 block row*/
  fprintf(stderr,"Now forming BLOCK IX (- Sm Y dZ0)\n");
  for (r=1;r<=N;r++)
  {
    pos_from_cln(&jj,&ii,r,NNx);
    if ((jj>=2) && (jj<=NNy-1))
      B[N0+r]=-1.0*lambda_y*(-2.0*I_dzcm[jj][ii]+I_dzcm[jj-1][ii]+I_dzcm[jj+1][ii]);
    else
      if (jj==1)
        B[N0+r]=-1.0*lambda_y*(-1.0*I_dzcm[jj][ii]+I_dzcm[jj+1][ii]);
      else
        B[N0+r]=-1.0*lambda_y*(I_dzcm[jj][ii]-1.0*I_dzcm[jj-1][ii]);
  }
  fprintf(stderr,"Success.\n");
  /*--------------------------------------------------------------------------------------------------*/

  N0+=N; /*Starting row number*/
  /*--------------------------------------------------------------------------------------------------*/
  fprintf(stderr,"Now forming smoothing blocks corresponding to upper physical layer.\n");
  for (kk=1;kk<=NNz;kk++)
  {
    for (r=1;r<=GL[kk].NC;r++)
    {
      pos_from_cln(&jj,&ii,GL[kk].CV[r],NNx);
      if ((ii>=GL[kk].ii_min+1) && (ii<=GL[kk].ii_max-1))
      {
        if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	    (UL_dscm[kk][jj][ii-1] < DUMMY_LIM) &&
	    (UL_dscm[kk][jj][ii+1] < DUMMY_LIM))
	  B[N0+r]=-1.0*lambda_x*(-2.0*UL_dscm[kk][jj][ii]+UL_dscm[kk][jj][ii-1]+UL_dscm[kk][jj][ii+1]);
      }
      else
      {
        if (ii==GL[kk].ii_min) {
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj][ii+1] < DUMMY_LIM))
	    B[N0+r]=-1.0*lambda_x*(-1.0*UL_dscm[kk][jj][ii]+UL_dscm[kk][jj][ii+1]); }
	else {
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj][ii-1] < DUMMY_LIM))
	    B[N0+r]=-1.0*lambda_x*(UL_dscm[kk][jj][ii]-1.0*UL_dscm[kk][jj][ii-1]); }
      }
    }
    /*Smoothing block. Sm Y*/
    N0+=GL[kk].NC;
    for (r=1;r<=GL[kk].NC;r++)
    {
      pos_from_cln(&jj,&ii,GL[kk].CV[r],NNx);
      if ((jj>=GL[kk].jj_min+1) && (jj<=GL[kk].jj_max-1))
      {
        if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	    (UL_dscm[kk][jj-1][ii] < DUMMY_LIM) &&
	    (UL_dscm[kk][jj+1][ii] < DUMMY_LIM))
	  B[N0+r]=-1.0*lambda_y*(-2.0*UL_dscm[kk][jj][ii]+UL_dscm[kk][jj-1][ii]+UL_dscm[kk][jj+1][ii]);
      }
      else
      {
        if (jj==GL[kk].jj_min) {
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj+1][ii] < DUMMY_LIM))
	    B[N0+r]=-1.0*lambda_y*(-1.0*UL_dscm[kk][jj][ii]+UL_dscm[kk][jj+1][ii]); }
	else {
	  if ((UL_dscm[kk][jj][ii] < DUMMY_LIM) &&
	      (UL_dscm[kk][jj-1][ii] < DUMMY_LIM))
	    B[N0+r]=-1.0*lambda_y*(UL_dscm[kk][jj][ii]-1.0*UL_dscm[kk][jj-1][ii]); }
      }
    }
    N0+=GL[kk].NC;
  }
  fprintf(stderr,"Success.\n");
  /*--------------------------------------------------------------------------------------------------*/

  /*Additional equations for smooth. constraints in vertical direction*/
  if ((NNz>2) && (lambda_z>0))
  {
    fprintf(stderr,"Additional equations for smoothing in vertical direction...");
    i=0;
    for (kk=2;kk<=NNz-1;kk++) /*For internal layers only*/
    //for (kk=2;kk<=4;kk++) /*For internal layers only*/
    {
      for (r=1;r<=GL[kk].NC;r++)
      {
        pos_from_cln(&jj,&ii,GL[kk].CV[r],NNx);
	if ((find_long(GL[kk-1].CV,GL[kk-1].NC,GL[kk].CV[r])>0) &&
	    (find_long(GL[kk+1].CV,GL[kk+1].NC,GL[kk].CV[r])>0)) /*Cells over and below exist...*/
	  if ((UL_dscm[kk][jj][ii]<DUMMY_LIM) &&
	      (UL_dscm[kk-1][jj][ii]<DUMMY_LIM) &&
	      (UL_dscm[kk+1][jj][ii]<DUMMY_LIM)) /*...and non-empty*/
	  {
	    i++;
	    B[N0+r]=-1.0*lambda_z*(-2.0*UL_dscm[kk][jj][ii]+UL_dscm[kk-1][jj][ii]+UL_dscm[kk+1][jj][ii]);
	  }
      }
      N0+=GL[kk].NC; /*Shift down to the next grid layer equations*/
    }
    if (i!=Num_SMZ_cells)
    {
      fprintf(stderr,"Num_SMZ_cells test failed. Num_SMZ_cells=%ld i=%ld\n",Num_SMZ_cells,i);
      exit(EXIT_FAILURE);
    }
    fprintf(stderr,"Success.\n");
  }
  fprintf(stderr,"Right-hand side done.\n");

  fprintf(stdout,"BiCsb class constructor...\n");
  BiCsb<VALUETYPE, INDEXTYPE> bicsb((unsigned)matrix_count, (unsigned)Num_E, (unsigned)Num_U, &(row[0]), &(column[0]), &(value[0]), 
				     gl_nworkers, forcelogbeta);

  /*fo = fopen("test_row_column.txt","wt");
  for (k=0;k<matrix_count;k++)
    fprintf(fo,"%ld %ld %ld %g\n",k,row[k],column[k],value[k]);
  fclose(fo);
  */
				     
  fprintf(stderr,"Preparing data for LSQR...");
 
  func_lsqr->mat_vec_prod=aprod_bicsb;
  prod_lsqr=(void*)(&bicsb);

  in_lsqr->num_rows=Num_E;
  in_lsqr->num_cols=Num_U;
  in_lsqr->damp_val=beta;
  in_lsqr->rel_mat_err=1e-5;
  in_lsqr->rel_rhs_err=1e-5;
  in_lsqr->cond_lim=1e7;
  in_lsqr->max_iter=Num_E+Num_U+5;

  in_lsqr->lsqr_fp_out=stderr;

  for (i=1;i<=Num_E;i++)
    in_lsqr->rhs_vec->elements[i-1]=B[i];
  for (i=1;i<=Num_U;i++)
    in_lsqr->sol_vec->elements[i-1]=0.0;
  fprintf(stderr,"Success.\n");

  fprintf(stderr,"LSQR...");
  lsqr(in_lsqr,out_lsqr,work_lsqr,func_lsqr,prod_lsqr);
  fprintf(stderr,"Done.\n");

  fprintf(stderr,"LSQR diagnostics:\n");
  fprintf(stderr,"term_flag = %ld\n",out_lsqr->term_flag);
  fprintf(stderr,"Number of iterations = %ld\n",out_lsqr->num_iters);
  fprintf(stderr,"Estimate of the cond A+d*I = %e\n",out_lsqr->mat_cond_num);
  fprintf(stderr,"||Ax-b||^2 = %e\n",out_lsqr->resid_norm);
  fprintf(stderr,"||A'Ax-A'b||^2 = %e\n",out_lsqr->mat_resid_norm);
  fprintf(stderr,"||x-x0||^2 = %e\n",out_lsqr->sol_norm);

  Err=fvector(1,Num_U);

  fprintf(stderr,"Now getting results and error estimates from LSQR output vector...");
  for (i=1;i<=Num_U;i++)
  {
    X[i]=out_lsqr->sol_vec->elements[i-1]; /*X now contains solution*/
    Err[i]=out_lsqr->std_err_vec->elements[i-1]; /*Err now contains error estimates*/
  }
  fprintf(stderr,"Success.\n");

  fprintf(stderr,"Computing residual vector...");
  for (i=1;i<=Num_E;i++)
    in_lsqr->rhs_vec->elements[i-1]=0.0;
  aprod_bicsb(0,out_lsqr->sol_vec,in_lsqr->rhs_vec,prod_lsqr);
  /*in_lsqr->rhs_vec now contains modelled values of data normalized by r.m.s. errors*/
  fprintf(stderr,"Success.\n");
  
  if (second_inversion)
  {
   /*--------------Write out residual vectors and statistics---------------------*/ 
   /*----------------------------------------------------------------------------*/
     fprintf(stderr,"Saving residuals and statistics...");
   fprintf(fstat,"Residual vector norms:\n             L2          Chi^2    \n");
   rL2=rChi2=0.0; /*Total residual norms go here*/

   //fo=fopen("res_rw.dat","wt");
   fo=dir_file_open(output_path,"res_rw.dat","wt");
   y=z=0.0;
   N0=0;
   for (m=1;m<=Nl;m++)
   {
     x=in_lsqr->rhs_vec->elements[N0+m-1];
     ss=dtl[m]-x*stl[m];
     fprintf(fo,"%ld %13g %13g %13g\n",m,dtl[m],x*stl[m],ss);
     z+=ss*ss;
     y+=ss*ss/(stl[m]*stl[m]);
   }
   fclose(fo);
   rL2+=z;
   rChi2+=y;
   fprintf(fstat,"RW   : %13g %13g\n",z,y/Nl);

   //fo=fopen("res_hw.dat","wt");
   fo=dir_file_open(output_path,"res_hw.dat","wt");
   y=z=0.0;
   N0+=Nl;
   for (m=1;m<=Nr;m++)
   {
     x=in_lsqr->rhs_vec->elements[N0+m-1];
     ss=dt[m]-x*st[m];
     fprintf(fo,"%ld %13g %13g %13g\n",m,dt[m],x*st[m],ss);
     z+=ss*ss;
     y+=ss*ss/(st[m]*st[m]);
   }
   fclose(fo);
   rL2+=z;
   rChi2+=y;
   fprintf(fstat,"HW   : %13g %13g\n",z,y/Nr);
   fprintf(fstat,"Sum  : %13g %13g\n",rL2,rChi2/(Nl+Nr));

   //fo=fopen("res_wells.dat","wt");
   fo=dir_file_open(output_path,"res_wells.dat","wt");
   y=z=0.0;
   N0+=Nr;
   for (m=1;m<=Nw;m++)
   {
     x=in_lsqr->rhs_vec->elements[N0+m-1];
     ss=B[N0+m]-x;
     fprintf(fo,"%ld %13g %13g %13g\n",m,B[N0+m],x,ss);
     z+=ss*ss;
     y+=ss*ss;
   }
   fclose(fo);
   rL2+=z;
   rChi2+=y;
   fprintf(fstat,"WLS  : %13g %13g\n",z,y/Nw);
   fprintf(fstat,"Sum  : %13g %13g\n",rL2,rChi2/(Nl+Nr+Nw));

   //fo=fopen("res_sm_si_ds.dat","wt");
   fo=dir_file_open(output_path,"res_sm_si_ds.dat","wt");
   y=z=0.0;
   N0+=Nw;
   for (m=1;m<=2*N;m++)
   {
     x=in_lsqr->rhs_vec->elements[N0+m-1];
     ss=B[N0+m]-x;
     fprintf(fo,"%ld %13g %13g %13g\n",m,B[N0+m],x,ss);
     z+=ss*ss;
     y+=ss*ss;
   }
   fclose(fo);
   rL2+=z;
   rChi2+=y;
   fprintf(fstat,"SMS  : %13g %13g\n",z,y/(2*N));
   fprintf(fstat,"Sum  : %13g %13g\n",rL2,rChi2/(Nl+Nr+Nw+2*N));

   //fo=fopen("res_sm_i_dz.dat","wt");
   fo=dir_file_open(output_path,"res_sm_i_dz.dat","wt");
   y=z=0.0;
   N0+=2*N;
   for (m=1;m<=2*N;m++)
   {
     x=in_lsqr->rhs_vec->elements[N0+m-1];
     ss=B[N0+m]-x;
     fprintf(fo,"%ld %13g %13g %13g\n",m,B[N0+m],x,ss);
     z+=ss*ss;
     y+=ss*ss;
   }
   fclose(fo);
   rL2+=z;
   rChi2+=y;
   fprintf(fstat,"SMZ  : %13g %13g\n",z,y/(2*N));
   fprintf(fstat,"Sum  : %13g %13g\n",rL2,rChi2/(Nl+Nr+Nw+4*N));

   //fo=fopen("res_ul.dat","wt");
   fo=dir_file_open(output_path,"res_ul.dat","wt");
   N0+=2*N;
   for (kk=1;kk<=NNz;kk++)
   {
     fprintf(fo,"Layer %ld\n",kk);

     fprintf(fo,"SMXY\n");
     y=z=0.0;
     for (m=1;m<=2*GL[kk].NC;m++)
     {
       x=in_lsqr->rhs_vec->elements[N0+m-1];
       ss=B[N0+m]-x;
       fprintf(fo,"%ld %13g %13g %13g\n",m,B[N0+m],x,ss);
       z+=ss*ss;
       y+=ss*ss;
     }
     rL2+=z;
     rChi2+=y;
     fprintf(fstat,"SMXY : %13g %13g\n",z,y/(2*GL[kk].NC));
     fprintf(fstat,"Sum  : %13g %13g\n",rL2,rChi2/(N0+2*GL[kk].NC));
     N0+=2*GL[kk].NC;
   }

   if ((NNz>2) && (lambda_z>0))
   {
     fprintf(fo,"SMZ\n");
     y=z=0.0;
     for (m=1;m<=Num_SMZ_E;m++)
     {
       x=in_lsqr->rhs_vec->elements[N0+m-1];
       ss=B[N0+m]-x;
       fprintf(fo,"%ld %13g %13g %13g\n",m,B[N0+m],x,ss);
       z+=ss*ss;
       y+=ss*ss;
     }
     rL2+=z;
     rChi2+=y;
     fprintf(fstat,"WDC  : %13g %13g\n",z,y/Num_SMZ_E);
     fprintf(fstat,"Sum  : %13g %13g\n",rL2,rChi2/(N0+Num_SMZ_E));
   }
   fclose(fo);

   fprintf(stderr,"Success\n");
   /*!!!!!!---END: Write out residual vectors and statistics---------------!!!!!!*/
   /*----------------------------------------------------------------------------*/

   /*----------------------------------------------------------------------------*/
   /*----------------------WRITING OUT WAVELETS STRUCTURE------------------------*/
   /*----------------------------------------------------------------------------*/

   //fo=fopen(WAVELETS_STRUCTURE_FILE,"wt");
   fo=dir_file_open(output_path,WAVELETS_STRUCTURE_FILE,"wt");
   fprintf(fo,"Wavelet functions used in the model construction:\n");
   fprintf(fo,"1 - %ld: subinterface velocity coefficients\n",w_s1_N);
   for (wq=1;wq<=w_s1_N;wq++)
   {
     find_haar_indices(wq,WV_1,&wl,&wm,&wn,&wk);
     fprintf(fo,"%5d %5ld %5ld %5ld %5ld\n",wq,wl,wm,wn,wk);
   }
   N0=w_s1_N;
   fprintf(fo,"%ld - %ld: interface depth coefficients\n",N0+1,N0+w_z1_N);
   for (wq=1;wq<=w_z1_N;wq++)
   {
     find_haar_indices(wq,WV_Z,&wl,&wm,&wn,&wk);
     fprintf(fo,"%5ld %5ld %5ld %5ld %5ld\n",N0+wq,wl,wm,wn,wk);
   }
   N0+=w_z1_N;
   for (kk=1,c=0;kk<=NNz;kk++) c+=GL[kk].w_ncf;
   fprintf(fo,"%ld - %ld: UL velocity coefficients that are in %ld slices (from top to bottom)\n",N0,N0+c,NNz);
   for (kk=1;kk<=NNz;kk++)
   {
     fprintf(fo,"Slice %ld: %ld-%ld\n",kk,N0+1,N0+GL[kk].w_ncf);
     for (wq=1;wq<=GL[kk].w_ncf;wq++)
     {
       find_haar_indices(wq,GL[kk].WV,&wl,&wm,&wn,&wk);
       fprintf(fo,"%5ld %5ld %5ld %5ld %5ld\n",N0+wq,wl,wm,wn,wk);
     }
     N0+=GL[kk].w_ncf;
   }
   fclose(fo);
  }
  
  /*----------------------------------------------------------------------------*/
  /*----------------------ESTIMATING RESOLUTION MATRIX--------------------------*/
  /*----------------------------------------------------------------------------*/
  if (!second_inversion)
  {
   
   UL_RV = new float***[NNz];
   for (kk=1;kk<=NNz;kk++)
   {
     UL_RV[kk-1] = new float**[GL[kk].WI_J];
     alloc_float_3d(UL_RV+kk-1,GL[kk].WI_J,l2n(GL[kk].WI_J-wl),l2n(GL[kk].WI_J-wl));
   }
   alloc_float_3d(&SI_RV,W_J,l2n(W_J-wl),l2n(W_J-wl));
   alloc_float_3d(&Z_RV,W_J,l2n(W_J-wl),l2n(W_J-wl));
   
   for (kk=1;kk<=NNz;kk++)
   {
     for (wl=GL[kk].WI_J;wl>=1;wl--)
       for (wm=0;wm<=l2n(GL[kk].WI_J-wl)-1;wm++)
         for (wn=0;wn<=l2n(GL[kk].WI_J-wl)-1;wn++)
	   UL_RV[kk-1][wl-1][wm][wn]=DUMMY_VALUE;
   }
   for (wl=W_J;wl>=1;wl--)
     for (wm=0;wm<=l2n(W_J-wl)-1;wm++)
       for (wn=0;wn<=l2n(W_J-wl)-1;wn++)
       {
	   SI_RV[wl-1][wm][wn]=DUMMY_VALUE;
	   Z_RV[wl-1][wm][wn]=DUMMY_VALUE;
       } 
    
   /*Resolution matrix - to output in GMT format*/
   alloc_float_matrix(&RM,Num_U,Num_U);

   fprintf(stderr,"Estimating resolution matrix\n");

   /*Output of the markup to plot the resolution matrix with GMT*/
   //fo=fopen(RESOLUTION_MARKUP_FILE,"wt");
   fo=dir_file_open(output_path,RESOLUTION_MARKUP_FILE,"wt");
   N0=0;
   if (w_s1_N!=0)
   {
     fprintf(fo,"> -W3,0/0/0\n");
     fprintf(fo,"%8f %8f\n",(float)1,-(float)(N0+w_s1_N));
     fprintf(fo,"%8f %8f\n",(float)Num_U,-(float)(N0+w_s1_N));
     fprintf(fo,"> -W3,0/0/0\n");
     fprintf(fo,"%8f %8f\n",(float)N0+w_s1_N,-(float)1);
     fprintf(fo,"%8f %8f\n",(float)N0+w_s1_N,-(float)Num_U);
   }
   N0+=w_s1_N;
   if (w_z1_N!=0)
   {
     fprintf(fo,"> -W3,0/0/0\n");
     fprintf(fo,"%8f %8f\n",(float)1,-(float)(N0+w_z1_N));
     fprintf(fo,"%8f %8f\n",(float)Num_U,-(float)(N0+w_z1_N));
     fprintf(fo,"> -W3,0/0/0\n");
     fprintf(fo,"%8f %8f\n",(float)N0+w_z1_N,-(float)1);
     fprintf(fo,"%8f %8f\n",(float)N0+w_z1_N,-(float)Num_U);
   }
   N0+=w_z1_N;
   for (kk=1;kk<=NNz;kk++)
   {
    if (GL[kk].w_ncf!=0)
    {
     fprintf(fo,"> -W3,0/255/0\n");
     fprintf(fo,"%8f %8f\n",(float)1,-(float)(N0+GL[kk].w_ncf));
     fprintf(fo,"%8f %8f\n",(float)Num_U,-(float)(N0+GL[kk].w_ncf));
     fprintf(fo,"> -W3,0/255/0\n");
     fprintf(fo,"%8f %8f\n",(float)N0+GL[kk].w_ncf,-(float)1);
     fprintf(fo,"%8f %8f\n",(float)N0+GL[kk].w_ncf,-(float)Num_U);
    }
    N0+=GL[kk].w_ncf;
   }
   fclose(fo);
   /*-----------------------------------------------------------*/

   /*--Memory for the R-E by-row difference norm--*/
   REdiff = new float[Num_U];
   ND_R_norm = new float[Num_U];
   /*---------------------------------------------*/
   
   //fo=fopen(RESOLUTION_MATRIX_FILE,"wt");
   fo=dir_file_open(output_path,RESOLUTION_MATRIX_FILE,"wt");
   fprintf(fo,"Resolution matrix structure:\n");
   fprintf(fo,"Rows 1 - %ld: subinterface velocity coefficients\n",w_s1_N);
   fprintf(fo,"Rows %ld - %ld: interface depth coefficients\n",w_s1_N+1,w_s1_N+w_z1_N);
   fprintf(fo,"UL velocity coefficients are in %ld slices\n",NNz);
   for (kk=1,N0=w_s1_N+w_z1_N;kk<=NNz;kk++) {fprintf(fo,"%ld: %ld-%ld\n",kk,N0+1,N0+GL[kk].w_ncf); N0+=GL[kk].w_ncf;}
   fprintf(fo,"\n");
   for (i=1;i<=Num_U;i++) fprintf(fo,"%8ld ",i);
   fprintf(fo,"\n");
   R_max=-MAXFLOAT;
   R_min=MAXFLOAT;
  
   for (i=1;i<=Num_U;i++)
   {
    fprintf(stderr,"Row %ld of %ld\n",i,Num_U);
    fprintf(fo,"%5ld ",i);
    /*Put the i'th column of the Frechet derivatives matrix as the unknowns*/
   
    for (j=1;j<=Nl+Nr;j++)
    {
      element_found=0;
      for (k=0;k<matrix_count;k++)
      {
        if ((row[k]==(j-1)) && (column[k]==(i-1)))
	{
	  in_lsqr->rhs_vec->elements[j-1]=value[k];
	  element_found=1;
	  break;
	}
      }
      if (!element_found) in_lsqr->rhs_vec->elements[j-1]=0.0;
    }
    /*The regularisation equations corresponding RHS remain unchanged*/
    for (j=Nl+Nr+1;j<=Num_E;j++)
      in_lsqr->rhs_vec->elements[j-1]=B[j];

    for (j=1;j<=Num_U;j++)
      in_lsqr->sol_vec->elements[j-1]=0.0;

    fprintf(stderr,"LSQR...");
    lsqr(in_lsqr,out_lsqr,work_lsqr,func_lsqr,prod_lsqr);
    fprintf(stderr,"Done.\n");

    fprintf(stderr,"LSQR diagnostics:\n");
    fprintf(stderr,"term_flag = %ld\n",out_lsqr->term_flag);
    fprintf(stderr,"Number of iterations = %ld\n",out_lsqr->num_iters);
    fprintf(stderr,"Estimate of the cond A+d*I = %e\n",out_lsqr->mat_cond_num);
    fprintf(stderr,"||Ax-b||^2 = %e\n",out_lsqr->resid_norm);
    fprintf(stderr,"||A'Ax-A'b||^2 = %e\n",out_lsqr->mat_resid_norm);
    fprintf(stderr,"||x-x0||^2 = %e\n",out_lsqr->sol_norm);

    ND_R_norm[i-1]=REdiff[i-1]=0;
    for (j=1;j<=Num_U;j++) {
      fprintf(fo,"%8.4f ",out_lsqr->sol_vec->elements[j-1]);
      RM[i-1][j-1]=out_lsqr->sol_vec->elements[j-1];
      if (R_max<RM[i-1][j-1]) R_max=RM[i-1][j-1];
      if (R_min>RM[i-1][j-1]) R_min=RM[i-1][j-1];
      REdiff[i-1]+=pow(RM[i-1][j-1]-delta_kron(i-1,j-1),2.0);
      if (i-1 != j-1) ND_R_norm[i-1]+=pow((double)RM[i-1][j-1],2.0);
    }
    REdiff[i-1]=sqrt(REdiff[i-1]/(Num_U));
    ND_R_norm[i-1]=sqrt(ND_R_norm[i-1]/(Num_U-1));
    fprintf(fo,"\n");
   }

   /*Output of resolution matrix to file for GMT plotting*/
   gmthead.nx=gmthead.ny=Num_U;
   gmthead.node_offset=1;
   gmthead.x_min=0;
   gmthead.y_min=-Num_U;
   gmthead.x_max=Num_U;
   gmthead.y_max=0;
   gmthead.z_min=R_min;
   gmthead.z_max=R_max;
   gmthead.x_inc=1.0;
   gmthead.y_inc=1.0;
   gmthead.z_scale_factor=1.0;
   gmthead.z_add_offset=0.0;
   strcpy(gmthead.x_units,"\0");
   strcpy(gmthead.y_units,"\0");
   strcpy(gmthead.z_units,"\0");
   strcpy(gmthead.title,"Resolution matrix");
   strcpy(gmthead.command,argv[0]);
   strcpy(gmthead.remark,"\0");

   //fo_gmt=fopen(RESOLUTION_MATRIX_GMT_FILE,"wb");
   fo_gmt=dir_file_open(output_path,RESOLUTION_MATRIX_GMT_FILE,"wb");
   fwrite(&gmthead,sizeof(GMT_native_header),1,fo_gmt);
   for (i=0;i<=Num_U-1;i++)
     fwrite(RM[i],sizeof(float),Num_U,fo_gmt);
   fclose(fo_gmt);
   /*----------------------------------------------------*/

   /*-----Output of the R diagonal elements with the corresponding hit counts and fisher statistics-----*/
   fi=dir_file_open(output_path,FISHER_SI_FILE,"rt");
   // fo_resol=file_open(RESOLUTION_SI_FILE,"wt");
   fo_resol=dir_file_open(output_path,RESOLUTION_SI_FILE,"wt");
   N0=0;
   for (i=1;i<=w_s1_N;i++)
   {
     fscanf(fi,"%d%ld%ld%ld%ld%d%g%g%g",&wq,&wl,&wm,&wn,&wk,&W_Nr,&SumL,&y,&z);
     if (WV_1[N0+i].l!=wl || WV_1[N0+i].m!=wm || WV_1[N0+i].n!=wn || WV_1[N0+i].k!=wk)
     {
       fprintf(stderr,"Incompatible wavelet indices for wq=%ld (file/WV): l %ld/%d m %ld/%d n %ld/%d k %ld/%d\n",
       		      N0+i,wl,WV_1[N0+i].l,wm,WV_1[N0+i].m,wn,WV_1[N0+i].n,wk,WV_1[N0+i].k);
     }
     fprintf(fo_resol,"%ld %ld %ld %ld %ld %d %g %g %g %g %g %g %g\n"
     		     ,N0+i,wl,wm,wn,wk,W_Nr,y,z,REdiff[N0+i-1],RM[N0+i-1][N0+i-1],ND_R_norm[N0+i-1],RM[N0+i-1][N0+i-1]/ND_R_norm[N0+i-1],SumL);
    switch (resolution_value) 
    {
       case ('r'): SI_RV[wl-1][wm][wn]=RM[N0+i-1][N0+i-1]; break;
       case ('s'): SI_RV[wl-1][wm][wn]=REdiff[N0+i-1]; break;
       case ('n'): SI_RV[wl-1][wm][wn]=ND_R_norm[N0+i-1]; break;
       case ('p'): SI_RV[wl-1][wm][wn]=RM[N0+i-1][N0+i-1]/ND_R_norm[N0+i-1]; break;
       default:	;
     }
   }
   file_close(fo_resol);
   file_close(fi);
   N0+=w_s1_N;

   fi=dir_file_open(output_path,COV_Z_FILE,"rt");
   fo_resol=dir_file_open(output_path,RESOLUTION_Z_FILE,"wt");
   for (i=1;i<=w_z1_N;i++)
   {
     fscanf(fi,"%d%ld%ld%ld%ld%ld%ld%g",&wq,&wl,&wm,&wn,&wk,&W_Np_h,&W_Np_r,&z);
     if (WV_Z[N0+i].l!=wl || WV_Z[N0+i].m!=wm || WV_Z[N0+i].n!=wn || WV_Z[N0+i].k!=wk)
     {
       fprintf(stderr,"Incompatible wavelet indices for wq=%ld (file/WV): l %ld/%d m %ld/%d n %ld/%d k %ld/%d\n",
       		      N0+i,wl,WV_Z[N0+i].l,wm,WV_Z[N0+i].m,wn,WV_Z[N0+i].n,wk,WV_Z[N0+i].k);
     }
     fprintf(fo_resol,"%ld %ld %ld %ld %ld %d %g %g %g %g %g %g %g\n"
     		     ,N0+i,wl,wm,wn,wk,W_Nr,y,z,REdiff[N0+i-1],RM[N0+i-1][N0+i-1],ND_R_norm[N0+i-1],RM[N0+i-1][N0+i-1]/ND_R_norm[N0+i-1],SumL);
     
     switch (resolution_value) 
     {
       case ('r'): Z_RV[wl-1][wm][wn]=RM[N0+i-1][N0+i-1]; break;
       case ('s'): Z_RV[wl-1][wm][wn]=REdiff[N0+i-1]; break;
       case ('n'): SI_RV[wl-1][wm][wn]=ND_R_norm[N0+i-1]; break;
       case ('p'): Z_RV[wl-1][wm][wn]=RM[N0+i-1][N0+i-1]/ND_R_norm[N0+i-1]; break;
       default:	;
     }
   }
   file_close(fo_resol);
   file_close(fi);
   N0+=w_z1_N;

   //fo_resol=file_open(RESOLUTION_UPL_FILE,"wt");
   fo_resol=dir_file_open(output_path,RESOLUTION_UPL_FILE,"wt");
   fi=dir_file_open(output_path,FISHER_UPL_FILE,"rt");
   for (kk=1;kk<=NNz;kk++)
   {
    for (i=1;i<=GL[kk].w_ncf;i++)
    {
     fscanf(fi,"%d%ld%ld%ld%ld%d%g%g%g",&wq,&wl,&wm,&wn,&wk,&W_Nr,&SumL,&y,&z);
     if (GL[kk].WV[N0+i].l!=wl || GL[kk].WV[N0+i].m!=wm || GL[kk].WV[N0+i].n!=wn || GL[kk].WV[N0+i].k!=wk)
     {
       fprintf(stderr,"Incompatible wavelet indices for wq=%ld (file/WV): l %ld/%d m %ld/%d n %ld/%d k %ld/%d\n",
       		      N0+i,wl,GL[kk].WV[N0+i].l,wm,GL[kk].WV[N0+i].m,wn,GL[kk].WV[N0+i].n,wk,GL[kk].WV[N0+i].k);
     }
     /*
     fprintf(fo_resol,"%ld %ld %ld %ld %ld %d %g %g %g %g %g %g %g\n"
     		     ,N0+i,wl,wm,wn,wk,W_Nr,y,z,REdiff[N0+i-1],RM[N0+i-1][N0+i-1],ND_R_norm[N0+i-1],RM[N0+i-1][N0+i-1]/ND_R_norm[N0+i-1],SumL);
     */
     fprintf(fo_resol,"%ld %ld %ld %ld %ld %d %g %g %g %g %g %g %g\n"
     		     ,N0+i,wl,wm,wn,wk,W_Nr,y,z,REdiff[N0+i-1],RM[N0+i-1][N0+i-1],ND_R_norm[N0+i-1],RM[N0+i-1][N0+i-1]/REdiff[N0+i-1],SumL);
     switch (resolution_value)
     {
       case ('r'): UL_RV[kk-1][wl-1][wm][wn]=RM[N0+i-1][N0+i-1]; break;
       case ('s'): UL_RV[kk-1][wl-1][wm][wn]=REdiff[N0+i-1]; break;
       case ('n'): UL_RV[kk-1][wl-1][wm][wn]=ND_R_norm[N0+i-1]; break;
       /*case ('p'): UL_RV[kk-1][wl-1][wm][wn]=RM[N0+i-1][N0+i-1]/ND_R_norm[N0+i-1]; break;*/
       case ('p'): UL_RV[kk-1][wl-1][wm][wn]=RM[N0+i-1][N0+i-1]/REdiff[N0+i-1]; break;
       default:	;
     }
    }
    N0+=GL[kk].w_ncf;
   }
   file_close(fi);
   file_close(fo_resol);

   /*-----Output of the R diagonal elements with the corresponding hit counts and fisher statistics-----*/
   
   fclose(fo);
   fprintf(stderr,"Success");
  /*----------------------------------------------------------------------------*/
  /*-----------------END: ESTIMATING RESOLUTION MATRIX--------------------------*/
  /*----------------------------------------------------------------------------*/
  
   /*------------------Now rearranging wawelets basis according to resolution---*/
   /*------------------------estimated from resolution matrix-------------------*/
   for (kk=1;kk<=NNz;kk++)
   {
    fprintf(stderr,"Now rearranging wavelets basis, layer %ld...\n",kk);
    if (GL[kk].WI_N!=0)
    {
     GL[kk].w_N=GL[kk].w_ncf=0; //Nothing may be included

     for (wl=GL[kk].WI_J;wl>=1;wl--)
      for (wm=0;wm<=l2n(GL[kk].WI_J-wl)-1;wm++)
        for (wn=0;wn<=l2n(GL[kk].WI_J-wl)-1;wn++)
        {
	  if (resolution_test(resolution_value,resolution_threshold,UL_RV[kk-1][wl-1][wm][wn])==1) {
	    GL[kk].w_N++;
	    /*The largest wavelet support also corresponds to the scaling function*/
	    if ((wl==GL[kk].WI_J) && (scalefunc_include))
	    {
	      GL[kk].w_ncf++;
	      GL[kk].WV[GL[kk].w_ncf].l=wl;
	      GL[kk].WV[GL[kk].w_ncf].m=0;
	      GL[kk].WV[GL[kk].w_ncf].n=0;
	      GL[kk].WV[GL[kk].w_ncf].k=0; /*This stays for the scaling function*/
	    }

	    GL[kk].WV[GL[kk].w_ncf+1].l=GL[kk].WV[GL[kk].w_ncf+2].l=GL[kk].WV[GL[kk].w_ncf+3].l=wl;
	    GL[kk].WV[GL[kk].w_ncf+1].m=GL[kk].WV[GL[kk].w_ncf+2].m=GL[kk].WV[GL[kk].w_ncf+3].m=wm;
	    GL[kk].WV[GL[kk].w_ncf+1].n=GL[kk].WV[GL[kk].w_ncf+2].n=GL[kk].WV[GL[kk].w_ncf+3].n=wn;
	    GL[kk].WV[GL[kk].w_ncf+1].k=1; GL[kk].WV[GL[kk].w_ncf+2].k=2; GL[kk].WV[GL[kk].w_ncf+3].k=3;
            GL[kk].w_ncf+=3;
	    fprintf(stderr,"OK: %ld %ld %ld %ld %g %c %g\n",kk,wl,wm,wn,UL_RV[kk-1][wl-1][wm][wn],resolution_value,resolution_threshold);
	  }
	  else
	    fprintf(stderr,"BAD: %ld %ld %ld %ld %g %c %g\n",kk,wl,wm,wn,UL_RV[kk-1][wl-1][wm][wn],resolution_value,resolution_threshold);
        }
     fprintf(stderr,"Total %ld supports, %ld wavelet coefficients will be used in the inversion\n",
           GL[kk].w_N,GL[kk].w_ncf);
     fprintf(stderr,"Success.\n");
    }
    else
      NNz-=1;
   }

   fprintf(stderr,"Now rearranging SI wavelets basis...\n");

   w_s1_N=W_nl_1=0; /*W11: Nothing may be included*/
   for (wl=W_J;wl>=1;wl--)
   {
    if (write_si_pattern)
    {
      sprintf(buf,"si_pattern-l%ld.xy",wl);
      //fo_pattern=file_open(buf,"wt");
     fo_pattern=dir_file_open(output_path,buf,"wt");
    }
    for (wm=0;wm<=l2n(W_J-wl)-1;wm++)
      for (wn=0;wn<=l2n(W_J-wl)-1;wn++)
      {
     	if (resolution_test(resolution_value,resolution_threshold,SI_RV[wl-1][wm][wn])==1) {
	  W_nl_1++;
	  if ((wl==W_J) && (scalefunc_include)) /*The largest wavelet support also corresponds to the scaling function*/
	  {
	    w_s1_N++;
	    WV_1[w_s1_N].l=wl;
	    WV_1[w_s1_N].m=wm;
	    WV_1[w_s1_N].n=wn;
	    WV_1[w_s1_N].k=0; /*This stays for the scaling function*/
	  }
	    WV_1[w_s1_N+1].l=WV_1[w_s1_N+2].l=WV_1[w_s1_N+3].l=wl;
	    WV_1[w_s1_N+1].m=WV_1[w_s1_N+2].m=WV_1[w_s1_N+3].m=wm;
	    WV_1[w_s1_N+1].n=WV_1[w_s1_N+2].n=WV_1[w_s1_N+3].n=wn;
	    WV_1[w_s1_N+1].k=1; WV_1[w_s1_N+2].k=2; WV_1[w_s1_N+3].k=3;
            w_s1_N+=3;
	    
	    if (write_si_pattern) /*Save the pattern of the wavelet coverage in the sub-interface layer*/
	    {
	      fprintf(fo_pattern,"> -Wthick\n");
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimax-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimax-W_Si-1)*dd,y0+(wjmax-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+(wjmax-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"> -Wthin\n");
	      fprintf(fo_pattern,"%g %g\n",x0+((wimin-1+wimax)/2-W_Si-1)*dd,y0+(wjmin-1-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+((wimin-1+wimax)/2-W_Si-1)*dd,y0+(wjmax-W_Sj)*dd);
	      fprintf(fo_pattern,"> -Wthin\n");
	      fprintf(fo_pattern,"%g %g\n",x0+(wimin-1-W_Si-1)*dd,y0+((wjmin-1+wjmax)/2-W_Sj)*dd);
	      fprintf(fo_pattern,"%g %g\n",x0+(wimax-W_Si-1)*dd,y0+((wjmin-1+wjmax)/2-W_Sj)*dd);
	    }
	}
      }
    if (write_si_pattern) fclose(fo_pattern);
   }
   fprintf(stderr,"Total %ld supports, %ld wavelet coefficients will be used in the inversion\n",
          W_nl_1,w_s1_N);
   fprintf(stderr,"Success.\n");
   /*--------------------------------------------------------------*/
  
   wz_N=w_z1_N=0; /*Nothing may be included*/
   fprintf(stderr,"Now rearranging wavelets basis for the interface expansion...");
   for (wl=W_J;wl>=1;wl--)
    for (wm=0;wm<=l2n(W_J-wl)-1;wm++)
      for (wn=0;wn<=l2n(W_J-wl)-1;wn++)
      {
 	if (resolution_test(resolution_value,resolution_threshold,Z_RV[wl-1][wm][wn])==1) {
	  wz_N++;
	  if ((wl==W_J) && (scalefunc_include)) /*The largest wavelet support also corresponds to the scaling function*/
	  {
	    w_z1_N++;
	    WV_Z[w_z1_N].l=wl;
	    WV_Z[w_z1_N].m=wm;
	    WV_Z[w_z1_N].n=wn;
	    WV_Z[w_z1_N].k=0; 
	  }
	  WV_Z[w_z1_N+1].l=WV_Z[w_z1_N+2].l=WV_Z[w_z1_N+3].l=wl;
	  WV_Z[w_z1_N+1].m=WV_Z[w_z1_N+2].m=WV_Z[w_z1_N+3].m=wm;
	  WV_Z[w_z1_N+1].n=WV_Z[w_z1_N+2].n=WV_Z[w_z1_N+3].n=wn;
	  WV_Z[w_z1_N+1].k=1; WV_Z[w_z1_N+2].k=2; WV_Z[w_z1_N+3].k=3;
          w_z1_N+=3;
	}
      }
   fprintf(stderr,"Total %ld supports, %ld interface depth wavelet coefficients will be used\n",
           wz_N,w_z1_N);
   fprintf(stderr,"Success.\n");
   /*----------------------------------------------------------------------*/
  
   fprintf(stderr,"Freeing memory used for the resolution values...\n");
   for (kk=1;kk<=NNz;kk++)
     free_float_3d(UL_RV[kk-1],GL[kk].WI_J,l2n(GL[kk].WI_J-wl),l2n(GL[kk].WI_J-wl));
   delete[] UL_RV;
   free_float_3d(SI_RV,W_J,l2n(W_J-wl),l2n(W_J-wl));
   free_float_3d(Z_RV,W_J,l2n(W_J-wl),l2n(W_J-wl));
   
   delete[] REdiff;
   delete[] ND_R_norm;
   free_float_matrix(RM,Num_U,Num_U);
   
   fprintf(stderr,"Jumping back to matrix composition and inversion. Now with the new basis.\n");
   second_inversion=1;
   goto MATRIX_LABEL;
  }

  fprintf(stderr,"Freeing up memory...");
  free_bms(&CB);
  free_lsqr_mem(in_lsqr,out_lsqr,work_lsqr,func_lsqr);
  free_fvector(B,1,Num_E);
  fprintf(stderr,"Success.\n");

  fprintf(stderr,"Allocating memory for results...");
  alloc_float_3d(&VFINE,Nz,Ny,Nx);
  alloc_float_matrix(&Zintf,Ny,Nx);
  alloc_float_matrix(&UL_ft,Ny,Nx);
  fprintf(stderr,"Success.\n");

  fprintf(stderr,"Saving results...\n");

  //ferr=fopen("results_and_errors.dat","wt");
  ferr=dir_file_open(output_path,"results_and_errors.dat","wt");
  fprintf(ferr,"SI_slowness\n");
  for (wq=1;wq<=w_s1_N;wq++)
  {
    find_haar_indices(wq,WV_1,&wl,&wm,&wn,&wk);
    fprintf(ferr,"%5ld %5ld %5ld %5ld %g %g\n",wl,wm,wn,wk,X[wq],Err[wq]);
  }
  fprintf(ferr,"interface\n");
  for (wq=1;wq<=w_z1_N;wq++)
  {
    find_haar_indices(wq,WV_Z,&wl,&wm,&wn,&wk);
    fprintf(ferr,"%5ld %5ld %5ld %5ld %g %g\n",wl,wm,wn,wk,X[w_s1_N+wq],Err[w_s1_N+wq]);
  }
  N0=w_s1_N+w_z1_N;
  for (kk=1;kk<=NNz;kk++)
  {
    fprintf(ferr,"UPL slice %ld\n",kk);
    for (wq=1;wq<=GL[kk].w_ncf;wq++)
    {
      find_haar_indices(wq,GL[kk].WV,&wl,&wm,&wn,&wk);
      fprintf(ferr,"%5ld %5ld %5ld %5ld %g %g\n",wl,wm,wn,wk,X[N0+wq],Err[N0+wq]);
    }
    N0+=GL[kk].w_ncf;
  }
  fclose(ferr);

  fo=file_open(foS,"wb");
  //ferr=file_open("si_slowness_errors.2d","wb");
  ferr=dir_file_open(output_path,"si_slowness_errors.2d","wb");
  for (jj=1;jj<=NNy;jj++)
    for (ii=1;ii<=NNx;ii++)
    {
      for (wq=1,ess=ss=0.0;wq<=w_s1_N;wq++)
      {
        ss+=ev_haar(wq,WV_1,W_Si+ii,W_Sj+jj)*X[wq];
	ess+=pow((double)ev_haar(wq,WV_1,W_Si+ii,W_Sj+jj)*Err[wq],(double)2.0);
      }
      ess=sqrt(ess);
      fwrite(&ss,sizeof(float),1,fo);
      fwrite(&ess,sizeof(float),1,ferr);
      SI_dscm[jj][ii]=ss; /*From here now SI_dscm contains variation!!!!*/
    }
  fclose(ferr);
  fclose(fo);
  fo=file_open(foZ,"wb");
  //ferr=file_open("interface_errors.2d","wb");
  ferr=dir_file_open(output_path,"interface_errors.2d","wb");
  for (jj=1;jj<=NNy;jj++)
    for (ii=1;ii<=NNx;ii++)
    {
      for (wq=1,ess=ss=0.0;wq<=w_z1_N;wq++)
      {
        ss+=ev_haar(wq,WV_Z,W_Si+ii,W_Sj+jj)*X[w_s1_N+wq];
	ess+=pow((double)ev_haar(wq,WV_Z,W_Si+ii,W_Sj+jj)*Err[w_s1_N+wq],(double)2.0);
      }
      ess=sqrt(ess);
      fwrite(&ss,sizeof(float),1,fo);
      fwrite(&ess,sizeof(float),1,ferr);
      I_dzcm[jj][ii]=ss; /*From here now I_dzcm contains variation!!!!*/
    }
  fclose(ferr);
  fclose(fo);
  fprintf(stderr,"Success...\n");

  //write_2D_surface(I_dzcm,NNy,NNx,"I_dzcm.2d");
  
  N0=w_s1_N+w_z1_N;
  for (kk=1;kk<=NNz;kk++)
  {
    sprintf(buf,"ul%0ld_ds_err.2d",kk);
    //ferr=fopen(buf,"wb");
    ferr=dir_file_open(output_path,buf,"wb");
    sprintf(buf,"ul%0ld_ds.2d",kk);
    //fo=fopen(buf,"wb");
    fo=dir_file_open(output_path,buf,"wb");
    for (jj=1;jj<=NNy;jj++)
      for (ii=1;ii<=NNx;ii++)
      {
        wj=GL[kk].W_Sj+(jj-GL[kk].jj_min+1); wi=GL[kk].W_Si+(ii-GL[kk].ii_min+1);
        ess=ss=0.0;
	UL_dscm[kk][jj][ii]=DUMMY_VALUE;
	  for (wq=1;wq<=GL[kk].w_ncf;wq++)
	  {
            find_haar_supp(wq,GL[kk].WV,&wimin,&wimax,&wjmin,&wjmax);
	    if ((wi>=wimin) && (wi<=wimax) && (wj>=wjmin) && (wj<=wjmax))
	    {
	      ss+=ev_haar(wq,GL[kk].WV,wi,wj)*X[N0+wq];
	      ess+=pow((double)ev_haar(wq,GL[kk].WV,wi,wj)*Err[N0+wq],(double)2.0);
	    }
	  }
	  UL_dscm[kk][jj][ii]=ss;
	  ess=sqrt(ess);
	fwrite(&ess,sizeof(float),1,ferr);
        fwrite(&UL_dscm[kk][jj][ii],sizeof(float),1,fo);
      }
      fclose(ferr);
      fclose(fo);
      N0+=GL[kk].w_ncf;
  }
  /*From here now UL_dscm contains inverted slowness variation!!!!*/

  mod_dt=new float [Nl+Nr] - 1;
  
  fprintf(fstat,"Residual vector norms:\n             L2          Chi^2    \n");
  rL2=rChi2=0.0;
  if (Nl>0)
  {
   fprintf(stderr,"Computing model reflected time-residuals...\n");
   //fo=file_open("dtl_mod.dat","wt");
   fo=dir_file_open(output_path,"dtl_mod.dat","wt");
   y=z=0.0;
   for (m=1;m<=Nl;m++)
   {
    mod_dt[m]=0.0;
   }
   for (i=1;i<=N;i++)
   {
    pos_from_cln(&jj,&ii,i,NNx);
    for(it_p=WL_Cells[i-1].begin(); it_p!=WL_Cells[i-1].end(); it_p++) // for rays in cell
    {
        mod_dt[it_p->first]+=(it_p->second)*I_dzcm[jj][ii];
       //x+=WL[m][i]*I_dzcm[jj][ii];
    }
   }
    for (kk=1;kk<=NNz;kk++)
      for (jj=1;jj<=NNy;jj++)
        for (ii=1;ii<=NNx;ii++)
           for(it=Cells_reflected[kk-1][cln(jj,ii,NNx)-1].begin();
               it!=Cells_reflected[kk-1][cln(jj,ii,NNx)-1].end(); it++) // for rays in cell
	         mod_dt[it->first]+=UL_dscm[kk][jj][ii]*(it->second.l);
         // x+=UL_dscm[kk][jj][ii]*Cells_reflected[kk-1][cln(jj,ii,NNx)-1][m].l;
	 // x+=UL_dscm[kk][jj][ii]*URL[m][kk][cln(jj,ii,NNx)];

   for (m=1;m<=Nl;m++)
   {
    z+=(mod_dt[m]-dtl[m])*(mod_dt[m]-dtl[m]);
    y+=(mod_dt[m]-dtl[m])*(mod_dt[m]-dtl[m])/(stl[m]*stl[m]);
    fprintf(fo,"%10g %10g %10g\n",dtl[m],mod_dt[m],dtl[m]-mod_dt[m]);
   }

   rL2+=z;
   rChi2+=y;
   fprintf(fo,"residual L2 norm: %g; residual Chi^2 norm: %g\n",z,y/Nl);
   fprintf(fstat,"RW   : %13g %13g\n",z,y/Nl);
   fprintf(stderr,"Closing file...");
   fclose(fo);
   fprintf(stderr,"Closed\n");
   fprintf(stderr,"Success.\n");
  }
  
  fprintf(stderr,"Computing model headwaves time-residuals...\n");
  //fo=file_open("dth_mod.dat","wt");
  fo=dir_file_open(output_path,"dth_mod.dat","wt");
  y=z=0.0;
  for (m=1;m<=Nr;m++)//головные
  {
    mod_dt[m]=0.0;
  }
  for (i=1;i<=N;i++)
  {
    pos_from_cln(&jj,&ii,i,NNx);
    for(it=Cells_2D_bott[i-1].begin(); it!=Cells_2D_bott[i-1].end(); it++) // for rays in cell
    {
	  mod_dt[it->first]+=(it->second).l*SI_dscm[jj][ii];
    }

     for(it_p=W_Cells[i-1].begin(); it_p!=W_Cells[i-1].end(); it_p++) // for rays in cell
    {
	  mod_dt[it_p->first]+=(it_p->second)*I_dzcm[jj][ii];
    }

   // x+=Cells_2D_bott[i-1][m].l*SI_dscm[jj][ii]+W_Cells[i-1][m]*I_dzcm[jj][ii];
      //x+=L[m][i]*SI_dscm[jj][ii]+W[m][i]*I_dzcm[jj][ii];
  }
    for (kk=1;kk<=NNz;kk++)
      for (jj=1;jj<=NNy;jj++)
        for (ii=1;ii<=NNx;ii++)
          for(it=Cells_hw_top[kk-1][cln(jj,ii,NNx)-1].begin();
              it!=Cells_hw_top[kk-1][cln(jj,ii,NNx)-1].end(); it++) // for rays in cell
           {
	        mod_dt[it->first]+=UL_dscm[kk][jj][ii]*((it->second).l);
           }
         //x+=UL_dscm[kk][jj][ii]* Cells_hw_top[kk][cln(jj,ii,NNx)][m].l;
	    // x+=UL_dscm[kk][jj][ii]*UHL[m][kk][cln(jj,ii,NNx)];
  for (m=1;m<=Nr;m++)//головные
  {
    z+=(mod_dt[m]-dt[m])*(mod_dt[m]-dt[m]);
    y+=(mod_dt[m]-dt[m])*(mod_dt[m]-dt[m])/(st[m]*st[m]);
    fprintf(fo,"%10g %10g %10g\n",dt[m],mod_dt[m],dt[m]-mod_dt[m]);
  }

  rL2+=z;
  rChi2+=y;
  fprintf(fo,"residual L2 norm: %g; residual Chi^2 norm: %g\n",z,y/Nr);
  fprintf(fstat,"HW   : %13g %13g\n",z,y/Nr);
  fprintf(fstat,"R+H  : %13g %13g\n",rL2,rChi2/(Nl+Nr));
  fclose(fo);
  fprintf(stderr,"Success.\n");
  
  /*@@@@---3D-w14 BLOCK CHANGED START---@@@@*/
  /*------------------------Model back-interpolation on to fine grid-------------------*/
  /*------Read reference true velocity model-------------------------------------------*/
  fprintf(stderr,"Reading fine velocity model...");
  read_3D_volume(VFINE,Nz,Ny,Nx,fi_VelMod);
  fprintf(stderr,"Success.\n");
  /*------Read fine grids of layer top and interface*/
  fprintf(stderr,"Read surfaces...");
  read_2D_surface(UL_ft,Ny,Nx,fi_fiTop);
  read_2D_surface(Zintf,Ny,Nx,fi_fiBas);
  fprintf(stderr,"Success.\n");

  /*------Second step: interpolation of the slowness variations on to fine grid--------*/
  /*------kk,jj,ii - coordinates of the coarse grid cells vortices (numbered from 0 to NNx-1)*/
  /*------k,i,j - coordinates of the fine grid nodes (numbered from 0 to Nx)*/
  fprintf(stderr,"Interpolation of slowness variations on the fine grid...\n");

  fprintf(stderr,"First interpolate sub-interface slowness variation...");
  alloc_float_matrix(&SI_dsv_fine,Ny,Nx);
  /*3D-W14: new interpolation procedure with nodes centered in the centers of coarse grid cells*/
  for (j=0;j<Ny;j++)
  {
    jj=floor((h*(float)j)/dd-0.5)+1;
    if (jj<1) jj=1; if (jj>=NNy) jj=NNy-1;
    y=(float)j*h-((float)jj-0.5)*dd;
    for (i=0;i<Nx;i++)
    {
      ii=floor((h*(float)i)/dd-0.5)+1;
      if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;
      x=(float)i*h-((float)ii-0.5)*dd;
      SI_dsv_fine[j][i]=(1.0/(dd*dd))*
	          (SI_dscm[jj][ii]*(dd-x)*(dd-y)+
	           SI_dscm[jj][ii+1]*x*(dd-y)+
	           SI_dscm[jj+1][ii]*(dd-x)*y+
	           SI_dscm[jj+1][ii+1]*x*y);
    }
  }

  if ((Rsmooth_ds>0) && (Psmooth_ds>0))
  {
    fprintf(stderr,"\nAdditional smoothing of the sub-interface slowness variation...");
    smooth_interface(SI_dsv_fine,Ny,Nx,x0,y0,h,Rsmooth_ds,Psmooth_ds);
    fprintf(stderr,"Success.\n");
  }

  fprintf(stderr,"Success.\n");

  fprintf(stderr,"Now interpolate UPL slowness variations...");
  if (NNz!=1)
  {
   for (k=0;k<Nz;k++)
   {
    /*kk=(k-(uk_min-1))/k_gen_v; if (kk==NNz) kk--;*/
    kk=floor(h*(float)(k-uk_min+1-0.5)/ddz)+1;
    if (kk<1) kk=1; if (kk>=NNz) kk=NNz-1;
    z_node=z0+h*(float)k;
    z=h*(float)(k-uk_min+1)-((float)kk-0.5)*ddz;
    for (j=0;j<Ny;j++)
    {
      jj=floor((h*(float)j)/dd-0.5)+1;
      if (jj<1) jj=1; if (jj>=NNy) jj=NNy-1;
      y=(float)j*h-((float)jj-0.5)*dd;
      for (i=0;i<Nx;i++)
      {
        ii=floor((h*(float)i)/dd-0.5)+1;
        if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;
        x=(float)i*h-((float)ii-0.5)*dd;
	    if ((k>floor((UL_ft[j][i]-z0)/h)) && (k<floor((Zintf[j][i]-z0)/h))) /*Node inside reference UPL*/
	    {
	      
	    ss=(1.0/(dd*dd*ddz))*
	       (UL_dscm[kk][jj][ii]*(dd-x)*(dd-y)*(ddz-z)+
	        UL_dscm[kk][jj][ii+1]*x*(dd-y)*(ddz-z)+
	        UL_dscm[kk][jj+1][ii]*(dd-x)*y*(ddz-z)+
	        UL_dscm[kk+1][jj][ii]*(dd-x)*(dd-y)*z+
	        UL_dscm[kk+1][jj][ii+1]*x*(dd-y)*z+
	        UL_dscm[kk+1][jj+1][ii]*(dd-x)*y*z+
	        UL_dscm[kk][jj+1][ii+1]*x*y*(ddz-z)+
	        UL_dscm[kk+1][jj+1][ii+1]*x*y*z);
	      
	      /*ss=UL_dscm[kk][jj][ii];*/
	    }
	    else /*Node outside reference UPL*/
	    {
	      if ((k<=floor((UL_ft[j][i]-z0)/h))) ss=0.0;
	      else /*Node below the interface*/
	        ss=SI_dsv_fine[j][i];
	    }
	    dSFINE[k][j][i]=ss;
	    VFINE[k][j][i]=VFINE[k][j][i]/(1.0+VFINE[k][j][i]*ss);
      }
    }
   }
  }
  else //2D model in plane XY
  {
    fprintf(stderr,"2D case\n");
    for (j=0;j<Ny;j++)
    {
      jj=floor((h*(float)j)/dd-0.5)+1;
      if (jj<1) jj=1; if (jj>=NNy) jj=NNy-1;
      y=(float)j*h-((float)jj-0.5)*dd;
      for (i=0;i<Nx;i++)
      {
        ii=floor((h*(float)i)/dd-0.5)+1;
        if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;
        x=(float)i*h-((float)ii-0.5)*dd;

	    ss=(1.0/(dd*dd))*
	       (UL_dscm[NNz][jj][ii]*(dd-x)*(dd-y)+
	        UL_dscm[NNz][jj][ii+1]*x*(dd-y)+
	        UL_dscm[NNz][jj+1][ii]*(dd-x)*y+
	        UL_dscm[NNz][jj+1][ii+1]*x*y);
	
	//fprintf(stderr,"%d %d %g\n",jj,ii,ss);
	//fprintf(stderr,"%d\n",Nz);
 
	 for (k=0;k<Nz;k++)
         {
	    dSFINE[k][j][i]=ss;
	    VFINE[k][j][i]=VFINE[k][j][i]/(1.0+VFINE[k][j][i]*ss);
	    //fprintf(stderr,"%d %g %g\n",k,ss,VFINE[k][j][i]);
	 }
      }
    }
  }
  fprintf(stderr,"Success.\n");

  free_float_matrix(SI_dsv_fine,Ny,Nx);

  /*------Second step: correcting the model due to the shift of the interface----------*/
  fprintf(stderr,"Computing the new position of the interface on the fine grid...");
  /*Compute the new interface first*/
  alloc_float_matrix(&Zintf_new,Ny,Nx);
  for (j=0;j<Ny;j++)
  {
    jj=floor((h*(float)j)/dd-0.5)+1;
    if (jj<1) jj=1; if (jj>=NNy) jj=NNy-1;
    y=(float)j*h-((float)jj-0.5)*dd;
    for (i=0;i<Nx;i++)
    {
      ii=floor((h*(float)i)/dd-0.5)+1;
      if (ii<1) ii=1; if (ii>=NNx) ii=NNx-1;
      x=(float)i*h-((float)ii-0.5)*dd;
      ss=(1.0/(dd*dd))*
	  (I_dzcm[jj][ii]*(dd-x)*(dd-y)+
	   I_dzcm[jj][ii+1]*x*(dd-y)+
	   I_dzcm[jj+1][ii]*(dd-x)*y+
	   I_dzcm[jj+1][ii+1]*x*y);

      Zintf_new[j][i]=Zintf[j][i]+ss;
    }
  }
  fprintf(stderr,"Success.\n");
  /*@@@@---3D-w14 BLOCK CHANGED END---@@@@*/

  if ((Rsmooth>0) && (Psmooth>0))
  {
    fprintf(stderr,"Additional smoothing of the interface...");
    smooth_interface(Zintf_new,Ny,Nx,x0,y0,h,Rsmooth,Psmooth);
    fprintf(stderr,"Success.\n");
  }

  fprintf(stderr,"Correcting model due to the shift of the interface...");
  for (j=0;j<Ny;j++)
  {
    for (i=0;i<Nx;i++)
    {
      ss=Zintf_new[j][i]-Zintf[j][i];
      if (ss<0) /*Shift upward*/
      {
        k=(Zintf[j][i]-z0)/h; /*Old interface position k index (last node over the interface)*/
	kk=(Zintf[j][i]+ss-z0)/h; /*New interface position k index (last node over the interface)*/
	if (kk<0)
	{
	  fprintf(stderr,"WARNNIG: Interface position over the model top. Was set to top.");
	  kk=0;
	}
	for (kkk=k;kkk>kk;kkk--) VFINE[kkk][j][i]=VFINE[k+1][j][i];
      }
      else if (ss>0) /*Shift dounward*/
      {
        k=(Zintf[j][i]-z0)/h; /*Old interface position k index (last node over the interface)*/
	kk=(Zintf[j][i]+ss-z0)/h; /*New interface position k index (last node over the interface)*/
	if (kk>=Nz)
	{
	  fprintf(stderr,"WARNNIG: Interface position below the model bottom. Was set to bottom.");
	  kk=Nz-1;
	}
	for (kkk=k;kkk<=kk;kkk++) VFINE[kkk][j][i]=VFINE[k-1][j][i];
      }
    }
  }
  fprintf(stderr,"Success.\n");

  /*------Saving results---------------------------------------------------------------*/
  fprintf(stderr,"Saving fine grid models...");
  write_3D_volume(VFINE,Nz,Ny,Nx,output_path,"modified_velmod.3d");
  write_3D_volume(dSFINE,Nz,Ny,Nx,output_path,"ds.3d");
  write_2D_surface(Zintf_new,Ny,Nx,output_path,"modified_interface.2d");
  fprintf(stderr,"Success.\n");
  /*-----------------------------------------------------------------------------------*/

  /*@@@@---3D-w13 BLOCK CHANGED START---@@@@*/
  fprintf(stderr,"Freeing up memory...");
  free_float_matrix(Zintf_new,Ny,Nx);
  free_f3tensor(UL_dscm,1,NNz,1,NNy,1,NNx);
  free_float_3d(dSFINE,Nz,Ny,Nx);
  free_float_3d(VFINE,Nz,Ny,Nx);
  free_float_matrix(Zintf,Ny,Nx);
  free_float_matrix(UL_ft,Ny,Nx);
  free_matrix(SI_dscm,1,NNy,1,NNx);
  free_matrix(I_dzcm,1,NNy,1,NNx);
  fprintf(stderr,"Success.\n");


  delete[]  Cells_hw_top[0];
  delete[] Cells_hw_top;

  delete[] Cells_2D_bott;

  delete[] Cells_reflected[0];
  delete[] Cells_reflected; 

  delete[] W_Cells;
   
  if (Nl!=0) delete[] WL_Cells;


  /*@@@@---3D-w13 BLOCK CHANGED END---@@@@*/

  /*
  fprintf(stderr,"Computing model velocity-density relationship...\n");
  fo=file_open("corr_mod.dat","wt");
  y=0.0;
  for (i=1;i<=N;i++)
  {
    x=-alpha*a*X[i]+alpha*X[2*N+i];
    fprintf(fo,"%10g %10g %10g\n",x,B[Nl+Nr+Ng+i],X[2*N+i]/X[i]);
    y+=x*x;
  }
  fprintf(fo,"residual norm: %g\n",y);
  fclose(fo);
  fprintf(stderr,"Sccess.\n");
  */

  /*--------------------------------------------------------------*/
  fprintf(stderr,"Freeing up memory...\n");

  for (kk=1;kk<=NNz;kk++)
  {
    fmemfree((char*)(GL[kk].WV+1));
    fmemfree((char*)(GL[kk].CV+1));
  }
  fmemfree((char*)(GL+1));
  //for (i=1;i<=Nl;i++) free_matrix(URL[i],1,NNz,1,N);
  //fmemfree((char*)(URL+1));
  //for (i=1;i<=Nr;i++) free_matrix(UHL[i],1,NNz,1,N);
  //fmemfree((char*)(UHL+1));
  fmemfree((char*)(WV_Z+1));
  fmemfree((char*)(WV_1+1));
  free_fvector(X,1,Num_U);
  free_fvector(ws,1,Nw);
  free_fvector(wz,1,Nw);
  free_fvector(wx,1,Nw);
  free_fvector(wy,1,Nw);
  //if (Nl!=0) free_matrix(WL,1,Nr,1,N);
  free_fvector(stl,1,Nr);
  free_fvector(dtl,1,Nr);
  //free_matrix(W,1,Nr,1,N);
  free_fvector(st,1,Nr);
  free_fvector(dt,1,Nr);
  //free_matrix(L,1,Nr,1,N);
  free_fvector(T,1,N);
  free_fvector(Z,1,N);
  fprintf(stderr,"Success.\n");
  fprintf(stderr,"Program finished successfully. Good lack!\n");

  time(&tmE);
  ctime_r(&tmE,time_string);
  fprintf(fstat,"Program stops: "); fputs(time_string,fstat);
  i=(int)(tmE-tm0);
  j=i/3600;
  k=(i-j*3600)/60;
  fprintf(fstat,"Runtime: %ld-%ld-%ld\n",j,k,i-j*3600-k*60);
  fclose(fstat);

  return 0;
}



