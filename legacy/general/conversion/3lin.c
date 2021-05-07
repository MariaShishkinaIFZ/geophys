#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"

#define DUMMY_LIM 1e9
#define DUMMY_VALUE 1e10

/*Coordinate system: z-layers (1st index), y-rows (2nd index), x-columns (3rd index)*/
int lin3D
(float ***src, float ***dest, long Nz, long Ny, long Nx, float *zV, float *yV, float *xV,
 float z0, float y0, float x0, float dz, float dy, float dx, long nnx, long nny, long nnz)
{
  long k,j,i,kk,ii,jj,n;
  float z_nd,y_nd,x_nd,z,y,x,ddx,ddy,ddz;
  
  n=0;
  for (k=0;k<nnz;k++)
  {
    z_nd=z0+dz*k;
    if (z_nd<zV[0]) kk=-1; else
    for (kk=0;kk<Nz-1;kk++) if ((z_nd>=zV[kk]) && (z_nd<=zV[kk+1])) break;
    for (j=0;j<nny;j++)
    {
      y_nd=y0+dy*j;
      if (y_nd<yV[0]) jj=-1; else
      for (jj=0;jj<Ny-1;jj++) if ((y_nd>=yV[jj]) && (y_nd<=yV[jj+1])) break;
      for (i=0;i<nnx;i++)
      {
        x_nd=x0+dx*i;
	if (x_nd<xV[0]) ii=-1; else
	for (ii=0;ii<Nx-1;ii++) if ((x_nd>=xV[ii]) && (x_nd<=xV[ii+1])) break;
	if ((kk>=0) && (kk<Nz-1) && (jj>=0) && (jj<Ny-1) && (ii>=0) && (ii<Nx-1))
	{
	  z=z_nd-zV[kk]; ddz=zV[kk+1]-zV[kk];
	  y=y_nd-yV[jj]; ddy=yV[jj+1]-yV[jj];
	  x=x_nd-xV[ii]; ddx=xV[ii+1]-xV[ii];
	  if ((src[kk][jj][ii]<DUMMY_LIM) && 
	      (src[kk+1][jj][ii]<DUMMY_LIM) &&
	      (src[kk][jj+1][ii]<DUMMY_LIM) &&
	      (src[kk][jj][ii+1]<DUMMY_LIM) &&
	      (src[kk+1][jj+1][ii]<DUMMY_LIM) &&
	      (src[kk+1][jj][ii+1]<DUMMY_LIM) &&
	      (src[kk][jj+1][ii+1]<DUMMY_LIM) &&
	      (src[kk+1][jj+1][ii+1]<DUMMY_LIM))
	  {
	      
	     dest[k][j][i]=(1.0/(ddz*ddy*ddx))*
	       (src[kk][jj][ii]*(ddx-x)*(ddy-y)*(ddz-z)+
	        src[kk][jj][ii+1]*x*(ddy-y)*(ddz-z)+
	        src[kk][jj+1][ii]*(ddx-x)*y*(ddz-z)+
	        src[kk+1][jj][ii]*(ddx-x)*(ddy-y)*z+
	        src[kk+1][jj][ii+1]*x*(ddy-y)*z+
	        src[kk+1][jj+1][ii]*(ddx-x)*y*z+
	        src[kk][jj+1][ii+1]*x*y*(ddz-z)+
	        src[kk+1][jj+1][ii+1]*x*y*z);
	     n++;
	  }
	  else
	    dest[k][j][i]=DUMMY_VALUE;
	}
      }
    }
  }
  return n;
}

int main(int narg, char **argv)
{
  long Nx, Ny, Nz;
  long nnx, nny, nnz;
  long i,j,k;
  float x0,y0,z0;
  float h;
  float ***data;
  float ***res;
  float *x, *y, *z;
  FILE *fi, *fo, *fp;
  
  if (narg<3)
  {
    fprintf(stderr,"Insufficient parameters. Use:\n3lin <data file> <out file [.3d]>\n");
    exit(EXIT_FAILURE);
  }
  
  fp = file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&nnx,&nny,&nnz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
  
  fi = file_open(argv[1],"rt");
  fscanf(fi,"%ld%ld%ld",&Nx,&Ny,&Nz);
  x=(float*)fmemalloc(Nx,sizeof(float));
  y=(float*)fmemalloc(Ny,sizeof(float));
  z=(float*)fmemalloc(Nz,sizeof(float));
  alloc_float_3d(&data,Nz,Ny,Nx);
  alloc_float_3d(&res,nnz,nny,nnx);
  for (i=0;i<Nx;i++) fscanf(fi,"%g",x+i);
  for (i=0;i<Ny;i++) fscanf(fi,"%g",y+i);
  for (i=0;i<Nz;i++) fscanf(fi,"%g",z+i);
  
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      for (i=0;i<Nx;i++)
        fscanf(fi,"%g",data[k][j]+i);
  
  fclose(fi);
  
  fprintf(stderr,"Data read. Interpolation...");
  lin3D(data,res,Nz,Ny,Nx,z,y,x,z0,y0,x0,h,h,h,nnx,nny,nnz);
  fprintf(stderr,"Success.\n");
  
  write_3D_volume(res,nnz,nny,nnx,argv[2]);
  
  free_float_3d(res,nnz,nny,nnx);
  free_float_3d(data,Nz,Ny,Nx);
  fmemfree(x);
  fmemfree(y);
  fmemfree(z);
}
