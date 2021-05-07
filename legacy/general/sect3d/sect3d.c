/*Making sclices from the 3d file*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "submem.h"
#include "subio.h"
#include "gmt_nan.h"

#define GEOMFILE "model-geometry.gras"
#define STRLN 256

int gras_round(float x)
{
  
  return (x-floor(x)) >= 0.5 ? ceil(x) : floor(x);
}

int main(int narg, char **argv)
{
  long Nx, Ny, Nz;
  float x0, y0 , z0;
  int i, i_g, j, j_g, k, k_g;
  float x, x_g, y_g;
  float Xa, Ya, Xb, Yb;
  float Lprof,cos_a,sin_a;
  float dh;
  float *V, *S;
  char buf[STRLN];
  char fln_in[STRLN], fln_o[STRLN];
  
  struct
  {
    int nx;
    int ny;
    int node_offset;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
    double x_inc;
    double y_inc;
    double z_scale_factor;
    double z_add_offset;
    char x_units[80];
    char y_units[80];
    char z_units[80];
    char title[80];
    char command[512];
    char remark[1000];
  } ncg_h;
  
  FILE *fp, *fi, *fo;
  
  float GMT_f_NaN;

  GMT_make_fnan(GMT_f_NaN);
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&dh);
  fclose(fp);
  
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
  fgets(buf,STRLN,fp);
  if (fp!=stdin) fclose(fp);

  V=(float*)fmemalloc(Nx*Ny*Nz,sizeof(float));
  fi=file_open(fln_in,"rb");
  for (k=0;k<Nz;k++)
    for (i=0;i<Ny;i++)
      fread(V+Nx*Ny*k+Nx*i,sizeof(float),Nx,fi);
  fclose(fi);
  
  ncg_h.nx=floor(pow((Xa-Xb)*(Xa-Xb)+(Ya-Yb)*(Ya-Yb),0.5)/dh)+1;
  ncg_h.x_min=0;
  ncg_h.x_max=(ncg_h.nx-1)*dh;
  ncg_h.x_inc=dh;
  
  ncg_h.ny=Nz;
  ncg_h.y_min=-z0-0.5*dh-(Nz-1)*dh;
  ncg_h.y_max=-z0-0.5*dh;
  ncg_h.y_inc=dh;
  
  ncg_h.node_offset=0;
  ncg_h.z_scale_factor=1;
  ncg_h.z_add_offset=0;
  
  strcpy(ncg_h.x_units,"");
  strcpy(ncg_h.y_units,"");
  strcpy(ncg_h.z_units,"");
  strcpy(ncg_h.title,"Section");
  sprintf(ncg_h.command,"%s %s %s",argv[0],argv[1],argv[2]);
  sprintf(ncg_h.remark,"Section from file %s, Xa=%g Ya=%g Xb=%g Yb=%g",fln_in,Xa,Ya,Xb,Yb);
  
  Lprof=pow((Xa-Xb)*(Xa-Xb)+(Ya-Yb)*(Ya-Yb),0.5);
  cos_a=(Xb-Xa)/Lprof;
  sin_a=(Yb-Ya)/Lprof;
  
  S=(float*)fmemalloc(ncg_h.nx*ncg_h.ny,sizeof(float));
  
  ncg_h.z_min=1.7e38;
  ncg_h.z_max=-1.7e38;
  for (i=0;i<ncg_h.nx;i++)
  {
      x=ncg_h.x_min+i*dh;
      x_g=Xa+x*cos_a;
      y_g=Ya+x*sin_a;
      i_g=(int)((x_g-x0)/dh);
      j_g=(int)((y_g-y0)/dh);
      if (i_g<0 || i_g>=Nx || j_g<0 || j_g>=Ny)
      {
	fprintf(stderr,"x: %g cos_a: %g sin_a: %g %g %g %d %d\n",x,cos_a,sin_a,x_g,y_g,i_g,j_g);
	exit(EXIT_FAILURE);
      }
      for (j=0;j<ncg_h.ny;j++)
      {
	k_g=j;
	*(S+j*ncg_h.nx+i)=*(V+Nx*Ny*k_g+Nx*j_g+i_g);
	if (*(S+j*ncg_h.nx+i) < 1e9)
	{
	  if (ncg_h.z_min > (*(S+j*ncg_h.nx+i))) ncg_h.z_min=*(S+j*ncg_h.nx+i);
	  if (ncg_h.z_max < (*(S+j*ncg_h.nx+i))) ncg_h.z_max=*(S+j*ncg_h.nx+i);
	}
	else
	  *(S+j*ncg_h.nx+i)=GMT_f_NaN;
      }
  }
   
  fo=fopen(fln_o,"wb");
  
  if (fwrite ((void*)&(ncg_h.nx),3*sizeof(int),(size_t)1,fo) != 1 ||
      fwrite ((void*)&(ncg_h.x_min),sizeof(ncg_h)-((long)&(ncg_h.x_min)-(long)&(ncg_h.nx)),(size_t)1,fo) != 1)
  {
     fprintf(stderr,"Error writing grid header of gmt file");
     exit(EXIT_FAILURE);
  }
  for (j=0;j<ncg_h.ny;j++)
    fwrite(S+j*ncg_h.nx,ncg_h.nx,sizeof(float),fo);
  fclose(fo);
  
  fmemfree(S);
  fmemfree(V);
  return 0;
}
