#include <stdio.h>
#include <math.h>
#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 80 	/*Maximum length of the file name*/

#define DUMMY_VALUE 1e10

int main (int narg, char **argv)
{
  int Nx, Ny, Nz, i, j, k, kk;
  int Nlay;
  int ask;
  float *Zb, *Vt, *Vb;
  float ***Vmod, **Intf;
  float xc, yc, V0, gradV;
  float x0, y0, z0, h;
  float z, v;
  char outfile[FLNLN];
  FILE *fo, *fp;
  
  fp=fopen(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
  
  if (narg < 2) errquit("Insufficient parameters.\n");
  if (strcmp(argv[1],"-f")==0) /*Use parameters file with the name given as 2nd parameter*/
  {
    if (narg < 3) errquit("Parameters file name not given.\n");
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
  
  if (ask) fprintf(stdout,"Nlay > ");
  fscanf(fp,"%d",&Nlay);
  Zb=fmemalloc(Nlay+1,sizeof(float));
  Vt=fmemalloc(Nlay,sizeof(float));
  Vb=fmemalloc(Nlay,sizeof(float));
  Zb[0]=z0;
  for (i=0;i<Nlay;i++)
  {
    if (ask) fprintf(stdout,"%d layer: input Z_b V_t V_b > ");
    fscanf(fp,"%g%g%g",Zb+i+1,Vt+i,Vb+i);
  }
    
  alloc_float_3d(&Vmod,Nz,Ny,Nx);
  alloc_float_matrix(&Intf,Ny,Nx);
  
  for (k=0;k<Nz;k++)
  {
    z=z0+k*h;
    if (z<Zb[0]) kk=-1; else
    for (kk=0;kk<Nlay;kk++) if ((z>=Zb[kk]) && (z<Zb[kk+1])) break;
    if (kk==-1) fprintf(stderr,"Nodes with z=%g (k=%d) are above the model top.\n",z,k);
    if (kk==Nlay) fprintf(stderr,"Nodes with z=%g (k=%d) are below the model bottom.\n",z,k);
    if ((kk>=0) && (kk<Nlay))
      v=Vt[kk]+(Vb[kk]-Vt[kk])/(Zb[kk+1]-Zb[kk])*(z-Zb[kk]);
    else
      v=DUMMY_VALUE;
    for (j=0;j<Ny;j++)
    {
      for (i=0;i<Nx;i++)
      {
	Vmod[k][j][i]=v;
      }
    }
  }
  write_3D_volume(Vmod,Nz,Ny,Nx,"vel_ref.3d");  
  
  for (kk=1;kk<=Nlay;kk++)
  {
    sprintf(outfile,"intf%0d_ref.2d",kk);
    for (j=0;j<Ny;j++)
    {
      for (i=0;i<Nx;i++)
      {
	Intf[j][i]=Zb[kk];
      }
    }
    write_2D_surface(Intf,Ny,Nx,outfile);
  }
  
  free_float_matrix(Intf,Ny,Nx);
  free_float_3d(Vmod,Nz,Ny,Nx);
  fmemfree(Vb);
  fmemfree(Vt);
  fmemfree(Zb);
}
