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
  float *Vt, *Vb;
  float ***Vmod, **Intf,***top,***bot;
  float xc, yc, V0, gradV;
  float x0, y0, z0, h;
  float z, v;
  char outfile[FLNLN];
  char **botfiles;
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
  botfiles=fmemalloc(Nlay+1,sizeof(char*));
    for (i=0;i<=Nlay;i++)
      botfiles[i]=fmemalloc(FLNLN,sizeof(char));
  Vt=fmemalloc(Nlay,sizeof(float));
  Vb=fmemalloc(Nlay,sizeof(float));
  
  if (ask) fprintf(stdout,"1 layer: input <top depth file> > ");
  fscanf(fp,"%s",botfiles[0]);
  fprintf(stdout,"%s\n",botfiles[0]);
  for (i=0;i<Nlay;i++)
  {
    if (ask) fprintf(stdout,"%d layer: input <bottom depth file> V_t V_b > ",i);
    fscanf(fp,"%s%g%g",botfiles[i+1],Vt+i,Vb+i);
    fprintf(stdout,"%s %g %g\n",botfiles[i+1],Vt[i],Vb[i]);
  }
    
  alloc_float_3d(&Vmod,Nz,Ny,Nx);

  top=fmemalloc(Nlay,sizeof(float**));
  bot=fmemalloc(Nlay,sizeof(float**));
  for (i=0;i<Nlay;i++)
  {
    alloc_float_matrix(&(top[i]),Ny,Nx);
    alloc_float_matrix(&(bot[i]),Ny,Nx);
    read_2D_surface(top[i],Ny,Nx,botfiles[i]);
    read_2D_surface(bot[i],Ny,Nx,botfiles[i+1]);
  }
  
  for (k=0;k<Nz;k++)
  {
    z=z0+k*h;
    for (i=0;i<Ny;i++)
    {
      for (j=0;j<Nx;j++)
      {        
        if (z<top[0][i][j]) kk=-1; else      
        for (kk=0;kk<Nlay;kk++) if ((z>=top[kk][i][j]) && (z<bot[kk][i][j])) break;        
        if (kk==-1) fprintf(stderr,"Nodes with z=%g (k=%d) are above the model top.\n",z,k);        
        if (kk==Nlay) fprintf(stderr,"Nodes with z=%g (k=%d) are below the model bottom.\n",z,k);        
        if ((kk>=0) && (kk<Nlay))        
          v=Vt[kk]+(Vb[kk]-Vt[kk])/(bot[kk][i][j]-top[kk][i][j])*(z-top[kk][i][j]);        
        else        
          v=DUMMY_VALUE;        
    	Vmod[k][i][j]=v;        
      }    
    }
  }
  write_3D_volume(Vmod,Nz,Ny,Nx,"vel_ref.3d");  
  
  for (i=0;i<Nlay;i++)
  {
    free_float_matrix(top[i],Ny,Nx);
    free_float_matrix(bot[i],Ny,Nx);
    fmemfree(botfiles[i]);
  }
  fmemfree(top);
  fmemfree(bot);
  fmemfree(botfiles);
  
  free_float_3d(Vmod,Nz,Ny,Nx);
  fmemfree(Vb);
  fmemfree(Vt);
}
