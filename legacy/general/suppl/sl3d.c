/*Making sclices from the 3d file*/

#include <stdio.h>
#include <math.h>

#include "submem.h"
#include "subio.h"

#define GEOMFILE "model-geometry.gras"
#define STRLN 256

int gras_round(float x)
{
  
  return (x-floor(x)) >= 0.5 ? ceil(x) : floor(x);
}

int main(int narg, char **argv)
{
  long Nx, Ny, Nz;
  int i, j, k;
  unsigned pos; 
  float x0, y0, z0;
  float xs, ys, zs;
  unsigned is, js, ks;
  float dh;
  char cs;
  float *V;
  char buf[STRLN];
  char fln_in[STRLN], fln_o[STRLN];
  FILE *fp, *fi, *fo;
  
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
  fscanf(fp,"%c%g",&cs,&xs);
  switch (cs)
  {
    case 'X': if ((xs<x0) || (xs>x0+(Nx-1)*dh)) errquit("Slice position is outside of volume.\n"); break;
    case 'Y': ys=xs; if ((ys<y0) || (ys>y0+(Ny-1)*dh)) errquit("Slice position is outside of volume.\n"); 
              break;
    case 'Z': zs=xs; if ((zs<z0) || (zs>z0+(Nz-1)*dh)) errquit("Slice position is outside of volume.\n");
              break;
    default: errquit("Incorrect label for slice position. MUST be either X, Y or Z.\n");
  }
  fgets(buf,STRLN,fp);
  if (fp!=stdin) fclose(fp);
  
  V=(float*)fmemalloc(Nx*Ny*Nz,sizeof(float));
  fi=file_open(fln_in,"rb");
  for (k=0;k<Nz;k++)
    for (i=0;i<Ny;i++)
      fread(V+Nx*Ny*k+Nx*i,sizeof(float),Nx,fi);
  fclose(fi);
  
  fo=fopen(fln_o,"wb");
  switch (cs)
  {
    case 'Z':
      ks = gras_round((zs-z0)/dh);
      printf("Zk = %d\n",ks);
      for (i=0;i<Ny;i++)
        fwrite(V+Nx*Ny*ks+Nx*i,sizeof(float),Nx,fo);
      break;
    case 'Y':
      is = gras_round((ys-y0)/dh);
      printf("Yi = %d\n",is);
      for (k=Nz-1;k>=0;k--)
        fwrite(V+Nx*Ny*k+Nx*is,sizeof(float),Nx,fo);
      break;
    case 'X':
      js = gras_round((xs-x0)/dh);
      printf("Xj = %d\n",js);
      for (k=Nz-1;k>=0;k--)
        for (i=0;i<Ny;i++)
	  fwrite(V+Nx*Ny*k+Nx*i+js,sizeof(float),1,fo);
  }
  fclose(fo);
  
  fmemfree(V);
  return 0;
}
