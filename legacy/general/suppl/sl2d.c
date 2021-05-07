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


void read_2D_surface(float **Buf, long Ny, long Nx, const char *filename)
{
  long j;
  FILE *fi;

  fi=file_open(filename,"rb");
  for (j=0;j<Ny;j++)
    fread(Buf[j],sizeof(float),Nx,fi);
  fclose(fi);
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
  float **V;
  char buf[STRLN];
  char fln_in[STRLN], fln_o[STRLN];
  FILE *fp, *fi, *fo;
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&dh);
  fclose(fp);
  
  if (narg < 2) errquit("Insufficient parameters\nUSAGE: sl2d -f parfilename | -s\n");
  if (strcmp(argv[1],"-f")==0)
  {
    if (narg < 3) errquit("Insufficient parameters\nUSAGE: sl2d -f parfilename | -s\n");
    fp=file_open(argv[2],"rt");
  }
  else
    if (strcmp(argv[1],"-s")!=0) errquit("Unknown option\nUSAGE: sl2d -f parfilename | -s\n");
    else fp=stdin;
    
  fscanf(fp,"%s",fln_in); fgets(buf,STRLN,fp);
  fscanf(fp,"%s",fln_o); fgets(buf,STRLN,fp);
  fscanf(fp,"%c%g",&cs,&xs);
  switch (cs)
  {
    case 'X': if ((xs<x0) || (xs>x0+(Nx-1)*dh)) errquit("Slice position is outside of volume.\n"); break;
    case 'Y': ys=xs; if ((ys<y0) || (ys>y0+(Ny-1)*dh)) errquit("Slice position is outside of volume.\n"); 
              break;
    default: errquit("Incorrect label for slice position. MUST be either X, Y or Z.\n");
  }
  fgets(buf,STRLN,fp);
  if (fp!=stdin) fclose(fp);

  alloc_float_matrix(&V,Ny,Nx);
  read_2D_surface(V,Ny,Nx,fln_in);
  

  fo=fopen(fln_o,"wt");
  switch (cs)
  {
    case 'Y':
      is = gras_round((ys-y0)/dh);
      printf("Yi = %d\n",is);
      for (k=0;k<Nx;k++)
        fprintf(fo,"%g %g\n",x0+k*dh, V[is][k]*(-1));
    break;
    case 'X':
      js = gras_round((xs-x0)/dh);
      printf("Xj = %d\n",js);
      for (k=0;k<Ny;k++)
	fprintf(fo,"%g %g\n",y0+k*dh, V[k][is]*(-1));
    break;
  }
  fclose(fo);
  
  free_float_matrix(V,Ny,Nx);
  return 0;
}
