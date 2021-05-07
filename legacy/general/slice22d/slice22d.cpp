#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "gstdata.h"
#include "gstio.h"
#include "isdgio.h"
#include "submem.h"
#include "subio.h"

int main(int narg, char **argv)
{
  int Nr,Nc,Ny,i,j,k;
  double x1,x2,y1,y2,dx,dy;
  GSTgrid_type t;
  GSTgriddedData *D;
  float ***buff3D;
  FILE *fo;
  
  if (narg<4)
    GSTerror("Insuffisient papameters. Use:\ngmt2d <GMT grid file> <3d file> Ny\n");
  
  D=readGMT(argv[1]);
  
  D->size(Nr,Nc);
  D->getArea(x1,x2,y1,y2,t);
  D->getIntervals(dx,dy);
  switch (t)
  {
    case area:
      fprintf(stderr,"Warning: Values are attributred to grid cells.\n");
  }
  if (dx!=dy)
    fprintf(stderr,"dx: %lg != dy: %lg This data file does not appear to be a slice from 3D gras data file.\n",dx,dy);
  
  sscanf(argv[3],"%d",&Ny);
  
  alloc_float_3d(&buff3D,Nr,Ny,Nc);
  for (i=0;i<Nc;i++)
    for (k=0;k<Nr;k++)
      for (j=0;j<Ny;j++)
	buff3D[k][j][i]=D->data()[Nr-1-k][i];
      
  delete D;

  write_3D_volume(buff3D,Nr,Ny,Nc,".",argv[2]);
  
  fo = file_open("model-geometry.gras","wt");
  fprintf(fo,"%d %d %d\n",Nc,Ny,Nr);
  fprintf(fo,"%lg %lg %lg\n",y1,-dx*(Ny/2),-x2);
  fprintf(fo,"%lg",dx);
  fclose(fo);
}

