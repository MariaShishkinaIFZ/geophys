#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "gstdata.h"
#include "gstio.h"

int main(int narg, char **argv)
{
  int Nr,Nc,i,j;
  double x1,x2,y1,y2;
  GSTgrid_type t;
  GSTgriddedData *D;
  FILE *fo;
  
  if (narg<3)
    GSTerror("Insuffisient papameters. Use:\ngmt2d <GMT grid file> <2d file>\n");
  
  D=readGMT(argv[1]);
  
  D->size(Nr,Nc);
  D->getArea(x1,x2,y1,y2,t);
  fprintf(stdout,"%d over %d grid\n",Nr,Nc);
  fprintf(stdout,"Rows limits: %lg - %lg\n",x1,x2);
  fprintf(stdout,"Columns limits: %lg - %lg\n",y1,y2);
  switch (t)
  {
    case node:
      fprintf(stdout,"Values are attributed to grid nodes.\n"); break;
    case area:
      fprintf(stdout,"Values are attributred to grid cells.\n");
  }
  
  fo=fopen(argv[2],"wb");
  for (i=0;i<Nr;i++)
    fwrite(D->data()[i],sizeof(float),Nc,fo);
  fclose(fo);
}

