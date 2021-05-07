#include <cstdio>
#include "gstio.h"
#include "gstdata.h"
#include "gsterror.h"
#include "gmt_nan.h"

int main(int narg, char **argv)
{
  int Nr, Nc;
  int i,j,c;
  double s;
  GSTgriddedData *grd;
  
  
  if (narg < 3)
    GSTerror("Usage: pc_by_row <infile.gmt> <outfile.gmt>\n");
  
  grd=readGMT(argv[1]);
  grd->size(Nr,Nc);
  
  for (i=0;i<Nr;i++)
  {
    s=0.0; c=0;
    for (j=0;j<Nc;j++)
      if (!GMT_is_fnan(grd->data()[i][j]))
      {
        c++;
	s+=grd->data()[i][j];
      }
    if (c!=0)
    {
      fprintf(stderr,"%d %d %g\n",i,c,s);
      s/=c;
      for (j=0;j<Nc;j++)
	if (!GMT_is_fnan(grd->data()[i][j]))
	  grd->setValue(i,j,(grd->data()[i][j]-s)/s*100.0);
    }
  }
  writeGMT(argv[2],*grd);
  delete grd;
}
