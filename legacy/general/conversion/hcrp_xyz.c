#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "subio.h"

int main(int narg, char **argv)
{
  int Ncrp,m,nray;
  float x,y,z;
  FILE *fi, *fo;
  
  
  if (narg<3) errquit("USE: hcrp_xyz <crosspoints file> <xys file>\n");
  fi=file_open(argv[1],"rt");
  fo=file_open(argv[2],"wt");
  fscanf(fi,"%u",&Ncrp);
  fprintf(stderr,"Number of basement crosspoints pairs: %u\n",Ncrp);
  for (m=1;m<=Ncrp;m++)
  {
    fscanf(fi,"%u",&nray);
    /*First crosspoint*/
    fscanf(fi,"%*d%*d%g%g%g%*g%*g",&x,&y,&z);
    fprintf(fo,"%13g %13g %13g\n",x,y,z);
    fscanf(fi,"%*d%*d%g%g%g%*g%*g",&x,&y,&z);
    fprintf(fo,"%13g %13g %13g\n",x,y,z);
  }
  fclose(fo);
  fclose(fi);
  fprintf(stderr,"Crosspoints information read.\n");
  return EXIT_SUCCESS;
}
