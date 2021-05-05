#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "subio.h"

int main(int narg, char **argv)
{
  int Nr;
  int n, n_ray, n_seg, i, j, k, cr;
  int refl_flag;
  float x, y, z, xp, yp, zp;
  float l1, l2;
  char plane;
  FILE *fi, *fo;
  
  if (narg<4) errquit("USE: tracegmt <ray traces file> <GMT file (.xy)> X|Y|Z\n");
  
  fi=file_open(argv[1],"rt");
  fo=file_open(argv[2],"wt");
  
  sscanf(argv[3],"%c",&plane);
  
  fscanf(fi,"%d",&Nr);
  for (n=0;n<Nr;n++)
  {
    fscanf(fi,"%d",&n_ray);
    printf("%d ",n_ray);
    fscanf(fi,"%d%*d%*d%*d%*g%*g%*g%*g%*g",&i); printf("%d\n",i);
    refl_flag=1;
    fprintf(fo,"> -W2/255/0/0\n");
    do
    {
      xp=x; yp=y; zp=z;
      fscanf(fi,"%d%d%d%d%g%g%g%g%g%d",&n_seg,&i,&j,&k,&x,&y,&z,&l1,&l2,&cr);
      if (n_seg>0)
      { 
       if ((refl_flag==1) && (cr==0))
       {
         fprintf(fo,"> -W2/0/0/255\n");
	 switch (plane)
         {
           case ('X'): fprintf(fo,"%g %g\n",yp,-zp); break;
	   case ('Y'): fprintf(fo,"%g %g\n",xp,-zp); break;
	   case ('Z'): fprintf(fo,"%g %g\n",xp,yp); break;
	   default : fprintf(fo,"%g %g\n",xp,yp); break;
         }
	 refl_flag=0;
       }
       switch (plane)
       {
        case ('X'): fprintf(fo,"%g %g\n",y,-z); break;
	case ('Y'): fprintf(fo,"%g %g\n",x,-z); break;
	case ('Z'): fprintf(fo,"%g %g\n",x,y); break;
	default : fprintf(fo,"%g %g\n",x,y); break;
       }
      }
    }
    while (n_seg!=-1);  
  }
  
  fclose(fo);
  fclose(fi);
  
  return EXIT_SUCCESS;
}
