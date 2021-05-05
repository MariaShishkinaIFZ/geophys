#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "subio.h"

int flag_present(char **argv, int n0, int n, char flag)
{
  int i;
  char c;
  for (i=n0;i<n;i++)
  {
    if ((strlen(argv[i])>1) && (argv[i][0]=='-')) /*flag*/
    {
      sscanf(argv[i]+1,"%c",&c);
      if (c==flag) return i;
    }
  }
  return -1; /*flag not found*/
}

int main(int narg, char **argv)
{
  int Nr;
  int n, n_ray, n_seg, i, j, k, cr;
  int refl_flag;
  int si_only;
  float x, y, z, xp, yp, zp;
  float l1, l2;
  char plane;
  FILE *fi, *fo;
  
  if (narg<4) errquit("USE: tracegmt <ray traces file> <GMT file (.xy)> X|Y|Z [-C]\n");
  
  fi=file_open(argv[1],"rt");
  fo=file_open(argv[2],"wt");
  
  sscanf(argv[3],"%c",&plane);
  
  if (flag_present(argv,4,narg,'C')!=-1) si_only=1;
  else si_only=0;
  
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
       if ((refl_flag==1) && (cr==1))
       {
         if (!si_only)
	 {
	   fprintf(fo,"> -W2/0/0/255\n");
	   switch (plane)
           {
             case ('X'): fprintf(fo,"%g %g\n",yp,-zp); break;
	     case ('Y'): fprintf(fo,"%g %g\n",xp,-zp); break;
	     case ('Z'): fprintf(fo,"%g %g\n",xp,yp); break;
	     default : fprintf(fo,"%g %g\n",xp,yp); break;
           }
	 }
	 refl_flag=2;
       }
       if ((refl_flag==2) && (cr==3))
       {
         if (!si_only)
	 {
	   fprintf(fo,"> -W2/0/255/0\n");
	   switch (plane)
           {
             case ('X'): fprintf(fo,"%g %g\n",yp,-zp); break;
	     case ('Y'): fprintf(fo,"%g %g\n",xp,-zp); break;
	     case ('Z'): fprintf(fo,"%g %g\n",xp,yp); break;
	     default : fprintf(fo,"%g %g\n",xp,yp); break;
           }
	 }
	 refl_flag=3;
       }
       if ((!si_only) || (cr==2))
       {
         switch (plane)
         {
           case ('X'): fprintf(fo,"%g %g\n",y,-z); break;
	   case ('Y'): fprintf(fo,"%g %g\n",x,-z); break;
	   case ('Z'): fprintf(fo,"%g %g\n",x,y); break;
	   default : fprintf(fo,"%g %g\n",x,y); break;
         }
       }
      }
    }
    while (n_seg!=-1);  
  }
  
  fclose(fo);
  fclose(fi);
  
  return EXIT_SUCCESS;
}
