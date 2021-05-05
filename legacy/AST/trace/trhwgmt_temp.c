#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#include "subio.h"

int main(int narg, char **argv)
{
  int Nr;
  int n, n_ray, n_seg, i, j, k, cr;
  int refl_flag;
  int si_only=0;
  float x, y, z, xp, yp, zp;
  float l1, l2;
  char plane;
  FILE *fi, *fo;
  char *opts = "s:XYZ"; // доступные опции
  
  
  if (narg==1) errquit("USE: tracegmt -l< Xa,Ya,Xb,Yb>|X|Y|Z|C <ray traces file> <GMT file (.xy)>\n");
  
  while((opt = getopt(narg, argv, opts)) != -1) { // вызываем getopt пока она не вернет -1
        switch(opt) {
	  case 'l':
	    Xa = optarg[0];
	    Ya = optarg[1];
	    Xb = optarg[2];
	    Yb = optarg[3];
	    k = optind;
	    break;
	    
	  case'X':
	    plane=X;
	    k = optind;
	    break;
	    
	  case'Y':
	    plane=Y;
	    break;
	    
	  case'Z':
	    plane=Z;
	    k = optind;
	    break;
	    
	  case'C':
	    si_only=1
	    errquit("Option -C is not described now, we are sorry.\n");
	    k = optind;
	    break;
	    
	  default errquit("USE: tracegmt -l< Xa,Ya,Xb,Yb>|X|Y|Z|C <ray traces file> <GMT file (.xy)>\n");
	}
  }
  
  
  
  
  fi=file_open(argv[k+1],"rt");
  fo=file_open(argv[k+2],"wt");
  
  
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
