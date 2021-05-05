#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "subio.h"

#define TEMPFILE "rays.XXXXXX"
#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/
#define STRLN 256

int rwr_ray (FILE *fi, FILE *fo, float *Xaa, float *Yaa, float *Xbb, float *Ybb)
{
  long int nrec,i,j,k,nc,rayN;
  int ncell, raypos;
  float x,y,z,dt,st,ll,q;

  if (fscanf(fi,"%ld%ld%ld%ld%ld%g%g%g%g%g",&rayN,&nrec,&i,&j,&k,&x,&y,&z,&dt,&st)!=10)
  {
    if (feof(fi)) return -1;
    else return -2;
  }
  fprintf(fo,"%ld\n%7ld %5ld %5ld %5ld %13g %13g %13g %13g %13g\n",rayN,nrec,i,j,k,x,y,z,dt,st);
    *Xaa=x;
    *Yaa=y;
  ncell=0;
  do
  {
    nc=ncell;
    if (fscanf(fi,"%d%ld%ld%ld%g%g%g%g%g%d",&ncell,&i,&j,&k,&x,&y,&z,&ll,&q,&raypos)!=10)
      return -2;
    if (ncell!=-1)
    {
      *Xbb=x;
      *Ybb=y;     
    }
    fprintf(fo,"%10d %5ld %5ld %5ld %13g %13g %13g %13g %13g %3d\n",ncell,i,j,k,x,y,z,ll,q,raypos);
  }
  while (ncell!=-1);
  return nc;
}


int filter_ray (FILE *fi, FILE *fo, float Xa, float Ya, float Xb, float Yb, float h)
{
  float Xaa, Yaa, Xbb, Ybb;
  int n, Nr,R_count=0;
  FILE *ftt;
  char *sftt;
  char template[FLNLN];
  
  strcpy(template,TEMPFILE);
    sftt=mktemp(template);
  fscanf(fi,"%d",&Nr);
//   printf("Original: %g %g %g %g \n \n",Xa, Ya, Xb, Yb);
  for (n=0;n<Nr;n++)
  {
   ftt=file_open(sftt,"wt");
   rwr_ray(fi,ftt, &Xaa, &Yaa, &Xbb, &Ybb);
   fclose(ftt);
//   printf("Delta: %g %g %g %g \n", sqrt((Xaa-Xa)*(Xaa-Xa)+(Yaa-Ya)*(Yaa-Ya)), sqrt((Xbb-Xb)*(Xbb-Xb)+(Ybb-Yb)*(Ybb-Yb)), sqrt((Xaa-Xb)*(Xaa-Xb)+(Yaa-Yb)*(Yaa-Yb)), sqrt((Xbb-Xa)*(Xbb-Xa)+(Ybb-Ya)*(Ybb-Ya)));
   
//    if ((abs(Xaa-Xa)<h && abs(Yaa-Ya)<h && abs(Xb-Xb)<h && abs(Ybb-Yb)<h) || (abs(Xaa-Xb)<h && abs(Yaa-Yb)<h && abs(Xb-Xa)<h && abs(Ybb-Ya)<h))
  if ((sqrt((Xaa-Xa)*(Xaa-Xa)+(Yaa-Ya)*(Yaa-Ya))<h &&  sqrt((Xbb-Xb)*(Xbb-Xb)+(Ybb-Yb)*(Ybb-Yb))<h) || (sqrt((Xaa-Xb)*(Xaa-Xb)+(Yaa-Yb)*(Yaa-Yb))<h &&  sqrt((Xbb-Xa)*(Xbb-Xa)+(Ybb-Ya)*(Ybb-Ya))<h))
   {
     ftt=file_open(sftt,"rt");
     rwr_ray(ftt,fo, &Xaa, &Yaa, &Xbb, &Ybb);
     R_count++;
     fclose(ftt);
     remove(sftt);
   }
   remove(sftt);
  }
 printf("%d rays done.\n", R_count);
 return R_count;        
}


int main(int narg, char **argv)
{
  int Nr;
  int n, n_ray, n_seg, i, j, k, cr;
  int refl_flag;
  int si_only=0;
  float x, y, z, xp, yp, zp;
  float l1, l2;
  float Xa, Ya, Xb, Yb;
  char plane,opt;
  FILE *fi, *fo, *ft, *fp;
  char *opts = "l:XYZ"; // доступные опции
  char *sft;
  
  long Nx,Ny,Nz, k_gen, k_gen_v;
  float x0,y0,z0;
  float h; 
  
  
  if (narg==1) errquit("USE: tracegmt -l<Xa/Ya/Xb/Yb>|X|Y|Z|C <ray traces file> <GMT file (.xy)>\n");
  
  while((opt = getopt(narg, argv, opts)) != -1) { // вызываем getopt пока она не вернет -1
        switch(opt) {
	  case 'l':
	    sscanf(strtok(optarg,"/"),"%g",&Xa);
	    sscanf(strtok('\0',"/"),"%g",&Ya);
	    sscanf(strtok('\0',"/"),"%g",&Xb);
	    sscanf(strtok('\0',"/"),"%g",&Yb);
	    plane='l';
	    k = optind;
	    break;
	    
	  case 'X':
	    plane='X';
	    k = optind;
	    break;
	    
	  case'Y':
	    plane='Y';
	    k = optind;
	    break;
	    
	  case 'Z':
	    plane='Z';
	    k = optind;
	    break;
	    
	  case 'C':
	    si_only=1;
	    errquit("Option -C is not described now, we are sorry.\n");
	    k = optind;
	    break;
	    
	  default: errquit("USE: tracegmt -l<Xa/Ya/Xb/Yb>|X|Y|Z|C <ray traces file> <GMT file (.xy)>\n");
	}
  }
  
  fprintf(stderr,"Coordinates Xa= %g; Ya= %g; Xb= %g; Yb= %g\n",Xa,Ya,Xb,Yb);
  fprintf(stderr,"Plane: %c\n",plane);
  fprintf(stderr,"Input file: %s\nOutput file %s \n",argv[k],argv[k+1]);
  
  fi=file_open(argv[k],"rt");
  fo=file_open(argv[k+1],"wt"); 
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fscanf(fp,"%ld",&k_gen);
  fscanf(fp,"%ld",&k_gen_v);
  fclose(fp);

  char template[FLNLN];
  strcpy(template,TEMPFILE);

switch (plane) {
    case 'l':    
      sft=mktemp(template);
      ft=file_open(sft,"wt");
      Nr=filter_ray(fi,ft,Xa,Ya,Xb,Yb,h*k_gen*5);
      fprintf(stderr,"1Nr = %d.\n",Nr);  
      fclose(ft);
      fclose(fi);
      fi=file_open(sft,"rt");
    break;

//     case 'X||Y'||'Z'||'C': 
    default:
      fscanf(fi,"%d",&Nr); 
    break;
   }
   

//    fprintf(stderr,"Nr = %d.\n",Nr);  
  
  for (n=0;n<Nr;n++)
  {
//      fprintf(stderr,"Check point.\n");  
    fscanf(fi,"%d",&n_ray);
    //printf("%d ",n_ray);
    fscanf(fi,"%d%*d%*d%*d%*g%*g%*g%*g%*g",&i); // printf("%d\n",i);
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
// 	   fprintf(stderr,"1 Plane: %c\n",plane);
	   fprintf(fo,"> -W2/0/0/255\n");
	   switch (plane)
           {
             case ('X'): fprintf(fo,"%g %g\n",yp,-zp); break;
	     case ('Y'): fprintf(fo,"%g %g\n",xp,-zp); break;
	     case ('Z'): fprintf(fo,"%g %g\n",xp,yp); break;
	     case ('l'): fprintf(fo,"%g %g\n",sqrt(pow((Xa-xp),2)+pow((Ya-yp),2)),-zp); break;
	     default : fprintf(fo,"%g %g\n",xp,yp); break;
           }
	 }
	 refl_flag=2;
       }
       if ((refl_flag==2) && (cr==3))
       {
         if (!si_only)
	 {
// 	   fprintf(stderr,"2 Plane: %c\n",plane);
	   fprintf(fo,"> -W2/0/255/0\n");
	   switch (plane)
           {
             case ('X'): fprintf(fo,"%g %g\n",yp,-zp); break;
	     case ('Y'): fprintf(fo,"%g %g\n",xp,-zp); break;
	     case ('Z'): fprintf(fo,"%g %g\n",xp,yp); break;
	     case ('l'): fprintf(fo,"%g %g\n",sqrt(pow((Xa-xp),2)+pow((Ya-yp),2)),-zp); break;
	     default : fprintf(fo,"%g %g\n",xp,yp); break;
           }
	 }
	 refl_flag=3;
       }
       if ((!si_only) || (cr==2))
       {
         switch (plane)
         {
// 	   fprintf(stderr,"3 Plane: %c\n",plane);
           case ('X'): fprintf(fo,"%g %g\n",y,-z); break;
	   case ('Y'): fprintf(fo,"%g %g\n",x,-z); break;
	   case ('Z'): fprintf(fo,"%g %g\n",x,y); break;
	   case ('l'): fprintf(fo,"%g %g\n",sqrt(pow((Xa-x),2)+pow((Ya-y),2)),-z); break;
	   default : fprintf(fo,"%g %g\n",x,y); break;
         }
       }
      }
    }
    while (n_seg!=-1);  
  }
  
  fclose(fo);
  fclose(fi);
  remove(sft);
  
  return EXIT_SUCCESS;
}
