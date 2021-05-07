#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "submem.h"
#include "subio.h"
#include "iogras.h"

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/

#define EMPTY_VALUE 1e10
#define EMPTY_LIMIT 1e9

/*x0 - начало отсчёта, w - ширина клетки, t - полуширина переходной зоны*/
float ch1d(float x, float x0, float w, float t)
{
  long n = floor((x-x0)/w); /*номер клетки*/
  long a = 1;
  
  if (n==0)
  {
    if ((x-x0>=n*w+t) && (x-x0<=(n+1)*w-t)) return a;
    else
    {
      if (x-x0<n*w+t) return a*sin(M_PI*(x-x0-n*w)/(2.0*t));
      else return a*sin(M_PI*((n+1)*w-(x-x0))/(2.0*t));
    }
  }
  else
    return EMPTY_VALUE;
}

float ch1d1(float x, float x0)
{
  long n = floor(x-x0); /*номер клетки*/
  long a = 1;
  
  if (n==0)
  return a;
  else
    return EMPTY_VALUE;
}

int main(int narg, char **argv)
{
  char *par_file;
  long int Nx, Ny, Nz; 
  int i, j, k, p, pp, ppp;
  long Ix, Iy, Iz;
  float al, xrot, yrot;
  float x0, y0, z0, h;
  float cx0, cxw, cxt, cy0, cyw, cyt, cz0, czw, czt, cx1, cy1, cz1, cx0_rot, cy0_rot;
  float abs_vel, rel_anom;
  float ***V;
  float incl_X_rot, incl_Y_rot, incl_Z;
  char fout[FLNLN], vm_file[FLNLN];
  FILE *fp;
  char *opts = "rsf:"; //flags
  char opt;
  int r_flag=0, s_flag=0;
  
  if (narg==1) errquit("USE: ins_incl_rot_rel [-r] (-f <par_file> | -s) \n -r - require relative amplitude of inclusion (for negative anomaly use minus), else - absolute amplitude of inclusion \n -f <par_file> use parameteres from file \n -s - use stdin to input parameters \n");
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%ld%ld%ld",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
  
   while((opt = getopt(narg, argv, opts)) != -1) { // вызываем getopt пока она не вернет -1 
  
     switch(opt) {
    //case 'a': // if flag -a
     // a_flag=1;
     // k=optind;
   // break;
    
    case 'r': // 
       r_flag=1; 
       k = optind;
    break;
    
    case 'f': // }
      sscanf(optarg,"%s",par_file);
      if (!(fp=file_open(par_file,"rt"))) errquit("Parameters file name not given.\n");
      fmemfree(par_file);
      k = optind; 
    break;
	
    case 's': //
      s_flag=1;
      fp=stdin;
      k = optind;
    break;   }

 }
 
  if(!r_flag) {
  
  if (s_flag) fprintf(stdout,"Please enter velocity model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",vm_file);
  
  if (s_flag) fprintf(stdout,"Please enter x0, x1, X wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cx0,&cx1,&cxt); cxw=cx1-cx0;
  if (s_flag) fprintf(stdout,"Please enter y0, y1, Y wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cy0,&cy1,&cyt); cyw=cy1-cy0;
  if (s_flag) fprintf(stdout,"Please enter z0, z1, Z wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cz0,&cz1,&czt); czw=cz1-cz0;
  if (s_flag) fprintf(stdout,"Please enter alpha (XY plane rotation angle, degrees)  > \n");
  fscanf(fp,"%g",&al);
  al*=M_PI/180.0;
  if (s_flag) fprintf(stdout,"Please enter the absolute velocity of the inclusion > \n");
  fscanf(fp,"%g",&abs_vel);
  if (s_flag) fprintf(stdout,"Please enter chessed model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",fout);
  if (!s_flag) fclose(fp);
  
  alloc_float_3d(&V,Nz,Ny,Nx);
  read_3D_volume(V,Nz,Ny,Nx,vm_file);
  cx0_rot=cx0*cos(al)+cy0*-sin(al);
  cy0_rot=cx0*sin(al)+cy0*cos(al);
  
  for (k=0;k<Nz;k++)
    for (i=0;i<Ny;i++)
      for (j=0;j<Nx;j++)
      { Ix=ch1d1(x0+h*j,cx0);
	Iy=ch1d1(y0+h*i,cy0);
	Iz=ch1d1(z0+h*k,cz0);
	if (Ix<EMPTY_LIMIT && Iy<EMPTY_LIMIT && Iz<EMPTY_LIMIT)
        fprintf(stdout,"Rotation cx0, cy0, cz0 coordinates:\n %g %g %g\n",x0+h*j,y0+h*i,z0+h*k);
            for (p=k+1;p<Nz;p++)
             for (pp=i+1;pp<Ny;pp++)
               for (ppp=j+1;ppp<Nx;ppp++)
      {
	xrot=(x0+h*ppp-x0-h*j)*cos(al)+(y0+h*pp-y0-h*i)*-sin(al)+x0+h*j;
	yrot=(x0+h*ppp-x0-h*j)*sin(al)+(y0+h*pp-y0-h*i)*cos(al)+y0+h*i;
        incl_Z=ch1d(z0+h*p,cz0,czw,czt);
        incl_Y_rot=ch1d(yrot,x0+h*j,cyw,cyt);
        incl_X_rot=ch1d(xrot,y0+h*i,cxw,cxt);

        if (incl_X_rot<EMPTY_LIMIT && incl_Y_rot<EMPTY_LIMIT && incl_Z<EMPTY_LIMIT)
        V[p][pp][ppp]+=(abs_vel-V[p][pp][ppp])*incl_X_rot*incl_Y_rot*incl_Z;
      }
      
  write_3D_volume(V,Nz,Ny,Nx,fout);
  free_float_3d(V,Nz,Ny,Nx);
  exit (1);
}
      }
      
if(r_flag) {
  
  if (s_flag) fprintf(stdout,"Please enter velocity model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",vm_file);
  
  if (s_flag) fprintf(stdout,"Please enter x0, x1, X wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cx0,&cx1,&cxt); cxw=cx1-cx0;
  if (s_flag) fprintf(stdout,"Please enter y0, y1, Y wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cy0,&cy1,&cyt); cyw=cy1-cy0;
  if (s_flag) fprintf(stdout,"Please enter z0, z1, Z wide of transit to next cell  > \n");
  fscanf(fp,"%g%g%g",&cz0,&cz1,&czt); czw=cz1-cz0;
  if (s_flag) fprintf(stdout,"Please enter alpha (XY plane rotation angle, degrees)  > \n");
  fscanf(fp,"%g",&al);
  al*=M_PI/180.0;
  if (s_flag) fprintf(stdout,"Please enter the relative velocity of the inclusion > \n");
  fscanf(fp,"%g",&rel_anom);
  if (s_flag) fprintf(stdout,"Please enter chessed model filename *.3d (%d chars max) > \n",FLNLN);
  fscanf(fp,"%s",fout);
  if (!s_flag) fclose(fp);
  
  alloc_float_3d(&V,Nz,Ny,Nx);
  read_3D_volume(V,Nz,Ny,Nx,vm_file);
  //cx0_rot=cx0*cos(al)+cy0*-sin(al);
  //cy0_rot=cx0*sin(al)+cy0*cos(al);
  
  for (k=0;k<Nz;k++)
    for (i=0;i<Ny;i++)
      for (j=0;j<Nx;j++)
      {
	xrot=(x0+h*j)*cos(al)+(y0+h*i)*sin(al)-cx0;
	yrot=(x0+h*j)*-sin(al)+(y0+h*i)*cos(al)-cy0;
        incl_Z=ch1d(z0+h*k,cz0,czw,czt);
        incl_Y_rot=ch1d(yrot,cy0,cyw,cyt);
        incl_X_rot=ch1d(xrot,cx0,cxw,cxt);
	
        if (incl_X_rot<EMPTY_LIMIT && incl_Y_rot<EMPTY_LIMIT && incl_Z<EMPTY_LIMIT)
	  
        V[k][i][j]=V[k][i][j]*(1+rel_anom*incl_X_rot*incl_Y_rot*incl_Z);
      }
      
  write_3D_volume(V,Nz,Ny,Nx,fout);
  free_float_3d(V,Nz,Ny,Nx);
  exit (1);
}

}