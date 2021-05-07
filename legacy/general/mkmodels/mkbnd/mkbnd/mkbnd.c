#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "submem.h"
#include "subio.h"
#include <unistd.h>

#define GEOMFILE "model-geometry.gras"
#define FLNLN 256 	/*Maximum length of the file name*/



int main(int narg, char **argv)
{
  int Nx, Ny, Nz, i, j, k;
  char *par_file;
  float x0, y0, z0, h;
  float a1, a2, a3;
  float b1, b2, b3;
  float c1, c2, c3;
  float nx, ny, nz;
  float alfa, delta;
  float **Int;
  char out_file[FLNLN];
  FILE *fp;
  char *opts = "acf:s"; //flags
  char opt;
  int c_flag=0, f_flag=0, s_flag=0;
 

  
  if (narg==1) errquit("USE: mkbnd [-c] (-f <par_file> | -s) \n");
    
  
  fp=file_open(GEOMFILE,"rt");
  fscanf(fp,"%d%d%d",&Nx,&Ny,&Nz);
  fscanf(fp,"%g%g%g",&x0,&y0,&z0);
  fscanf(fp,"%g",&h);
  fclose(fp);
  
  //par_file = fmemalloc(FLNLN,sizeof(char));
  
  while((opt = getopt(narg, argv, opts)) != -1) { // вызываем getopt пока она не вернет -1 
  switch(opt) {
    //case 'a': // if flag -a
     // a_flag=1;
     // k=optind;
   // break;
    
    case 'c': // 
       c_flag=1; 
       k = optind;
    break;
    
    case 'f': // }
       f_flag=1;
       
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
  
  //if (f_flag) /*Use parameters file with the name given as 2nd parameter*/
  //{
    //fp=file_open(argv[k],"rt");
  //}
  
  
  //if (s_flag) fp=stdin;
  
  if(c_flag) {
  if (s_flag) fprintf(stdout,"Input x y z, - coordinates of 1st point \n");
  fscanf(fp,"%g%g%g",&a1,&a2,&a3);
  
  if (s_flag) fprintf(stdout,"Input x y z - coordinates of 2nd point \n");
  fscanf(fp,"%g%g%g",&b1, &b2, &b3);
  
  if (s_flag) fprintf(stdout,"Input x y z - coordinates of 3d point \n");
  fscanf(fp,"%g%g%g",&c1, &c2, &c3);
  
  if (s_flag) fprintf(stdout,"Input boundary model filename.2d (%d chars max) \n",FLNLN);
  fscanf(fp,"%s",out_file);
  
  if (!s_flag) fclose(fp);
  
     nx=(b2-a2)*(c3-a3)-(b3-a3)*(c2-a2);
  
     ny=(b1-a1)*(c3-a3)-(b3-a3)*(c1-a1);
  
     nz=(b1-a1)*(c2-a2)-(b2-a2)*(c1-a1);
 

  alloc_float_matrix(&Int,Ny,Nx);
 
  for (i=0; i<Ny; i++){
    for (j=0; j<Nx; j++){ 
	Int[i][j]=(nz*a3-nx*(j*h-a1)+ny*(i*h-a2))/nz;
    }
  }

  write_2D_surface(Int,Ny,Nx,out_file);

  free_float_matrix(Int,Ny,Nx);
}
  
  else {
    if (s_flag) fprintf(stdout,"Input x y z, - coordinates \n");
       fscanf(fp,"%g%g%g",&a1,&a2,&a3);
  
    if (s_flag) fprintf(stdout,"Input hade and strech angles \n");
       fscanf(fp,"%g%g",&alfa, &delta);
    
    if (s_flag) fprintf(stdout,"Input boundary model filename.2d (%d chars max) \n",FLNLN);
       fscanf(fp,"%s",out_file);
    
    if (!s_flag) fclose(fp);
 
    if (delta>=90)
    {
    fprintf(stderr,"Hade should be  in interval between 0 and 90 degrees!\n");
    exit(0);
    }

    if (delta<0)
    {
    fprintf(stderr,"Hade should be  in interval between 0 and 90 degrees!\n");
    exit(0);
    }

  alfa =(alfa*(M_PI/180.0));
  delta=(delta*(M_PI/180.0));

  nx=sin(delta)*sin(alfa);
  ny=sin(delta)*cos(alfa);
  nz=-cos(delta);

  alloc_float_matrix(&Int,Ny,Nx);
 
  for (i=0; i<Ny; i++){
    for (j=0; j<Nx; j++){ 
	Int[i][j]=(nz*a3-nx*(j*h-a1)-ny*(i*h-a2))/nz;
    }
  }

  write_2D_surface(Int,Ny,Nx,out_file);

  free_float_matrix(Int,Ny,Nx);
}


return 0;

}

  
  
  
  
