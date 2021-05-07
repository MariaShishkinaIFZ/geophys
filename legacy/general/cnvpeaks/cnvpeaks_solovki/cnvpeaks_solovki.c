
//Convertor desined for Kapustyan data 
//source.xyz 
//reciver.xyz 
//and pairs files for each source: rec_id time 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "subio.h"

#define STRLN 256
#define LN 160		/*Maximum length of string in file*/

int gras_round(float x)
{
  return x-floor(x) >= 0.5 ? ceil(x) : floor(x);
}

int main(int narg, char **argv)
{
  int found, c;
  float xs,ys,zs;
  float right_boarder, left_boarder;
  float right_cut, left_cut;
  float xp, yp, zp, tp;
  int id_p, id_tmp, prof_n;
  float d, az;
  char s[STRLN], string[LN];

  FILE *fi, *fd, *fo;

  if (narg<7)
  {
    printf("Insufficient parameters.\nUSAGE: cnvpeaks <reciver idXYZ file> <peaks data file> <peaks.tpk_file> xs ys zs \n");
    exit(1);
  }

  fd=file_open(argv[1],"rt");
  fi=file_open(argv[2],"rt");
  fo=file_open(argv[3],"wt");


  sscanf(argv[4],"%g",&xs);
  sscanf(argv[5],"%g",&ys);
  sscanf(argv[6],"%g",&zs);


  fgets(s,STRLN,fi);
  c=0;
  while (!feof(fi))
  {
    sscanf(s,"%d%g",&id_p,&tp);

	found=0;
	  do
	  {
	  fgets(string,LN,fd);
	  sscanf(string,"%d%d%g%g%g",&prof_n, &id_tmp,&xp,&yp,&zp); 
	  
	  if(id_tmp == id_p)
	    {
	    rewind(fd);
	    found=1;
	    }
	  }
	  while (!found && !feof(fd));
        if (!found) errquit("Reciver coordinates not found ");  
    
    //printf("%g %g %g\n", xp,yp,zp);
    d=sqrt((xs-xp)*(xs-xp)+(ys-yp)*(ys-yp)+(zs-zp)*(zs-zp));
    az=180.0/M_PI*atan2(xp-xs,yp-ys);


 fprintf(fo,"%5d  %5d %10g %10g %10g\n",prof_n,id_p,d,az,tp);

    c++;
    fgets(s,STRLN,fi);
  }
  printf("%d peaks done.\n",c);

  fclose(fo);
  fclose(fd);
  fclose(fi);
  return 1;
}
