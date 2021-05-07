
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
  int id_p, id_tmp, prof_n, prof_p;
  float d, az;
  char s[STRLN], string[LN], s_id[LN], id_s_tmp[LN];

  FILE *fi, *fd, *fd1, *fo;

  if (narg<5)
  {
    printf("Insufficient parameters.\nUSAGE: cnvpeaks <reciver's prof ID  X Y Z file> <source's ID  X Y Z file> <peaks data file> source_ID <OUT peaks.tpk_file>\n");
    exit(1);
  }

  fd=file_open(argv[1],"rt");
  fd1=file_open(argv[2],"rt");
  fi=file_open(argv[3],"rt");
  sscanf(argv[4],"%s",&s_id);
  fo=file_open(argv[5],"wt");


	found=0;
	  do
	  {
	  fgets(string,LN,fd1);
	  sscanf(string,"%s%g%g%g", &id_s_tmp,&xs,&ys,&zs); 
	  
// 	  printf("%s %s\n",s_id,id_s_tmp);
	  
	  if(strcmp(s_id,id_s_tmp)==0)
	    {
	    rewind(fd1);
	    found=1;
	    }
	  }
	  while (!found && !feof(fd1));
        if (!found) {printf("Coordinates of TFO (source) id: %s not found.\n",s_id); errquit("");}
    



  fgets(s,STRLN,fi);
  c=0;
  while (!feof(fi))
  {
    sscanf(s,"%d%d%g",&prof_p,&id_p,&tp);

	found=0;
	  do
	  {
	  fgets(string,LN,fd);
	  sscanf(string,"%d%d%g%g%g",&prof_n, &id_tmp,&xp,&yp,&zp); 
// 	  printf("%d %d \n",id_p,id_tmp);
	  
	  if(id_tmp == id_p)
	    {
	    rewind(fd);
	    found=1;
	    }
	  }
	  while (!found && !feof(fd));
        if (!found) {printf("Coordinates of RTO (reciver) id: %d not found.\n",id_p);} //errquit("");
    
    //printf("%g %g %g\n", xp,yp,zp);
    d=sqrt((xs-xp)*(xs-xp)+(ys-yp)*(ys-yp)+(zs-zp)*(zs-zp));
    az=180.0/M_PI*atan2(xp-xs,yp-ys);

    if (found) {fprintf(fo,"%5d  %5d %10g %10g %10g\n",prof_p,id_p,d,az,tp); c++; }
    fgets(s,STRLN,fi);
  }
  printf("%d peaks done.\n",c);

  fclose(fo);
  fclose(fd);
  fclose(fi);
  return 1;
}
