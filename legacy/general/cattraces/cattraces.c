#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <values.h>

#include "subio.h"
#include "submem.h"

#include "gmtio.h"

#define FLNLN 80 	/*Maximum length of the file name*/

int skip_ray(FILE *f)
{
  long int nrec,i,j,k,nc,rayN;
  int ncell, raypos;
  float x,y,z,dt,st,ll,q;

  if (fscanf(f,"%ld%ld%ld%ld%ld%g%g%g%g%g",&rayN,&nrec,&i,&j,&k,&x,&y,&z,&dt,&st)!=10)
  {
    //fprintf(stderr,"%ld\n%7ld %5ld %5ld %5ld %13g %13g %13g %13g %13g\n",rayN,nrec,i,j,k,x,y,z,dt,st);
    if (feof(f)) return -1;
    else return -2;
  }

  ncell=0;
  do
  {
    nc=ncell;
    if (fscanf(f,"%d%ld%ld%ld%g%g%g%g%g%d",&ncell,&i,&j,&k,&x,&y,&z,&ll,&q,&raypos)!=10)
      return -2;
  }
  while (ncell!=-1);
  return nc;
}

int rw_ray(FILE *fi, FILE *fo, long int nr)
{
  long int nrec,i,j,k,nc,rayN;
  int ncell, raypos;
  float x,y,z,dt,st,ll,q;

  if (fscanf(fi,"%ld%ld%ld%ld%ld%g%g%g%g%g",&rayN,&nrec,&i,&j,&k,&x,&y,&z,&dt,&st)!=10)
  {
    if (feof(fi)) return -1;
    else return -2;
  }
  fprintf(fo,"%ld\n%7ld %5ld %5ld %5ld %13g %13g %13g %13g %13g\n",nr,nrec,i,j,k,x,y,z,dt,st);

  ncell=0;
  do
  {
    nc=ncell;
    if (fscanf(fi,"%d%ld%ld%ld%g%g%g%g%g%d",&ncell,&i,&j,&k,&x,&y,&z,&ll,&q,&raypos)!=10)
      return -2;
    fprintf(fo,"%10d %5ld %5ld %5ld %13g %13g %13g %13g %13g %3d\n",ncell,i,j,k,x,y,z,ll,q,raypos);
  }
  while (ncell!=-1);
  return nc;
}

int rw_point(FILE *fi, FILE *fo, long int nr)
{
  long int i,j,pointN;
  float x,y,z,dt,st;
  
  
  if (fscanf(fi,"%ld%ld%ld%g%g%g%g%g",&pointN,&i,&j,&x,&y,&z,&dt,&st)!=8)
  {
    if (feof(fi)) return -1;
    else return -2;
  }
  fprintf(fo,"%4ld\n",nr);
  fprintf(fo,"%4ld %6ld %9g %9g %9g %13g %13g\n",i,j,x,y,z,dt,st);  
  
  if (fscanf(fi,"%ld%ld%g%g%g%g%g",&i,&j,&x,&y,&z,&dt,&st)!=7) return -2;
  fprintf(fo,"%4ld %6ld %9g %9g %9g %13g %13g\n",i,j,x,y,z,dt,st);
    
  return 2;
}

int main(int narg, char **argv)
{
  int Nf, Nf_H, Nf_D, N_rays, N0;
  char T_flag[FLNLN];
  int i, j, ih, id;
  int *NR;
  int ec;
  char **filenames, **pointfilenames;
  char fln[FLNLN];
  FILE *fi, *fo, *fl;

  if (narg<3)
  {
    fprintf(stderr,"Insufficient parameters.\n cattraces <filelist file> <outfile>\n");
    exit(EXIT_FAILURE);
  }
  fl=file_open(argv[1],"rt");
  Nf=0;
  Nf_H=0;
  Nf_D=0;
  
  while (!feof(fl))
  {
    fgets(fln,160,fl);
    if (!feof(fl))
    {
      Nf++;
      if (fln[0] == 'H') Nf_H++; 
      if (fln[0] == 'D') Nf_D++; 
    }
  }

  
  fprintf(stderr,"Paramets file consits of %d strings.\n %d are nominated as head waves files\n %d are nominated as diving waves files\n",Nf, Nf_H, Nf_D);
      
  if (narg<4&&Nf_H!=0)
  {
    fprintf(stderr,"Insufficient parameters\nNumber of files nominated as head waves is not 0 .\n cattraces <filelist file> <raytraces_outfile> <ref_points_outfile>\n");
    exit(EXIT_FAILURE);
  }
      
  rewind(fl);
  filenames=(char**)fmemalloc((Nf_H+Nf_D),sizeof(char*));
  pointfilenames=(char**)fmemalloc((Nf_H),sizeof(char*));  
  NR=(int*)fmemalloc((Nf_H+Nf_D),sizeof(int));
  for (i=0;i<(Nf_H+Nf_D);i++) filenames[i]=(char*)fmemalloc(FLNLN,sizeof(char));
  for (i=0;i<(Nf_H);i++) pointfilenames[i]=(char*)fmemalloc(FLNLN,sizeof(char));
  
  ih=0;
  id=0;
  //fprintf(stderr,"Reading filenames from file_list. . .\n");
  while (!feof(fl))
  {
    if (fscanf(fl,"%s",&T_flag)==1)
    {
     //fprintf(stderr,"%s\n",T_flag);
     if (T_flag[0] == 'H') 
      {
      fscanf(fl,"%s%s",filenames[ih],pointfilenames[ih]);
      //fprintf(stderr,"%d H type file_nanme %s is read from file_list\n",ih,filenames[ih]);
      ih++;
      }
    
     if (T_flag[0] == 'D') 
      {
      fscanf(fl,"%s",filenames[id+Nf_H]);
      //fprintf(stderr,"%d D type file_nanme %s is read from file_list, name is writen to %d position in memory\n",id,filenames[id+Nf_H],(id+Nf_H));
      id++;
      } 
    }
  }
  fclose(fl);
  
  N_rays=0;
  for (i=0;i<(Nf_H+Nf_D);i++)
  {
    
    fi=file_open(filenames[i],"rt");
    fprintf(stderr,"File: %s opened",filenames[i]);
    NR[i]=0;
    do
    {
      ec=skip_ray(fi);
      //fprintf(stderr,"%d %d\n",ec,NR[i]);
      if (ec==-2)
      {
        fprintf(stderr,"Error reading ray\n");
        exit(EXIT_FAILURE);
      }
      if (ec!=-1) NR[i]++;
    }
    while (ec!=-1);
    fclose(fi);
    N_rays+=NR[i];
    fprintf(stderr,", consists of %d rays\n",NR[i]);
  }
  fprintf(stderr,"All rays in all files are readable and there are %d rays \n",N_rays);

  fo=file_open(argv[2],"wt");

  fprintf(fo,"%d\n",N_rays);
  N0=1;
  for (i=0;i<(Nf_H+Nf_D);i++)
  {
    fi=file_open(filenames[i],"rt");
    fprintf(stderr,"File: %s opened ",filenames[i]);
    for (j=0;j<NR[i];j++)
      rw_ray(fi,fo,N0+j);
    N0+=NR[i];
    fprintf(stderr,"and written to: %s \n", argv[2]);
  }
  fclose(fo);
  
  if(Nf_H!=0)
  {
    fo=file_open(argv[3],"wt");
  
    N0=0;
    for (i=0;i<(Nf_H);i++)
    {
    N0+=NR[i];
    }
  
    fprintf(fo,"%d\n",N0);
    N0=1;
    for (i=0;i<(Nf_H);i++)
    {
      fi=file_open(pointfilenames[i],"rt");
      fprintf(stderr,"Point file: %s opened ",pointfilenames[i]);
      for (j=0;j<NR[i];j++)
        rw_point(fi,fo,N0+j);
      N0+=NR[i];
      fprintf(stderr,"and written to: %s \n", argv[3]);
    }
    fclose(fo);
    fprintf(stderr,"*******%d points are renumbered and writen to: %s *******\n", (N0-1), argv[3]);
  }
  
  
  for (i=0;i<(Nf_H+Nf_D);i++) fmemfree(filenames[i]);
  for (i=0;i<(Nf_H);i++) fmemfree(pointfilenames[i]);
  fmemfree(filenames);
  fmemfree(pointfilenames);
  fmemfree(NR);
  
  fprintf(stderr,"*******%d rays are renumbered and writen to: %s *******\n", N_rays, argv[2]);
  
  return 1;
}
