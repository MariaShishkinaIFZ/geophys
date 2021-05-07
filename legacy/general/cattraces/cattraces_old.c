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

int main(int narg, char **argv)
{
  int Nf, N_rays, N0;
  int i, j;
  int *NR;
  int ec;
  char **filenames;
  char fln[FLNLN];
  FILE *fi, *fo, *fl;

  if (narg<3)
  {
    fprintf(stderr,"Insufficient parameters.\n cattraces <filelist file> <outfile>\n");
    exit(EXIT_FAILURE);
  }
  fl=file_open(argv[1],"rt");
  Nf=0;
  fscanf(fl,"%s",fln);
  while (!feof(fl))
  {
    Nf++;
    fscanf(fl,"%s",fln);
  }
  
  rewind(fl);
  filenames=(char**)fmemalloc(Nf,sizeof(char*));
  NR=(int*)fmemalloc(Nf,sizeof(int));
  for (i=0;i<Nf;i++) filenames[i]=(char*)fmemalloc(FLNLN,sizeof(char));
  
  for (i=0;i<Nf;i++)  fscanf(fl,"%s",filenames[i]);
  fclose(fl);
  
  fo=file_open(argv[2],"wt");
  N_rays=0;
  for (i=0;i<Nf;i++)
  {
    fi=file_open(filenames[i],"rt");
    fprintf(stderr,"%s opened ",filenames[i]);
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
    fprintf(stderr,"%d rays\n",NR[i]);
  }
  
  fprintf(fo,"%d\n",N_rays);
  N0=1;
  for (i=0;i<Nf;i++)
  {
    fi=file_open(filenames[i],"rt");
    fprintf(stderr,"%s opened ",filenames[i]);
    for (j=0;j<NR[i];j++)
      rw_ray(fi,fo,N0+j);
    N0+=NR[i];
    fprintf(stderr," written\n");
  }
  fclose(fo);
  
  for (i=0;i<Nf;i++)
    fmemfree(filenames[i]);
  fmemfree(filenames);
  fmemfree(NR);
  return 1;
}
