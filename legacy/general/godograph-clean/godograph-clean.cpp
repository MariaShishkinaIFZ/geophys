#include <cstdio>
#include <vector>
#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "subio.h"
#include "submem.h"

#define STRLN 256

using namespace std;

std::vector<float> *average(std::vector<float> *vec, int N, int L)
{
  int i,j,c;
  float x;
  std::vector<float> *vec_av = new vector<float>;

  for (i=0; i<N; i++)
  {
    x=0;
    c=0;
    for (j=i-L/2;j<=i+L/2;j++)
    {
      if (j>=0 && j<N)
      {
        x+=(*vec)[j];
        c++;
      }
    }
    if (c>0) x/=c;
    vec_av->push_back(x);
  }
  return vec_av;
}

int is_slash(char ch) {
    char *sc, *sc0;
    sc = sc0 = (char*)fmemalloc(4,sizeof(char));
    sc[0]=' ';
    sc[1]=',';
    sc[2]='\t';
    sc[3]='\0';
    while( *sc ) {
         if(*(sc++) == ch)
             return 1;
    }
    fmemfree(sc0);
    return 0;
}
 
int  words(const char* str) {
   int len = 0;
   char ch = 0;
   do {
       if( is_slash(*str) || *str == '\0') {
            if( ! is_slash(ch))
               ++len;
       }
       ch = *str;
   } while(*str++ != '\0');
   return  len;
}

char *word_find(char* str, int n) {
   int len = 0;
   char ch = 0;
   if (n==0) return str;
   do {
       if( is_slash(*str) || *str == '\0') {
            if( ! is_slash(ch))
               ++len;
	    if (len==n) return str;
       }
       ch = *str;
   } while(*str++ != '\0');
   return  str;
}

int main(int narg, char **argv)
{
  
  int col,i,k;
  int L;
  int strnum,num_clear;
  int nfields;
  float T,x;
  char *opts = "C:L:T:"; //flags
  char opt;
  char buf[STRLN];
  char **cols;
  char *ch;
  int C_flag=0, T_flag=0, L_flag=0;
  vector<float> data;
  vector<float> *data_averaged;
  FILE *fi;
  
  if (narg<4) 
    errquit("USE: godograph-clean -C<column number> -L<window width> [-T<threshold>] <godograph data file>\n");
  
  while((opt = getopt(narg, argv, opts)) != -1) { // вызываем getopt пока она не вернет -1 
  switch(opt) {
    case 'C': //
      C_flag=1;
      sscanf (optarg,"%d",&col);
      k=optind;
    break;
    
    case 'L': //
       L_flag=1;
       sscanf (optarg,"%d",&L);
       k = optind;
    break;
    
    case 'T': //
      T_flag=1;
      sscanf(optarg,"%g",&T);
      k = optind;
    break; 
  }
 }
 
  if ((!L_flag) || (!C_flag))
    errquit("USE: godograph-clean -C<column number> -L<window width> [-T<threshold>] <godograph data file>\n");
    
  fi=file_open(argv[k],"rt");
  
  fgets(buf,STRLN,fi); strnum=1;
  while (!feof(fi))
  {
    if (words(buf)<col+1)
    {
      fprintf(stderr,"Data format error at string %d of %s: %d fields found, %d expected\n"
                    ,strnum,argv[k],words(buf),col+1);
      exit(EXIT_FAILURE);
    }
    if (sscanf(word_find(buf,col),"%g",&x)!=1)
    {
      fprintf(stderr,"Data format error at string %d of %s: %d column is: "
                    ,strnum,argv[k],col);
      fputs(buf,stderr);
      exit(EXIT_FAILURE);
    }
    data.push_back(x);
    fgets(buf,STRLN,fi); strnum++;
  }
  strnum--;
  
  data_averaged=average(&data,strnum,L);
  
  num_clear=0;
  
  rewind(fi);
  
  for (i=0;i<strnum;i++)
  {
    if (!T_flag)  //averaging only
    {
      fgets(buf,STRLN,fi);
      ch=strchr(buf,'\n');
      if (ch!=NULL) (*ch)='\0';
      fputs(buf,stdout);
      fprintf(stdout," %g\n",(*data_averaged)[i]);
    }
    else
    {
      fgets(buf,STRLN,fi);
      if (fabs(data[i]-(*data_averaged)[i])<T)
      {
        fputs(buf,stdout);
	num_clear++;
      }
    }
  }
  
  if (T_flag) 
    fprintf(stderr,"%d of %d points are written (%g\%).\n",num_clear,strnum,(float)num_clear/strnum*100.0);
  
  delete data_averaged;
  fclose(fi);
}


