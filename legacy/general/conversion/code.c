#include <stdio.h>

int main(int narg, char **argv)
{
  int i;
  unsigned char ch;
  FILE *fi, *fo;
  
  fo=fopen("table.txt","wt");
  for(i=0;i<256;i++)
  {
    fprintf(fo,"%5d%c\n",i,(unsigned char)i);
  }
  fclose(fo);
  return 0;
}
