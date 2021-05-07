#include <stdio.h>

int main(int narg, char **argv)
{
  unsigned char ch;
  FILE *fi, *fo;
  
  fi=fopen(argv[1],"rt");
  fo=fopen(argv[2],"wt");
  ch=fgetc(fi);
  while (!feof(fi))
  {
    if (((int)ch!=10) && (ch!='#')) fputc(ch,fo);
    else if (ch=='#') fputc((char)10,fo);
    ch=fgetc(fi);
  }
  fclose(fo);
  fclose(fi); 
  return 1;
}
