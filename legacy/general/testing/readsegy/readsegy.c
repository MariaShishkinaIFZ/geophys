#include<stdio.h>

const int N = 480;  //количество трасс
 
struct trace{
    int num;     
    int Xtrans;  // пункт взрыва
    int Xrec;    //пункт приема
    double OGT;  //огт
};
 
int convert4(int value)
{
    int res;
    
    res = value << 24;
    res |= (value <<  8) & 0xFF0000;
    res |= (value >>  8) & 0xFF00;
    res |= (value >> 24) & 0xFF;
    
    return res;
}
 
short convert2(short value)
{
  short res;
  res = value << 8;
  value = (value >> 8) & 0xFF;
  res |= value;
 
  return res;
}
 
int main()
{    
     struct trace Traces[N];
     int temp[1500];
     unsigned long i, buf, j;   
     short k, m, metka;          
     long L;                     
     double buf2, current_ogt;
     unsigned char c;
 
     FILE *fp = fopen("Sample_i4.sgy", "rb");
     FILE *fp2 = fopen("result.sgy", "wb+");
     
     fseek(fp, 3220, 0);
     fread(&k, 4, 1, fp);
     k = convert2(k);
     L = k*4;
     
     fseek(fp, 3224, 0);
     fread(&m, 4, 1, fp);
     m = convert2(m);
     
     
     for(i = 0; i < N; i++)
     {
        (Traces[i]).num = i+1;
        
        fseek(fp, 3672 + i*(L + 240), 0);    
        fread(&(Traces[i].Xtrans), 4, 1, fp);
        (Traces[i]).Xtrans = convert4((Traces[i]).Xtrans);
        
        fseek(fp, 3680 + i*(L + 240), 0);
        fread(&(Traces[i].Xrec),   4, 1, fp);
        (Traces[i]).Xrec = convert4((Traces[i]).Xrec);
        
        Traces[i].OGT = 1.0*(Traces[i].Xtrans + Traces[i].Xrec)/2;
    }
 
//------------------------------------------------------------------------------
 
    do
    {
      metka = 0;
      
      for(i = 0; i < N-1; i++)
            if(Traces[i].OGT > Traces[i+1].OGT)
            {
                     buf = Traces[i].num;
                     Traces[i].num = Traces[i+1].num;
                     Traces[i+1].num = buf;
                     buf = Traces[i].Xtrans;
                     Traces[i].Xtrans = Traces[i+1].Xtrans;
                     Traces[i+1].Xtrans = buf;
                     buf = Traces[i].Xrec;
                     Traces[i].Xrec = Traces[i+1].Xrec;
                     Traces[i+1].Xrec = buf;
                     buf2 = Traces[i].OGT;
                     Traces[i].OGT = Traces[i+1].OGT;
                     Traces[i+1].OGT = buf2;
                     
                     metka = 1;
            }
      }
    while(metka);
    
    
    
    
//------------------------------------------------------------------------------
    for(i = 0; i < N; i++)
    {
        printf("N = %d\n", (Traces[i]).num);
        printf("Xtrans = %d\n", (Traces[i]).Xtrans);
        printf("Xres = %d\n", (Traces[i]).Xrec);
        printf("OGT = %f\n", (Traces[i]).OGT);
        printf("\n");
    }
    
    printf("%d\n", L);
    printf("%d\n", m);
 
//------------------------------------------------------------------------------    
fseek(fp, 0, 0);  
 
for(i = 0; i < 3600; i++)
{
      fread(&c, 1, 1, fp);
      fwrite(&c,1,1,fp2);
}
 
//----------------------------------------копирование заголовка
 
for(j = 0; j < 1500; j++)    
      temp[j] = 0;      // î÷èñòêà òåìï ìàññèâà
 
current_ogt = Traces[0].OGT;
 
i = 0;
 
while ( i < N)
{
  if(Traces[i].OGT == current_ogt)
      {
          for(j = 0; j < 1500; j++)
          {
                fseek(fp, 3600 + (Traces[i].num-1)*(6000+240) + 240 + 4*j, 0);
                fread(&buf, 4, 1, fp);
                temp[j] += buf;
          }
      i++;
      continue;
      }  
    
           
      fseek(fp, 3600 + (Traces[i-1].num-1)*(6000+240), 0);      
           
      for(j = 0; j < 240; j++)
    {
        fread(&c, 1, 1, fp);                         
      fwrite(&c, 1, 1,fp2);
      }
      
      for(j = 0; j < 1500; j++)                        
      {
          fwrite(&(temp[j]), 4, 1,fp2);
      }
           
      for(j = 0; j < 1500; j++)           
          temp[j] = 0;
                 
      current_ogt = Traces[i].OGT;
}
 
fseek(fp, 3600 + (Traces[N-1].num-1)*(6000+240), 0);      
           
for(j = 0; j < 240; j++)
{
  fread(&c, 1, 1, fp);                         
  fwrite(&c, 1, 1,fp2);
}
      
for(j = 0; j < 1500; j++)                        
{
  fwrite(&(temp[j]), 4, 1,fp2);
}
 
fclose (fp);
fclose (fp2);
return 0;
}