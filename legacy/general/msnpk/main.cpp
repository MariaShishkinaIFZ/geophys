#include <QtCore>

#include "gstio.h"
#include "gsterror.h"
#include "gstdata.h"

float gasdev(long *idum);

int main(int argc, char* argv[])
{
   GSTscatterData *dataP, *compP;
   int n[10];
   int i,j,num_coinc;
   long idum=(-13);
   char field_names[10];
   char s[GST_STRLN];

   if (argc<3)
   {
     fprintf(stderr,"msnpk <6-cols data file [.tpk]> <8-cols computed peaks file from getpeaks [.dat]>\n");
     sprintf(s,"Selects computed time corresponding fod data peaks and adds noise according to given rms error\n");
     GSTerror(s);
   }

   n[0]=0; n[1]=1; n[2]=2; n[3]=3; n[4]=4; n[5]=5;
   field_names[0]='p'; field_names[1]='n'; field_names[2]='d'; field_names[3]='a'; field_names[4]='t'; field_names[5]='e';
   dataP=ReadColumns(argv[1],6,6,n,field_names," ",0,false);

   n[0]=0; n[1]=7;
   field_names[0]='n'; field_names[1]='t';
   compP=ReadColumns(argv[2],8,2,n,field_names," ",0,false);

   idum=-rand();

   num_coinc=0;
   for (i=0;i<compP->numPnt();i++)
   {
     j=0;
     do
     {
       if (compP->f(i,'n')==dataP->f(j,'n'))
       {
         fprintf(stdout,"%5d %5d %13g %13g %13g %13g\n",(int)dataP->f(j,'p'),(int)dataP->f(j,'n'),dataP->f(j,'d'),dataP->f(j,'a'),
                         compP->f(i,'t')+dataP->f(j,'e')*gasdev(&idum),dataP->f(j,'e'));
         num_coinc++;
         break;
       }
       j++;
     }
     while (j<dataP->numPnt());
   }
   fprintf(stderr,"N data points = %d Found = %d\n",(int)(dataP->numPnt()),num_coinc);
}
