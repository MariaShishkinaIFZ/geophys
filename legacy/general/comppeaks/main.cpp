#include <QtCore>

#include "gstio.h"
#include "gsterror.h"
#include "gstdata.h"

int main(int argc, char* argv[])
{
   GSTscatterData *dataP, *compP;
   double L2, xi2;
   int n[10];
   int i,j,num_coinc;
   char field_names[10];
   char s[GST_STRLN];

   if (argc<3)
   {
     sprintf(s,"comppeacks <6-cols data file [.tpk]> <4-cols computed peaks file [.tp]>\n");
     GSTerror(s);
   }

   n[0]=0; n[1]=1; n[2]=2; n[3]=3;
   field_names[0]='n'; field_names[1]='d'; field_names[2]='a'; field_names[3]='t';
   compP=ReadColumns(argv[2],4,4,n,field_names," ",0,false);
   n[0]=1; n[1]=4; n[2]=5;
   field_names[0]='n'; field_names[1]='t'; field_names[2]='e';
   dataP=ReadColumns(argv[1],6,3,n,field_names," ",0,false);

   L2=xi2=0.0;
   num_coinc=0;
   for (i=0;i<compP->numPnt();i++)
   {
     j=0;
     do
     {
       if (compP->f(i,'n')==dataP->f(j,'n'))
       {
         fprintf(stdout,"%5d %13g %13g %13g %13g %13g %13g\n",(int)compP->f(i,'n'),compP->f(i,'d'),compP->f(i,'a'),
                         compP->f(i,'t'),dataP->f(j,'t'),compP->f(i,'t')-dataP->f(j,'t'),dataP->f(j,'e'));
         L2+=pow(compP->f(i,'t')-dataP->f(j,'t'),2.0);
         xi2+=pow((compP->f(i,'t')-dataP->f(j,'t'))/dataP->f(j,'e'),2.0);
         num_coinc++;
         break;
       }
       j++;
     }
     while (j<dataP->numPnt());
   }
   fprintf(stderr,"N=%d L2=%13g xi2=%13g\n",num_coinc,pow(L2,0.5),pow(xi2/num_coinc,0.5));
}
