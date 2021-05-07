#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <values.h>

#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QString>


#include "gstdata.h"
#include "gsterror.h"
#include "gstdefs.h"

#include "gmt_nan.h"

typedef struct
  {
    int nx;
    int ny;
    int node_offset;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
    double x_inc;
    double y_inc;
    double z_scale_factor;
    double z_add_offset;
    char x_units[80];
    char y_units[80];
    char z_units[80];
    char title[80];
    char command[320];
    char remark[160];
  } gmt_h;



GSTgriddedData *readGMT(char *fileName)
{
  FILE *fi;
  gmt_h H;
  int i;
  GSTgriddedData *A;
  GSTgrid_type t;
  char str[GST_STRLN];


  if((fi=fopen(fileName,"rb"))==NULL)
  {
    sprintf(str,"Can't open file %s",fileName);
    GSTerror(str);
  }

  // fread(&H,sizeof(gmt_h),1,fi);
  /* Because gmt_h is not 64-bit aligned we must read it in parts */
  if (fread ((void *)&(H.nx), 3*sizeof (int), (size_t)1, fi) != 1 ||
      fread ((void *)&(H.x_min),sizeof(gmt_h)-((long)&(H.x_min)-(long)&(H.nx)),(size_t)1,fi) != 1)
    GSTerror("Error reading grid header in readGMT\n");

  if(H.node_offset==0)t=node;
  else t=area;

  fprintf(stderr,"At readGMT NX: %d, NY: %d, type: %d, X: %lg-%lg, Y: %lg-%lg\n",H.ny,H.nx,H.node_offset,H.y_min,H.y_max,H.x_min,H.x_max);
  fprintf(stderr,"At readGMT Z: %lg - %lg, x_inc: %lg, y_inc: %lg, z_scale_factor: %lg, z_add_offset: %lg\n",
                              H.z_min,H.z_max,H.x_inc,H.y_inc,H.z_scale_factor,H.z_add_offset);
  fprintf(stderr,"At readGMT Command: "); fputs(H.command,stderr); fprintf(stderr,"\n");
  fprintf(stderr,"At readGMT Title: "); fputs(H.title,stderr); fprintf(stderr,"\n");


  A = new GSTgriddedData(H.ny,H.nx,H.y_min,H.y_max,H.x_min,H.x_max,t);

  for (i=H.ny-1;i>=0;i--)
     if (fread(A->data()[i],sizeof(float),H.nx,fi)!=(unsigned)H.nx)
     {
      sprintf(str,"Read error while reding file %s",fileName);
      GSTerror(str);
     }

  fclose(fi);
 return A;
}


 int writeGMT (char *FileName, GSTgriddedData &A, double z_scale_factor=1, double z_add_offset=0,
  char *x_units="\0", char *y_units="\0",char *z_units="\0", char *title="\0" )


{
  FILE *fi;
  int i,j;
  GSTgrid_type t;
  gmt_h H;
  char str[GST_STRLN];
  float GMT_f_NaN;

  GMT_make_fnan(GMT_f_NaN);

  if((fi=fopen(FileName,"wb"))==NULL)
  {
    sprintf(str,"Can't open file %s",FileName);
    GSTerror(str);
  }

  A.size(H.ny,H.nx);
  A.getArea(H.y_min, H.y_max, H.x_min, H.x_max,t);
  if(t==node) H.node_offset=0;
  if(t==area) H.node_offset=1;

  A.getIntervals(H.y_inc, H.x_inc);


   for (i=0;(i<H.ny*H.nx) && (A.data()[0][i]>=GST_EMPTY_LIMIT);i++) ;
   H.z_min = MAXFLOAT;
   H.z_max= -MAXFLOAT;


     for (i=0;i<H.ny;i++)
     {
      for (j=0; j<H.nx; j++)
      {
       if (A.data()[i][j]<GST_EMPTY_LIMIT)
       {
         if(A.data()[i][j]> H.z_max) H.z_max=(double)A.data()[i][j];
         if(A.data()[i][j]< H.z_min) H.z_min=(double)A.data()[i][j];
       }
       else
         A.data()[i][j]=GMT_f_NaN;
      }
     }




  H.z_scale_factor=z_scale_factor;
  H.z_add_offset=z_add_offset;
  strcpy( H.x_units,x_units);
  strcpy( H.y_units,y_units);
  strcpy( H.z_units,z_units);
  strcpy( H.title,title);
  strcpy(H.command,"\0");
  strcpy(H.remark,"\0");


   //A.getComment(H.remark);


  //fwrite(&H,sizeof(gmt_h),1,fi);
  if (fwrite ((void*)&(H.nx),3*sizeof(int),(size_t)1,fi) != 1 ||
      fwrite ((void*)&(H.x_min),sizeof(gmt_h)-((long)&(H.x_min)-(long)&(H.nx)),(size_t)1,fi) != 1)
    GSTerror("Error writing grid header in writeGMT");

   for (i=H.ny-1; i>=0 ; i--)
     fwrite(A.data()[i],sizeof(float),H.nx ,fi) ;

  fclose(fi);
  return 1;
}

GSTgriddedData *readgrd(char *fileName, GSTgrid_type t = node)
{
  FILE *fi;
  char d[4];
  int count;
  int  Nr, Nc;
  int i,j;
  double X1, X2, Y1,Y2;
  GSTgriddedData *A;
  float Max;
  char str[GST_STRLN];

  if((fi=fopen(fileName,"rt"))==NULL)
  {
    sprintf(str,"Can't open file %s",fileName);
    GSTerror(str);
  }

  count=fscanf(fi,"%c%c%c%c",d+0,d+1,d+2,d+3);
  if(count!=4)
  {
    sprintf(str,"File %s format error",fileName);
    GSTerror(str);
  }

  if(d[0]!='D' || d[1]!='S' || d[2]!='A' || d[3]!='A')
  {
    sprintf(str,"File %s does not appear to be a Surfer ASCII grid",fileName);
    GSTerror(str);
  }

  count=fscanf(fi,"%d%d",&Nc,&Nr);
  if(count!=2)
  {
    sprintf(str,"File %s scanf error",fileName);
    GSTerror(str);
  }

  if(Nr<=0)
  {
    sprintf(str," in file %s Nr<=0 ",fileName);
    GSTerror(str);
  }

  if(Nc<=0)
  {
    sprintf(str," in file %s Nc<=0 ",fileName);
    GSTerror(str);
  }

  count=fscanf(fi,"%lg%lg%lg%lg",&X1,&X2,&Y1,&Y2);
  if(count!=4)
  {
    sprintf(str,"File %s scanf error",fileName);
    GSTerror(str);
  }
  printf("\nx1=%lg x2=%lg y1=%lg y2=%lg Nc=%d Nr=%d\n",X1,X2,Y1,Y2,Nc,Nr);

  //if(X2<=X1)
  //{
   // sprintf(str," in file %s X2<=X1 ",fileName);
   // GSTerror(str);
  //}

  if(Y2<=Y1)
  {
    sprintf(str," in file %s Y2<=Y1 ",fileName);
    GSTerror(str);
  }

  count=fscanf(fi,"%g%g",&Max,&Max);
  if(count!=2)
  {
    sprintf(str,"File %s scanf error",fileName);
    GSTerror(str);
  }

  A = new GSTgriddedData(Nr,Nc,X1,X2,Y1,Y2,t);
  for (i=0;i<Nr;i++)
  {
   for (j=0; j<Nc; j++)
    {
    count=fscanf(fi,"%g",A->data()[i]+j);
    if(count!=1)
    {
    sprintf(str,"File %s scanf error",fileName);
    GSTerror(str);
     }

    }
  }
  fclose(fi);
 return A;
}

void writegrd(char *FileName, GSTgriddedData &A)
{
  FILE *fi;
  int  Nr, Nc;
  int i,j;
  double X1, X2, Y1,Y2;
  GSTgrid_type t;
  char str[GST_STRLN];

  float Max,Min;

  if((fi=fopen(FileName,"wt"))==NULL)
  {
    sprintf(str,"Can't open file %s",FileName);
    GSTerror(str);
  }

  fprintf(fi,"DSAA\n");
  A.size(Nr,Nc);
  A.getArea(X1,X2,Y1,Y2,t);


  for (i=0;(i<Nr*Nc) && (A.data()[0][i]>=GST_EMPTY_LIMIT);i++) ;
  Min =(double)A.data()[0][i];
  Max=(double)A.data()[0][i];

  for (i=0;i<Nr;i++)
  {
   for (j=0; j<Nc; j++)
   {
    if (A.data()[i][j]<GST_EMPTY_LIMIT)
    {
      if(A.value(i,j)>Max)Max=A.value(i,j);
      if(A.value(i,j)<Min) Min=A.value(i,j);
    }
   }
  }
  fprintf(fi,"%d %d\n",Nc,Nr);
  fprintf(fi,"%g %g\n %g %g\n", X1,X2,Y1,Y2);
  fprintf(fi,"%g %g\n",Min, Max);

   for (i=0;i<Nr;i++)
  {
   for (j=0; j<Nc; j++)
    {
      fprintf(fi,"%g ",A.value(i,j));
    }
    fprintf(fi,"\n");
  }

  fclose(fi);
}

void writeXYZ(char *FileName, GSTgriddedData &A)
{
  FILE *fi;
  int  Nr, Nc;
  int i,j;
  double X1, X2, Y1,Y2;
  double dx, dy;
  GSTgrid_type t;
  char str[GST_STRLN];

  float X,Y;

  if((fi=fopen(FileName,"wt"))==NULL)
  {
    sprintf(str,"Can't open file %s",FileName);
    GSTerror(str);
  }

  A.size(Nr,Nc);
  A.getArea(X1,X2,Y1,Y2,t);
  A.getIntervals(dx,dy);

 if(t==node)
  {
   for (i=0;i<Nr;i++)
   {
    for (j=0; j<Nc; j++)
     {
      X=X1+dx*i;
      Y=Y1+dy*j;
     fprintf(fi,"%10g %10g %10g \n",X,Y,A.value(i,j));

     }
   }
  }

 if(t==area)
  {
   for (i=0;i<Nr;i++)
   {
    for (j=0; j<Nc; j++)
     {
      X=X1+dx*0.5+dx*i;
      Y=Y1+dy*0.5+dy*j;
      fprintf(fi,"%10g %10g %10g \n",X,Y,A.value(i,j));
     }
   }
  }

  fclose(fi);
}

GSTgriddedData *readESRIasc(char *fileName)
{
  int nr, nc, i, j;
  int npar;
  char nodata[GST_STRLN];
  double xllc, yllc, clsz;
  char s1[GST_STRLN];
  FILE *fi;

  if((fi=fopen(fileName,"rt"))==NULL)
  {
    sprintf(s1,"Can't open file %s",fileName);
    GSTerror(s1);
  }

  npar=0;
  do
  {
    i=fscanf(fi,"%s",s1);
    if (i==1)
    {
      if (strcmp(s1,"nrows")==0)
      {
        if (fscanf(fi,"%d",&nr)!=1) GSTerror("Data format error at nrows");
	npar++;
      }
      if (strcmp(s1,"ncols")==0)
      {
        if (fscanf(fi,"%d",&nc)!=1) GSTerror("Data format error at ncols");
	npar++;
      }
      if (strcmp(s1,"xllcorner")==0)
      {
        if (fscanf(fi,"%lg",&xllc)!=1) GSTerror("Data format error at xllcorner");
	npar++;
      }
      if (strcmp(s1,"yllcorner")==0)
      {
        if (fscanf(fi,"%lg",&yllc)!=1) GSTerror("Data format error at yllcorner");
	npar++;
      }
      if (strcmp(s1,"cellsize")==0)
      {
        if (fscanf(fi,"%lg",&clsz)!=1) GSTerror("Data format error at cellsize");
	npar++;
      }
      if (strcmp(s1,"NODATA_value")==0)
      {
        if (fscanf(fi,"%s",nodata)!=1) GSTerror("Data format error at NODATA_value");
	npar++;
      }
    }
    else
      if (!feof(fi)) GSTerror("Data format error");
  }
  while ((npar<6) && (!feof(fi)));

  if (npar<6) GSTerror("Insufficient parameters in ESRI ascii raster data file");

  GSTgriddedData *D = new GSTgriddedData(nr,nc,xllc,xllc+nr*clsz,yllc,yllc+nc*clsz,node);

  for (i=0;i<nr;i++)
    for (j=0;j<nc;j++)
    {
      if (fscanf(fi,"%s",s1)!=1) GSTerror("Data format error at data reading");
      if (strcmp(s1,nodata)==0) D->data()[i][j]=GST_EMPTY_VALUE;
      else if (sscanf(s1,"%g",D->data()[i]+j)!=1) GSTerror("Data format error at data reading");
    }
  fclose(fi);
  return D;
}


float GSTtoFloat(QString &s, bool *ok=0)
{
  bool f = true;

  s.toFloat(&f);
  if(f)
   {
    if (ok!=0) *ok=true;
    return s.toFloat();
   }
  else
  {
   if (s.toLower()=="nan")
    {
     if (ok!=0) *ok=true;
     return GST_EMPTY_VALUE;
     }
   else
   {
    if (ok!=0) *ok=false;
    return 0.0;
    }
  }
}

int CountColumns(char *fileName, char *divid, QStringList &names, int &headline)
{
//если первая строка коментарий- сравнеие форматов происходит у второй и третьей строки
/*функция прденозначена для анализа колоночных данных.
  возвращает число колонок в файле fileName,
  в переменную *divid записывается "разделитель" с помощью которого разделены значения в файле
  возможные разделители " " "," ";" , headline-количество строк вначаеле файла которые, не нужно
  анализировать, для анализа можно оставлять одну строку заголовка ; желательно, чтобы она
  состояла из того же  числа колонок, что и сам файл а колонки были разделены тем же
  разделителем, что и сам файл, если заголовок найден, то он будет записан в переменную &names)
  */

 bool ok;
 int i;
 int sz;
 char str[GST_STRLN];
 QString line1, line2;
 QStringList cols1, cols2;
 QFile file(fileName);


  if (!file.open(QIODevice::ReadOnly))
  {
      sprintf(str,"Can't open file %s",fileName);
      GSTerror(str);
  }

    QTextStream in(&file);
    for(i=0;i<headline;i++) in.readLine();

    line1 = in.readLine();

    line1 = line1.simplified();

    cols1 = line1.split(",");
    sz=cols1.size();


    strcpy(divid,",");

    if (cols1.size() == 1)
    {
      cols1 = line1.split(" ");
      strcpy(divid," ");
      if (cols1.size() == 1)
      {
        cols1 = line1.split(";");
	strcpy(divid,";");
      }
    }


     sz=cols1.size();
     ok=true; i=0;
     do
     {
       GSTtoFloat(cols1[i],&ok);
       i++;
     }
     while ((ok) && (i<sz));

     if(!ok)
     {
        headline=headline+1;
        names=cols1;
        GSTwarning("Header line found \n");

	line1=(in.readLine()).simplified();
        cols1 = line1.split(",");

	strcpy(divid,",");
        if (cols1.size() == 1)
        {
          cols1 = line1.split(" ");
          strcpy(divid," ");
          if (cols1.size() == 1)
          {
           cols1 = line1.split(";");
	   strcpy(divid,";");

         }
        }
	sz=cols1.size();
	ok=true; i=0;
          do
            {
	     GSTtoFloat(cols1[i],&ok);
             i++;
            }
           while ((ok) && (i<sz));


          if(!ok)
           {
             sprintf(str,"format of the second line in file %s not float\n",fileName);
             GSTerror(str);
           }

      }

   line2 = in.readLine();
   line2 = line2.simplified();
   cols2 = line2.split(divid);

   if(cols1.size()!= cols2.size())
   {
      sprintf(str,"different format of lines in file %s\n",fileName);
      GSTerror(str);
   }

   if (sz==1)divid[0]='\0';

    ok=true; i=0;
    do
     {
     GSTtoFloat(cols2[i],&ok);
     i++;
      }
     while ((ok) && (i<sz));


     if(!ok)
      {
       sprintf(str,"not float format in file %s\n",fileName);
       GSTerror(str);
      }

  file.close();
  return sz;
}

GSTscatterData *ReadColumns
(const char *fileName, int M, int N, int *n, char *name,const char *divid, int headline=0, bool strict = false)
{
 /*функция предназначена для чтения файла состоящего из колонок, в GSTscatterData

  M-читсло колонок в файле N- число колонок которые необходимо прочитать
  *n массив размерностью N,с номерами колонок, которые необходимао заприсать,
  *name- массив размерностью N с именами тороые будут даны записываемым колонкам,
  *divid- разделитель используемый в файле, headline- число срок отведенное под комментрарий,
  т.е. число строк которые необходимо пропустить
  strict - устанавливает требование к точности соответствия числа читаемых и преобразуемых к числовому виду
  полей заданному N. Если strict=true, то несоответствие приводит к ошибке, если strict=false,
  то прочитываются и записываются первые n<=N читаемых полей, прочие игнорируются, недостающие поля
  структуры данных заполняются пустыми значениями.
   */

  QFile file(fileName);
  bool ok;
  QString line;
  QStringList cols;
  int actM;
  int c;
  int i;
  GSTscatterData *D;
  char str[GST_STRLN];
  GSTpointMap A;

  if(N>M)GSTerror("N>M");
  D=new GSTscatterData(name);

  if (!(file.open(QIODevice::ReadOnly)))
  {
   sprintf(str,"Can't open file %s",fileName);
   GSTerror(str);
  }

  QTextStream in(&file);
  for(i=0; i<headline;i++) in.readLine();

  c=0;

  for (i=0;i<N;i++)
  {
    if ((n[i]<0) || (n[i]>M))
      GSTerror("Column number out of range in function ReadColumns\n");
  }

  while (!in.atEnd())
  {
   line=in.readLine();
   line = line.simplified();
   cols = line.split(divid);
   c=c+1;

   actM=cols.size();
   if (actM!=M && strict)
   {
     sprintf(str,"number of columns in line %d, file %s file != M and strict rules are set",c,fileName);
     GSTerror(str);
   }

    for(i=0;i<N;i++)
    {
     ok=true;

     if (n[i] < actM)
       A[name[i]]=GSTtoFloat(cols[n[i]],&ok);
     else
       A[name[i]]=GST_EMPTY_VALUE;

     if(!ok)
     {
       if (strict) {sprintf(str,"Data format error in file %s and strict rules are set",fileName);
       GSTerror(str);}
       else A[name[i]]=GST_EMPTY_VALUE;
     }

    }
   D->pntAdd(A);
  }
 file.close();
 return D;
}


void  WriteColumns
(const char *fileName, GSTscatterData &D,char *names=NULL,const char *divid=",", int fieldwidth=0)
{

  int nf;
  int nP=D.numPnt();
  int i,j;
  GSTpointMap A;
  char str[GST_STRLN];
  std::string s;

  if(names==NULL)
    D.namesFLD(s);
  else
  {
    s.erase();
    s+=names;
  }
  nf=s.size();
  if(nf>D.numFld()) GSTwarning("Too more names of filds \n");

  QFile file(fileName);
  if(!file.open(QIODevice::WriteOnly))
  {
      sprintf(str,"Can't open file %s\n",fileName);
      GSTerror(str);
  }

 QTextStream out(&file);

 for (i=0;i<nf-1;i++)
   out << qSetFieldWidth(fieldwidth) << right << s[i] << divid;
 out << qSetFieldWidth(fieldwidth) << right << s[nf-1] << endl;


 for(i=0; i<nP; i++)
  {

    for (j=0;j<nf-1;j++)
      out << qSetFieldWidth(fieldwidth) << right << D.f(i,s[j]) << divid;
    out << qSetFieldWidth(fieldwidth) << right << D.f(i,s[nf-1]) << endl;
  }


 file.close();

}









