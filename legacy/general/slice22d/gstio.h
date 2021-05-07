#include <stdio.h>
#include <math.h>

#include <QFile>
#include <QTextStream>
#include <QStringList>
#include <QString>


#include "gstdata.h"
#include "gsterror.h"
#include "gstdefs.h"


GSTgriddedData *readGMT(char *fileName);

 int writeGMT(char *FileName, GSTgriddedData &A, double z_scale_factor=1, double z_add_offset=0,
  char *x_units="\0", char *y_units="\0",char *z_units="\0", char *title="\0" );
  
/*параметры по умолчанию: z_scale_factor- число,на которое будут умноежны все значения Z; z_add_offset-число,которое будет добавлено ко всем значениям Z; x_units- размерность x; y_units- размерность у;z_units-размерность z;  title- заголовок*/     
 

GSTgriddedData *readgrd(char *fileName, GSTgrid_type t = node);
 
int writegrd(char *FileName, GSTgriddedData &A); 
  
int writeXYZ(char *FileName, GSTgriddedData &A);
 
GSTgriddedData *readESRIasc(char *fileName);


void  WriteColumns
(const char *fileName, GSTscatterData &D,char *names=NULL,const char *divid=",", int fieldwidth=0);

/*записывает  в "колоночный файл" fileName, данные из GSTscatterData D.
В переменной names передаются имена полей, которые необходимо записать. По умолчанию будут 
записаны все поля GSTscatterData. если встретиться поле, которого нет в GSTscatterData то
будут записаны пустые значения.  
 первая строка файла
 коментарий- название столбцов (полей GSTscatterData &D) divid-разделитель которым будут
 разделены значения в файле, fieldwidth- ширина поля отведенная для каждого значения*/
 
 float GSTtoFloat(QString &s, bool *ok=0);
 /* переводит QString &s в float формат если встречается значение "NAN" возвращает 
  GST_EMPTY_VALUE, если bool *ok передан, то в случае успешного преобразования возвращается 
  true, а в случае не удачи false  */
  
  int CountColumns(char *fileName, char *divid, QStringList &names, int &headline);
/*функция прденозначена для анализа колоночных данных.
  возвращает число колонок в файле fileName,
  в переменную *divid возвращается найденный программой "разделитель" с помощью которого
  разделены значения в файле   возможные разделители " " "," ";" 
  (после последнего значения в строке не должен стоять разделитель),
  headline-количество строк вначаеле файла которые, не нужно  анализировать, для анализа можно
  оставлять одну строку заголовка ; желательно, чтобы она  состояла из того же  числа колонок,
  что и сам файл а колонки были разделены тем же  разделителем, что и сам файл, если заголовок
  найден, то он будет записан в переменную &names, а  headline=headline+1
  функция проверяет соответсвие формату в 1 и 2 строки данных )
  */
  
  GSTscatterData *ReadColumns
(const char *fileName, int M, int N, int *n, char *name, const char *divid, int headline=0);
 /*функция предназначена для чтения файла состоящего из колонок, в GSTscatterData
 
  M-читсло колонок в файле N- число колонок которые необходимо записать, 
  *n массив размерностью N,с номерами колонок(нумерация колонок начинается с 0),
   которые необходимо записать,
  *name- массив размерностью N с именами тороые будут даны записываемым колонкам,
  *divid- разделитель используемый в файле, headline- число срок отведенных под комментрарий,
   т.е. число строк которые необходимо пропустить
   */
