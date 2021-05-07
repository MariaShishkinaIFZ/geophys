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
  
/*��������� �� ���������: z_scale_factor- �����,�� ������� ����� �������� ��� �������� Z; z_add_offset-�����,������� ����� ��������� �� ���� ��������� Z; x_units- ����������� x; y_units- ����������� �;z_units-����������� z;  title- ���������*/     
 

GSTgriddedData *readgrd(char *fileName, GSTgrid_type t = node);
 
int writegrd(char *FileName, GSTgriddedData &A); 
  
int writeXYZ(char *FileName, GSTgriddedData &A);
 
GSTgriddedData *readESRIasc(char *fileName);


void  WriteColumns
(const char *fileName, GSTscatterData &D,char *names=NULL,const char *divid=",", int fieldwidth=0);

/*����������  � "���������� ����" fileName, ������ �� GSTscatterData D.
� ���������� names ���������� ����� �����, ������� ���������� ��������. �� ��������� ����� 
�������� ��� ���� GSTscatterData. ���� ����������� ����, �������� ��� � GSTscatterData ��
����� �������� ������ ��������.  
 ������ ������ �����
 ����������- �������� �������� (����� GSTscatterData &D) divid-����������� ������� �����
 ��������� �������� � �����, fieldwidth- ������ ���� ���������� ��� ������� ��������*/
 
 float GSTtoFloat(QString &s, bool *ok=0);
 /* ��������� QString &s � float ������ ���� ����������� �������� "NAN" ���������� 
  GST_EMPTY_VALUE, ���� bool *ok �������, �� � ������ ��������� �������������� ������������ 
  true, � � ������ �� ����� false  */
  
  int CountColumns(char *fileName, char *divid, QStringList &names, int &headline);
/*������� ������������� ��� ������� ���������� ������.
  ���������� ����� ������� � ����� fileName,
  � ���������� *divid ������������ ��������� ���������� "�����������" � ������� ��������
  ��������� �������� � �����   ��������� ����������� " " "," ";" 
  (����� ���������� �������� � ������ �� ������ ������ �����������),
  headline-���������� ����� �������� ����� �������, �� �����  �������������, ��� ������� �����
  ��������� ���� ������ ��������� ; ����������, ����� ���  �������� �� ���� ��  ����� �������,
  ��� � ��� ���� � ������� ���� ��������� ��� ��  ������������, ��� � ��� ����, ���� ���������
  ������, �� �� ����� ������� � ���������� &names, �  headline=headline+1
  ������� ��������� ����������� ������� � 1 � 2 ������ ������ )
  */
  
  GSTscatterData *ReadColumns
(const char *fileName, int M, int N, int *n, char *name, const char *divid, int headline=0);
 /*������� ������������� ��� ������ ����� ���������� �� �������, � GSTscatterData
 
  M-������ ������� � ����� N- ����� ������� ������� ���������� ��������, 
  *n ������ ������������ N,� �������� �������(��������� ������� ���������� � 0),
   ������� ���������� ��������,
  *name- ������ ������������ N � ������� ������ ����� ���� ������������ ��������,
  *divid- ����������� ������������ � �����, headline- ����� ���� ���������� ��� ������������,
   �.�. ����� ����� ������� ���������� ����������
   */
