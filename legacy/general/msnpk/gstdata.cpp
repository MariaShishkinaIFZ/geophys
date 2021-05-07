/*Версия от 29.01.06, не полностью тестированная*/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gstdata.h"
#include "gsterror.h"
#include "gstdefs.h"
#include "gstmath.h"

GSTgriddedData::GSTgriddedData
(int N_row, int N_col, double X_min, double X_max, double Y_min, double Y_max, GSTgrid_type type)
{
  if ((N_row<=0) || (N_col<=0)) GSTerror("Error constructing GSTgriddedData object");
  Nr=N_row;
  Nc=N_col;
  
  Xmin=X_min;
  Xmax=X_max;
  Ymin=Y_min;
  Ymax=Y_max;
  
  grid_type=type;
  
  comment[0]='\0';
  
  switch (type)
  {
    case (node): 
      if (Nr!=1)
        dX=(Xmax-Xmin)/(Nr-1);
      else
        dX=GST_EMPTY_VALUE;
      if (Nc!=1)
	dY=(Ymax-Ymin)/(Nc-1);
      else
        dY=GST_EMPTY_VALUE;
      break;
    case (area):
      dX=(Xmax-Xmin)/Nr;
      dY=(Ymax-Ymin)/Nc;
  }
  V = new float*[Nr];
  if (V==NULL) GSTerror("Insufficient memory");
  V[0] = new float[Nr*Nc];
  if (V[0]==NULL) GSTerror("Insufficient memory");
  for (int i=1;i<Nr;i++) V[i]=V[i-1]+Nc;
}

GSTgriddedData::GSTgriddedData(const GSTgriddedData &rhs)
{
  int N_row, N_col;
  double x0, x1, y0, y1;
  GSTgrid_type t;
  
  rhs.size(N_row,N_col);
  rhs.getArea(x0,x1,y0,y1,t);
  
  Nr=N_row;
  Nc=N_col;
  
  Xmin=x0;
  Xmax=x1;
  Ymin=y0;
  Ymax=y1;
  
  grid_type=t;
  
  comment[0]='\0';
  
  switch (t)
  {
    case (node): 
      if (Nr!=1)
        dX=(Xmax-Xmin)/(Nr-1);
      else
        dX=GST_EMPTY_VALUE;
      if (Nc!=1)
	dY=(Ymax-Ymin)/(Nc-1);
      else
        dY=GST_EMPTY_VALUE;
      break;
    case (area):
      dX=(Xmax-Xmin)/Nr;
      dY=(Ymax-Ymin)/Nc;
  }
  V = new float*[Nr];
  if (V==NULL) GSTerror("Insufficient memory");
  V[0] = new float[Nr*Nc];
  if (V[0]==NULL) GSTerror("Insufficient memory");
  for (int i=1;i<Nr;i++) V[i]=V[i-1]+Nc;
  
  float **rhs_data=rhs.data();
  for (int i=0;i<Nr;i++)
    for (int j=0;j<Nc;j++)
      V[i][j]=rhs_data[i][j];
}

GSTgriddedData::~GSTgriddedData()
{
  delete[] V[0];
  delete[] V;
  //printf("GSTgriddedData destructor called\n");
}
void GSTgriddedData::size(int &N_row, int &N_col) const
{
  N_row=Nr; N_col=Nc;
}

void GSTgriddedData::setArea(double X_min, double X_max, double Y_min, double Y_max, GSTgrid_type type)
{
  Xmin=X_min;
  Xmax=X_max;
  Ymin=Y_min;
  Ymax=Y_max;
  grid_type=type;
  switch (type)
  {
    case (node): 
      if (Nr!=1)
        dX=(Xmax-Xmin)/(Nr-1);
      else
        dX=GST_EMPTY_VALUE;
      if (Nc!=1)
	dY=(Ymax-Ymin)/(Nc-1);
      else
        dY=GST_EMPTY_VALUE;
      break;
    case (area):
      dX=(Xmax-Xmin)/Nr;
      dY=(Ymax-Ymin)/Nc;
  }
}

int GSTgriddedData::getArea(double &X_min, double &X_max, double &Y_min, double &Y_max, GSTgrid_type &type) const
{
  X_min=Xmin;
  X_max=Xmax;
  Y_min=Ymin;
  Y_max=Ymax;
  type=grid_type;
  
  if (Xmin<GST_EMPTY_LIMIT) return 1;
  else return 0;
}

void GSTgriddedData::setType(GSTgrid_type type)
{
  if (type!=grid_type)
  grid_type=type;
  switch (grid_type)
  {
    case (node):
      if (Nr!=1)
        dX=(Xmax-Xmin)/(Nr-1);
      else
        dX=GST_EMPTY_VALUE;
      if (Nc!=1)
	dY=(Ymax-Ymin)/(Nc-1);
      else
        dY=GST_EMPTY_VALUE;
      break;
    case (area):
      dX=(Xmax-Xmin)/Nr;
      dY=(Ymax-Ymin)/Nc;
  }
}

int GSTgriddedData::getIntervals(double &d_X, double &d_Y) const
{
  d_X=dX;
  d_Y=dY;
  
  if (Xmin<GST_EMPTY_LIMIT) return 1;
  else return 0;
}

float GSTgriddedData::value(int r, int c) const
{
  if ((r>=0) && (r<Nr) && (c>=0) && (c<Nc))
    return V[r][c];
  else
    return GST_EMPTY_VALUE;
}

float GSTgriddedData::value(double x, double y) const
{
  int r,c;
  switch (grid_type)
  {
    case node:
      r = lround((x-Xmin)/dX);
      c = lround((y-Ymin)/dY);
    break;
    case area:
      r = (int)floor((x-Xmin)/dX);
      c = (int)floor((y-Ymin)/dY);
  }
  if ((r>=0) && (r<Nr) && (c>=0) && (c<Nc))
    return V[r][c];
  else
    return GST_EMPTY_VALUE;
}

void GSTgriddedData::setValue(float v)
{
  for (int i=0;i<Nr;i++)
    for (int j=0;j<Nc;j++)
      V[i][j]=v;
}

void GSTgriddedData::setValue(int r, int c, float v)
{
  if ((r>=0) && (r<Nr) && (c>=0) && (c<Nc))
    V[r][c]=v;
}

void GSTgriddedData::setValue(double x, double y, float v)
{
  int r,c;
  switch (grid_type)
  {
    case node:
      r = lround((x-Xmin)/dX);
      c = lround((y-Ymin)/dY);
    break;
    case area:
      r = (int)floor((x-Xmin)/dX);
      c = (int)floor((y-Ymin)/dY);
  }
  if ((r>=0) && (r<Nr) && (c>=0) && (c<Nc))
    V[r][c]=v;
}

char *GSTgriddedData::getComment(char *comm) const
{
  strcpy(comm,comment);
  return comm;
}

void GSTgriddedData::setComment(const char *comm)
{
  strncpy(comment,comm,GST_STRLN);
}

bool GSTgriddedData::haveEqualArea(GSTgriddedData &d) const
{
  int dNr, dNc;
  double x0, x1, y0, y1;
  GSTgrid_type t;
  d.size(dNr,dNc);
  d.getArea(x0,x1,y0,y1,t);
  
  return (Nr==dNr) && (Nc==dNc) && GSTdoublesEqual(Xmin,x0) &&
         GSTdoublesEqual(Xmax,x1) && GSTdoublesEqual(Ymin,y0) &&
	 GSTdoublesEqual(Ymax,y1);
}

bool GSTgriddedData::haveEqualGrid(GSTgriddedData &d) const
{
  int dNr, dNc;
  double x0, x1, y0, y1;
  GSTgrid_type t;
  d.size(dNr,dNc);
  d.getArea(x0,x1,y0,y1,t);
  
  return (Nr==dNr) && (Nc==dNc) && GSTdoublesEqual(Xmin,x0) &&
         GSTdoublesEqual(Xmax,x1) && GSTdoublesEqual(Ymin,y0) &&
	 GSTdoublesEqual(Ymax,y1) && (grid_type==t);
}

float** GSTgriddedData::data() const
{
  return V;
}

float GSTgriddedData::qmin() const
{
  float zmin = GST_EMPTY_VALUE;
  int i,j;
  for (i=0;i<Nr;i++)
    for (j=0;j<Nc;j++)
      if ((V[i][j]<GST_EMPTY_LIMIT) && (V[i][j]<zmin)) 
        zmin=V[i][j];
  return zmin;
}

float GSTgriddedData::qmin(float x0, float x1, float y0, float y1) const
{
  float zmin = GST_EMPTY_VALUE;
  int i,j, i0, i1, j0, j1;
  
  if (grid_type==node)
  {
    i0=(int)floor((x0-Xmin)/dX)+1;
    i1=(int)floor((x1-Xmin)/dX);
    j0=(int)floor((y0-Ymin)/dY)+1;
    j1=(int)floor((y1-Ymin)/dY);
  }
  else
  {
    i0=(int)floor((x0-Xmin-0.5*dX)/dX)+1;
    i1=(int)floor((x1-Xmin-0.5*dX)/dX);
    j0=(int)floor((y0-Ymin-0.5*dY)/dY)+1;
    j1=(int)floor((y1-Ymin-0.5*dY)/dY);
  }
  if (i0<0) i0=0; if (i0>=Nr) return GST_EMPTY_VALUE;
  if (j0<0) j0=0; if (j0>=Nc) return GST_EMPTY_VALUE;
  if (i1<0) return GST_EMPTY_VALUE; if (i1>=Nr) i1=Nr-1;
  if (j1<0) return GST_EMPTY_VALUE; if (j1>=Nc) j1=Nc-1;
  
  for (i=i0;i<=i1;i++)
    for (j=j0;j<=j1;j++)
      if ((V[i][j]<GST_EMPTY_LIMIT) && (V[i][j]<zmin)) 
        zmin=V[i][j];
  return zmin;
}

float GSTgriddedData::qmax() const
{
  float zmax = -GST_EMPTY_VALUE;
  int i,j;
  for (i=0;i<Nr;i++)
    for (j=0;j<Nc;j++)
      if ((V[i][j]<GST_EMPTY_LIMIT) && (V[i][j]>zmax)) 
        zmax=V[i][j];
  return zmax;
}

float GSTgriddedData::qmax(float x0, float x1, float y0, float y1) const
{
  float zmax = -GST_EMPTY_VALUE;
  int i,j, i0, i1, j0, j1;
  
  if (grid_type==node)
  {
    i0=(int)floor((x0-Xmin)/dX)+1;
    i1=(int)floor((x1-Xmin)/dX);
    j0=(int)floor((y0-Ymin)/dY)+1;
    j1=(int)floor((y1-Ymin)/dY);
  }
  else
  {
    i0=(int)floor((x0-Xmin-0.5*dX)/dX)+1;
    i1=(int)floor((x1-Xmin-0.5*dX)/dX);
    j0=(int)floor((y0-Ymin-0.5*dY)/dY)+1;
    j1=(int)floor((y1-Ymin-0.5*dY)/dY);
  }
  if (i0<0) i0=0; if (i0>=Nr) return GST_EMPTY_VALUE;
  if (j0<0) j0=0; if (j0>=Nc) return GST_EMPTY_VALUE;
  if (i1<0) return GST_EMPTY_VALUE; if (i1>=Nr) i1=Nr-1;
  if (j1<0) return GST_EMPTY_VALUE; if (j1>=Nc) j1=Nc-1;
  
  for (i=i0;i<=i1;i++)
    for (j=j0;j<=j1;j++)
      if ((V[i][j]<GST_EMPTY_LIMIT) && (V[i][j]>zmax)) 
        zmax=V[i][j];
  return zmax;
}


float GSTgriddedData::qmidl(float x0, float x1, float y0, float y1) const
{
  float zmidl = 0.0;
  int count=0;
  int i,j, i0, i1, j0, j1;
  
  if (grid_type==node)
  {
    i0=(int)floor((x0-Xmin)/dX)+1;
    i1=(int)floor((x1-Xmin)/dX);
    j0=(int)floor((y0-Ymin)/dY)+1;
    j1=(int)floor((y1-Ymin)/dY);
  }
  else
  {
    i0=(int)floor((x0-Xmin-0.5*dX)/dX)+1;
    i1=(int)floor((x1-Xmin-0.5*dX)/dX);
    j0=(int)floor((y0-Ymin-0.5*dY)/dY)+1;
    j1=(int)floor((y1-Ymin-0.5*dY)/dY);
  }
  if (i0<0) i0=0; if (i0>=Nr) return GST_EMPTY_VALUE;
  if (j0<0) j0=0; if (j0>=Nc) return GST_EMPTY_VALUE;
  if (i1<0) return GST_EMPTY_VALUE; if (i1>=Nr) i1=Nr-1;
  if (j1<0) return GST_EMPTY_VALUE; if (j1>=Nc) j1=Nc-1;
  
  for (i=i0;i<=i1;i++)
    for (j=j0;j<=j1;j++)
      if (V[i][j]<GST_EMPTY_LIMIT )
      { 
        zmidl=zmidl+V[i][j];
	count=count+1;
      }	
  zmidl=zmidl/count;    
  return zmidl;
}


void GSTgriddedData::minimize()
{
  resize(1,1);
}

void GSTgriddedData::resize(int N_row, int N_col)
{
  if ((N_row<=0) || (N_col<=0)) GSTerror("Error resizing GSTgriddedData object");
  Nr=N_row;
  Nc=N_col;
  
  comment[0]='\0';
  
  switch (grid_type)
  {
    case (node): 
      if (Nr!=1)
        dX=(Xmax-Xmin)/(Nr-1);
      else
        dX=GST_EMPTY_VALUE;
      if (Nc!=1)
	dY=(Ymax-Ymin)/(Nc-1);
      else
        dY=GST_EMPTY_VALUE;
      break;
    case (area):
      dX=(Xmax-Xmin)/Nr;
      dY=(Ymax-Ymin)/Nc;
  }
  
  delete[] V[0]; delete[] V;
  
  V = new float*[Nr];
  if (V==NULL) GSTerror("Insufficient memory");
  V[0] = new float[Nr*Nc];
  if (V[0]==NULL) GSTerror("Insufficient memory");
  for (int i=1;i<Nr;i++) V[i]=V[i-1]+Nc;
  
}

void GSTgriddedData::resize
(int N_row, int N_col, double X_min, double X_max, double Y_min, double Y_max, GSTgrid_type type)
{
  if ((N_row<=0) || (N_col<=0)) GSTerror("Error resizing GSTgriddedData object");
  Nr=N_row;
  Nc=N_col;
  
  Xmin=X_min;
  Xmax=X_max;
  Ymin=Y_min;
  Ymax=Y_max;
  
  grid_type=type;
  
  comment[0]='\0';
  
  switch (type)
  {
    case (node): 
      if (Nr!=1)
        dX=(Xmax-Xmin)/(Nr-1);
      else
        dX=GST_EMPTY_VALUE;
      if (Nc!=1)
	dY=(Ymax-Ymin)/(Nc-1);
      else
        dY=GST_EMPTY_VALUE;
      break;
    case (area):
      dX=(Xmax-Xmin)/Nr;
      dY=(Ymax-Ymin)/Nc;
  }
  
  delete[] V[0]; delete[] V;
  
  V = new float*[Nr];
  if (V==NULL) GSTerror("Insufficient memory");
  V[0] = new float[Nr*Nc];
  if (V[0]==NULL) GSTerror("Insufficient memory");
  for (int i=1;i<Nr;i++) V[i]=V[i-1]+Nc;
}

void GSTgriddedData::resize(const GSTgriddedData &spec)
{
  spec.size(Nr,Nc);
  spec.getArea(Xmin,Xmax,Ymin,Ymax,grid_type);
  spec.getIntervals(dX,dY);
  spec.getComment(comment);
    
  delete[] V[0]; delete[] V;
  
  V = new float*[Nr];
  if (V==NULL) GSTerror("Insufficient memory");
  V[0] = new float[Nr*Nc];
  if (V[0]==NULL) GSTerror("Insufficient memory");
  for (int i=1;i<Nr;i++) V[i]=V[i-1]+Nc;  
}

GSTgriddedData &GSTgriddedData::operator=(const GSTgriddedData &rhs)
{
  if (this==&rhs) return *this;
  
  resize(rhs);
  float **pnt = rhs.data();
  for (int i=0;i<Nr;i++)
    for (int j=0;j<Nc;j++)
      V[i][j]=pnt[i][j];
      
  return *this;
}

GSTgriddedData operator+(const GSTgriddedData &op1, const GSTgriddedData &op2)
{
  int i, j, nr1, nc1, nr2, nc2;
  
  op1.size(nr1,nc1);
  op2.size(nr2,nc2);
  
  GSTgriddedData tmp(op1);
  
  if ((nr1!=nr2) || (nc1!=nc2))
    GSTerror("Attempt to use operator '+' with the gridded data operands of different size");
  
  for (i=0;i<nr1;i++)
    for (j=0;j<nc1;j++)
      tmp.V[i][j]=op1.V[i][j]+op2.V[i][j];
      
  return tmp;
}

GSTgriddedData operator+(float op1, const GSTgriddedData &op2)
{
  int i, j;
  GSTgriddedData tmp(op2);
    
  for (i=0;i<op2.Nr;i++)
    for (j=0;j<op2.Nc;j++)
      tmp.V[i][j]=op2.V[i][j]+op1;
      
  return tmp;
}

GSTgriddedData operator+(const GSTgriddedData &op1, float op2)
{
  int i, j;
  GSTgriddedData tmp(op1);
    
  for (i=0;i<op1.Nr;i++)
    for (j=0;j<op1.Nc;j++)
      tmp.V[i][j]=op1.V[i][j]+op2;
      
  return tmp;
}

GSTgriddedData operator-(const GSTgriddedData &op1, const GSTgriddedData &op2)
{
  int i, j, nr1, nc1, nr2, nc2;
  
  op1.size(nr1,nc1);
  op2.size(nr2,nc2);
  
  GSTgriddedData tmp(op1);
  
  if ((nr1!=nr2) || (nc1!=nc2))
    GSTerror("Attempt to use operator '-' with the gridded data operands of different size");
  
  for (i=0;i<nr1;i++)
    for (j=0;j<nc1;j++)
      tmp.V[i][j]=op1.V[i][j]-op2.V[i][j];
      
  return tmp;
}

GSTgriddedData operator-(float op1, const GSTgriddedData &op2)
{
  int i, j;
  GSTgriddedData tmp(op2);
    
  for (i=0;i<op2.Nr;i++)
    for (j=0;j<op2.Nc;j++)
      tmp.V[i][j]=op2.V[i][j]-op1;
      
  return tmp;
}

GSTgriddedData operator-(const GSTgriddedData &op1, float op2)
{
  int i, j;
  GSTgriddedData tmp(op1);
    
  for (i=0;i<op1.Nr;i++)
    for (j=0;j<op1.Nc;j++)
      tmp.V[i][j]=op1.V[i][j]-op2;
      
  return tmp;
}

GSTgriddedData operator*(const GSTgriddedData &op1, const GSTgriddedData &op2)
{
  int i, j, nr1, nc1, nr2, nc2;
  
  op1.size(nr1,nc1);
  op2.size(nr2,nc2);
  
  GSTgriddedData tmp(op1);
  
  if ((nr1!=nr2) || (nc1!=nc2))
    GSTerror("Attempt to use operator '*' with the gridded data operands of different size");
  
  for (i=0;i<nr1;i++)
    for (j=0;j<nc1;j++)
      tmp.V[i][j]=op1.V[i][j]*op2.V[i][j];
      
  return tmp;
}

GSTgriddedData operator*(float op1, const GSTgriddedData &op2)
{
  int i, j;
  GSTgriddedData tmp(op2);
    
  for (i=0;i<op2.Nr;i++)
    for (j=0;j<op2.Nc;j++)
      tmp.V[i][j]=op2.V[i][j]*op1;
      
  return tmp;
}

GSTgriddedData operator*(const GSTgriddedData &op1, float op2)
{
  int i, j;
  GSTgriddedData tmp(op1);
    
  for (i=0;i<op1.Nr;i++)
    for (j=0;j<op1.Nc;j++)
      tmp.V[i][j]=op1.V[i][j]*op2;
      
  return tmp;
}

GSTgriddedData operator/(const GSTgriddedData &op1, const GSTgriddedData &op2)
{
  int i, j, nr1, nc1, nr2, nc2;
  
  op1.size(nr1,nc1);
  op2.size(nr2,nc2);
  
  GSTgriddedData tmp(op1);
  
  if ((nr1!=nr2) || (nc1!=nc2))
    GSTerror("Attempt to use operator '/' with the gridded data operands of different size");
  
  for (i=0;i<nr1;i++)
    for (j=0;j<nc1;j++)
    if (GSTnonZero(op2.V[i][j]))
      tmp.V[i][j]=op1.V[i][j]/op2.V[i][j];
    else
    {
      tmp.V[i][j]=GST_EMPTY_VALUE;
      GSTwarning("Attempt to divide by zero");
    }
      
  return tmp;
}

GSTgriddedData operator/(float op1, const GSTgriddedData &op2)
{
  int i, j;
  GSTgriddedData tmp(op2);
    
  for (i=0;i<op2.Nr;i++)
    for (j=0;j<op2.Nc;j++)
    if (GSTnonZero(op2.V[i][j]))
      tmp.V[i][j]=op1/op2.V[i][j];
    else
    {
      tmp.V[i][j]=GST_EMPTY_VALUE;
      GSTwarning("Attempt to divide by zero");
    }
      
  return tmp;
}

GSTgriddedData operator/(const GSTgriddedData &op1, float op2)
{
  int i, j;
  GSTgriddedData tmp(op1);
    
  for (i=0;i<op1.Nr;i++)
    for (j=0;j<op1.Nc;j++)
    if (GSTnonZero(op2))
      tmp.V[i][j]=op1.V[i][j]/op2;
    else
    {
      tmp.V[i][j]=GST_EMPTY_VALUE;
      GSTwarning("Attempt to divide grid by zero constant");
    }
      
  return tmp;
}

/*------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------*/
/*-------------------------------Realisation of GSTscatterData memebers---------------------------------*/
/*------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------------------*/

GSTscatterData::GSTscatterData(const char *fields, float x_const, float y_const, float z_const)
{
  const_x_=x_const;
  const_y_=y_const;
  const_z_=z_const;
  Nf_=0;
  Np_=0;
  keys_.clear();
 _Br=1;
 _Bc=1;
 _sort=0;
 
 //_Num=new int*[1];
 //_Num[0]=new int[1];
  
  if (strchr(fields,'x')!=NULL) { 
    GSTdataVector *vec = new GSTdataVector;
    V['x']=vec;
    xvar_=true;
    Nf_++;
    keys_.push_back('x');
    
    _Nmin['x']=-1;
    _Nmax['x']=-1;
  }
  else xvar_=false;
  
  
  if (strchr(fields,'y')!=NULL) { 
    GSTdataVector *vec = new GSTdataVector;
    V['y']=vec;
    yvar_=true;
    Nf_++;
    keys_.push_back('y');
    _Nmin['y']=-1;
    _Nmax['y']=-1;
  }
  else yvar_=false;
  
  if (strchr(fields,'z')!=NULL) { 
    GSTdataVector *vec = new GSTdataVector;
    V['z']=vec;
    zvar_=true;
    Nf_++;
    keys_.push_back('z');
    _Nmin['z']=-1;
    _Nmax['z']=-1;
  }
  else zvar_=false;
  
  
  char ch;
  for (unsigned i=0;i<strlen(fields);i++) {
    ch=fields[i];
    if ((ch!='x') && (ch!='y') && (ch!='z'))
    {
      GSTdataVector *vec = new GSTdataVector;
      V[ch]=vec;
      _Nmin[ch]=-1;
      _Nmax[ch]=-1;
      
      Nf_++;
      keys_.push_back(ch);
    }
  }
}

GSTscatterData::~GSTscatterData()
{
 
  GSTscatterMap::iterator it;
  for (it=V.begin();it!=V.end();it++)
  {
    it->second->clear();
    delete it->second;
  }
    
  V.clear();
  Np_=Nf_=0;  
  
  if(_sort)
   {
    delete[] _Num[0];
    delete[] _Num;    
   }
}

float GSTscatterData::x(int i) const
{
  if ((i>=0) && (i<Np_))
    if (xvar_) return ((V.find('x'))->second)->at(i);
    else return const_x_;
  else return GST_EMPTY_VALUE;
}

float GSTscatterData::y(int i) const
{
  if ((i>=0) && (i<Np_))
    if (yvar_) return ((V.find('y'))->second)->at(i);
    else return const_y_;
  else return GST_EMPTY_VALUE;
}

float GSTscatterData::z(int i) const
{
  if ((i>=0) && (i<Np_))
    if (zvar_) return ((V.find('z'))->second)->at(i);
    else return const_z_;
  else return GST_EMPTY_VALUE;
}

void GSTscatterData::xyz(int i, float &xv, float &yv, float &zv) const
{
  xv = x(i); yv = y(i); zv = z(i);
}

void GSTscatterData::xyz(int i, GSTxyzType &p) const
{
  p.x = x(i); p.y = y(i); p.z = z(i);
}

int GSTscatterData::setX(int i, float xv)
{
  if ((i>=0) && (i<Np_))
    if (xvar_)
     {
      if(_Nmin['x']!=-1)
       {
        if(xv <= ((*V['x'])[_Nmin['x']])) _Nmin['x']=i;
	else if(_Nmin['x']==i) _Nmin['x']=-1;
       }
      if(_Nmax['x']!=-1 )
       {
        if((xv >= ((*V['x'])[_Nmax['x']])) && (xv < GST_EMPTY_LIMIT)) _Nmax['x']=i;
	else if(_Nmax['x']==i) _Nmax['x']=-1;
       } 
 
     
     (*V['x'])[i]=xv;
      return 1;     
     }
  return 0;
}

int GSTscatterData::setY(int i, float yv)
{
  if ((i>=0) && (i<Np_))
    if (yvar_) 
     {
      if(_Nmin['y']!=-1 )
       {
        if(yv<=((*V['y'])[_Nmin['y']])) _Nmin['y']=i;
	else  if(_Nmin['y']==i) _Nmin['y']=-1;
       }
      if(_Nmax['y']!=-1 )
       {
        if((yv>=((*V['y'])[_Nmax['y']])) && (yv<GST_EMPTY_LIMIT)) _Nmax['y']=i;
	else  if(_Nmax['y']==i) _Nmax['y']=-1;
       }  
     
     
     
     (*V['y'])[i]=yv; return 1;}
  return 0;
}

int GSTscatterData::setZ(int i, float zv)
{
  if ((i>=0) && (i<Np_))
    if (zvar_)
     {
      if(_Nmin['z']!=-1 )
       {
        if(zv<=((*V['z'])[_Nmin['z']])) _Nmin['z']=i;
	else  if(_Nmin['z']==i) _Nmin['z']=-1;
       }
      if(_Nmax['z']!=-1 )
       {
        if((zv>=((*V['z'])[_Nmax['z']])) && (zv<GST_EMPTY_LIMIT )) _Nmax['z']=i;
	else  if(_Nmax['z']==i) _Nmax['z']=-1;
       } 
      
      (*V['z'])[i]=zv; return 1;}
  return 0;
}

int GSTscatterData::setF(int i, char key, float v)
{
  if ((i>=0) && (i<Np_))
    if (V.find(key)!=V.end())
     {
      if(_Nmin[key]!=-1 )
       {
        if(v<=((*V[key])[_Nmin[key]])) _Nmin[key]=i;
	else  if(_Nmin[key]==i) _Nmin[key]=-1;
       }
      if(_Nmax[key]!=-1 )
       {
        if((v>=((*V[key])[_Nmax[key]])) && (v<GST_EMPTY_LIMIT )) _Nmax[key]=i;
	else  if(_Nmax[key]==i) _Nmax[key]=-1;
       } 
      
      (*V[key])[i]=v; return 1;}
  return 0;
}

int GSTscatterData::setF( char key, float v)
{
  int i;
    if (V.find(key)!=V.end())
     {
      for(i=0;i<Np_;i++)      
       (*V[key])[i]=v;
       
       _Nmax[key]=-1;
       _Nmin[key]=-1;
       return 1;
      }
  return 0;
}


int GSTscatterData::setXYZ(int i, float xv, float yv, float zv)
{
  int n=0;
  if ((i>=0) && (i<Np_))
  {
    if (xvar_)     
     {
      if(_Nmin['x']!=-1 )
       {
        if(xv<=((*V['x'])[_Nmin['x']])) _Nmin['x']=i;
	else  if(_Nmin['x']==i) _Nmin['x']=-1;
       }
      if(_Nmax['x']!=-1 )
       {
        if((xv>=((*V['x'])[_Nmax['x']])) && (xv<GST_EMPTY_LIMIT)) _Nmax['x']=i;
	else  if(_Nmax['x']==i) _Nmax['x']=-1;
       }       
            
      (*V['x'])[i]=xv; n++;
     }
    if (yvar_)
     {
      if(_Nmin['y']!=-1 )
       {
        if(yv<=((*V['y'])[_Nmin['y']])) _Nmin['y']=i;
	else  if(_Nmin['y']==i) _Nmin['y']=-1;
       }
      if(_Nmax['y']!=-1 )
       {
        if((yv>=((*V['y'])[_Nmax['y']])) && (yv<GST_EMPTY_LIMIT)) _Nmax['y']=i;
	else  if(_Nmax['y']==i) _Nmax['y']=-1;
       }     
      
      (*V['y'])[i]=yv; n++;
      }
      
    if (zvar_) 
    { 
     if(_Nmin['z']!=-1 )
       {
        if(zv<=((*V['z'])[_Nmin['z']])) _Nmin['z']=i;
	else  if(_Nmin['z']==i) _Nmin['z']=-1;
       }
      if(_Nmax['z']!=-1 )
       {
        if((zv>=((*V['z'])[_Nmax['z']])) && (zv<GST_EMPTY_LIMIT )) _Nmax['z']=i;
	else  if(_Nmax['z']==i) _Nmax['z']=-1;
       }  
    
    (*V['z'])[i]=zv; n++;
    }
  }
  return n;
}

int GSTscatterData::setXYZ(int i, GSTxyzType &p)
{
  int n=0;
  if ((i>=0) && (i<Np_))
  {
    if (xvar_)
     { 
       if(_Nmin['x']!=-1)
        {
         if(p.x <= ((*V['x'])[_Nmin['x']])) _Nmin['x']=i;
	 else if(_Nmin['x']==i) _Nmin['x']=-1;
        }
       if(_Nmax['x']!=-1 )
        {
         if((p.x >= ((*V['x'])[_Nmax['x']])) && (p.x < GST_EMPTY_LIMIT)) _Nmax['x']=i;
	 else if(_Nmax['x']==i) _Nmax['x']=-1;
       } 
     (*V['x'])[i]=p.x; n++;
     }
    if (yvar_)
     {
      if(_Nmin['y']!=-1 )
       {
        if(p.y<=((*V['y'])[_Nmin['y']])) _Nmin['y']=i;
	else  if(_Nmin['y']==i) _Nmin['y']=-1;
       }
      if(_Nmax['y']!=-1 )
       {
        if((p.y>=((*V['y'])[_Nmax['y']])) && (p.y<GST_EMPTY_LIMIT)) _Nmax['y']=i;
	else  if(_Nmax['y']==i) _Nmax['y']=-1;
       }     
      
      (*V['y'])[i]=p.y; n++;
     }
    if (zvar_)
     {
      if(_Nmin['z']!=-1 )
       {
        if(p.z<=((*V['z'])[_Nmin['z']])) _Nmin['z']=i;
	else  if(_Nmin['z']==i) _Nmin['z']=-1;
       }
      if(_Nmax['z']!=-1 )
       {
        if((p.z>=((*V['z'])[_Nmax['z']])) && (p.z<GST_EMPTY_LIMIT )) _Nmax['z']=i;
	else  if(_Nmax['z']==i) _Nmax['z']=-1;
       }  
     
     (*V['z'])[i]=p.z; n++;
     }
  }
  return n;
}

void GSTscatterData::allF(int i, GSTpointMap &pm) const
{
  GSTscatterMap::const_iterator it;
  pm.clear();
  
  if ((i>=0) && (i<Np_))
  {
    for (it = V.begin();it!=V.end();it++)
      pm[it->first]=(*(it->second))[i];
  }
}

float GSTscatterData::f(int i, char key) const
{
  GSTscatterMap::const_iterator it;
  if ((i>=0) && (i<Np_))
  {
    it = V.find(key);
    if (it!=V.end())
      return (*(it->second))[i];
  }
  return GST_EMPTY_VALUE;
}

int GSTscatterData::setAllF(int i, const GSTpointMap &pm)
{
  GSTpointMap::const_iterator it;
  int n=0;
    
  if ((i>=0) && (i<Np_))
  {
    for (it = pm.begin();it!=pm.end();it++)
      if (V.find(it->first)!=V.end())
      {
      if(_Nmin[it->first]!=-1 )
       {
        if(it->second <=((*V[it->first])[_Nmin[it->first]])) _Nmin[it->first]=i;
	else  if(_Nmin[it->first]==i) _Nmin[it->first]=-1;
       }
      if(_Nmax[it->first]!=-1 )
       {
        if((it->second >=((*V[it->first])[_Nmax[it->first]])) && (it->second < GST_EMPTY_LIMIT ))
	  _Nmax[it->first]=i;
	else  if(_Nmax[it->first]==i) _Nmax[it->first]=-1;
       }       
      
        (*V[it->first])[i]=it->second;
        n++;
      }
  }
  return n;
}

int GSTscatterData::pntDel(int i) 
{
  GSTscatterMap::const_iterator it;
  if ((i>=0) && (i<Np_))
  {
    for (it=V.begin();it!=V.end();it++)
     { 
      if(_Nmax[it->first]==i) _Nmax[it->first]=-1;
      else if(_Nmax[it->first]>i) _Nmax[it->first]=_Nmax[it->first]-1;
      if(_Nmin[it->first]==i) _Nmin[it->first]=-1;
      else if(_Nmin[it->first]>i) _Nmin[it->first]=_Nmin[it->first]-1;
      
      V[it->first]->erase(V[it->first]->begin()+i);
      
     }
    Np_--; 
    return 1;
  }
  return 0;
}

int GSTscatterData::pntAdd(GSTpointMap &pm)
{
  int n=0;
  GSTscatterMap::iterator it;
  for (it=V.begin();it!=V.end();it++)
    if (pm.find(it->first)!=pm.end())
    {
    
      V[it->first]->push_back(pm[it->first]);
      n++;
      //min-max
      if((Np_==0) || (_Nmin[it->first]==-2))
       {
        if(pm[it->first]< GST_EMPTY_LIMIT) {_Nmin[it->first]=Np_; _Nmax[it->first]=Np_;}
	else{_Nmin[it->first]=-2; _Nmax[it->first]=-2;}       
       }
      else
       {
        if((pm[it->first]< GST_EMPTY_LIMIT) && (_Nmin[it->first]>=0) && 
	   ( pm[it->first]<(*V[it->first])[_Nmin[it->first]])) _Nmin[it->first]=Np_;
	   
	if((pm[it->first]< GST_EMPTY_LIMIT) && (_Nmax[it->first]>=0) && 
	   ( pm[it->first]>(*V[it->first])[_Nmax[it->first]])) _Nmax[it->first]=Np_;        
	}
	      
     }
     else
     {
       V[it->first]->push_back(GST_EMPTY_VALUE);
     }
  Np_++;
  return n;
}

int GSTscatterData::numPnt() const
{
  return Np_;
}

int GSTscatterData::numFld() const
{
  return Nf_;
}

void GSTscatterData::resize(int N, float val) 
{
  GSTscatterMap::iterator it;
  for (it=V.begin();it!=V.end();it++)
  {
    if(N>Np_ && val>(*V[it->first])[_Nmax[it->first]])  _Nmax[it->first]=Np_;
    if(N>Np_ && val<(*V[it->first])[_Nmin[it->first]])  _Nmin[it->first]=Np_;
    if(N<Np_ && (_Nmax[it->first]+1)>N) _Nmax[it->first]=-1;
    if(N<Np_ && (_Nmin[it->first]+1)>N) _Nmin[it->first]=-1;
    it->second->resize(N,val);
  }  
  Np_=N;
}

void GSTscatterData::reserve(int N)
{
  GSTscatterMap::iterator it;
  for (it=V.begin();it!=V.end();it++)
    it->second->reserve(N);
}

GSTdataVector* GSTscatterData::field(char key)
{
  if (V.find(key)!=V.end())
    return V[key];
  return NULL;
}

float* GSTscatterData::field_plain(char key)
{
  if (V.find(key)!=V.end())
    return &((*(V[key]))[0]);
  return NULL;
}

int GSTscatterData::fieldDel(char key)
{
  if (V.find(key)!=V.end())
  {
    switch (key)
    {
      case 'x': xvar_=false; const_x_=GST_EMPTY_VALUE; break;
      case 'y': yvar_=false; const_y_=GST_EMPTY_VALUE; break;
      case 'z': zvar_=false; const_z_=GST_EMPTY_VALUE;
    }
    V[key]->clear();
    delete V[key];
    V.erase(key);
    _Nmin.erase(key);
    _Nmax.erase(key);
    return 1;
  }
  else
  {
    switch (key)
    {
      case 'x': xvar_=false; const_x_=GST_EMPTY_VALUE; return 1;
      case 'y': yvar_=false; const_y_=GST_EMPTY_VALUE; return 1;
      case 'z': zvar_=false; const_z_=GST_EMPTY_VALUE; return 1;
      default: return 0;
    }
  }
}

int GSTscatterData::fieldIns(char key, float val)
{
  if (V.find(key)!=V.end()) return 0;
  GSTdataVector *ptr = new GSTdataVector;
  V[key]=ptr;
  V[key]->resize(Np_,val);
  _Nmin[key]=-1;
  _Nmax[key]=-1;
  Nf_++;
  return 1;
}

int GSTscatterData:: namesFLD(std::string &name) const
{
  name.erase();
  GSTscatterMap::const_iterator it;
  for (it=V.begin();it!=V.end();it++)
    name+=it->first;
  return Nf_;  
}

void GSTscatterData::chknp()
{
  GSTscatterMap::iterator it;
  for (it=V.begin();it!=V.end();it++)
    if (Np_!=(int)it->second->size())
      GSTerror("Size problem...\n"); 
}

float GSTscatterData::maxX(int *n)
{
 if (xvar_)
  {
   if(_Nmax['x']>=0) { if(n!=NULL) *n=_Nmax['x'];  return (*V['x'])[_Nmax['x']];}
   else
   { 
   
    float max=(*V['x'])[0];
    _Nmax['x']=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V['x'])[i]> max){ max=(*V['x'])[i]; _Nmax['x']=i; }
    }
    if(n!=NULL) *n=_Nmax['x']; 
    return max;
   }
  } 
 else {*n=0; return const_x_;} 
}

float GSTscatterData::minX(int *n)
{
 if (xvar_)
  {
   if(_Nmin['x']>=0) { if(n!=NULL) *n=_Nmin['x'];  return (*V['x'])[_Nmin['x']];}
   else
   { 
   
    float min=(*V['x'])[0];
    _Nmin['x']=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V['x'])[i]< min){ min=(*V['x'])[i]; _Nmin['x']=i; }
    }
    if(n!=NULL) *n=_Nmin['x']; 
    return min;
   }
  } 
 else{ *n=0; return const_x_;} 
}

float GSTscatterData::maxY(int *n)
{
 if (yvar_)
  {
   if(_Nmax['y']>=0) { if(n!=NULL) *n=_Nmax['y'];  return (*V['y'])[_Nmax['y']];}
   else
   { 
   
    float max=(*V['y'])[0];
    _Nmax['y']=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V['y'])[i]> max){ max=(*V['y'])[i]; _Nmax['y']=i; }
    }
    if(n!=NULL) *n=_Nmax['y']; 
    return max;
   }
  } 
 else {*n=0; return const_y_;} 
}


float GSTscatterData::minY(int *n)
{
 if (yvar_)
  {
   if(_Nmin['y']>=0) { if(n!=NULL) *n=_Nmin['y'];  return (*V['y'])[_Nmin['y']];}
   else
   { 
   
    float min=(*V['y'])[0];
    _Nmin['y']=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V['y'])[i]< min){ min=(*V['y'])[i]; _Nmin['y']=i; }
    }
    if(n!=NULL) *n=_Nmin['y']; 
    return min;
   }
  } 
 else  {*n=0; return const_y_;} 
}

float GSTscatterData::maxZ(int *n)
{
 if (zvar_)
  {
   if(_Nmax['z']>=0) { if(n!=NULL) *n=_Nmax['z'];  return (*V['z'])[_Nmax['z']];}
   else
   { 
   
    float max=(*V['z'])[0];
    _Nmax['z']=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V['z'])[i]> max){ max=(*V['z'])[i]; _Nmax['z']=i; }
    }
    if(n!=NULL) *n=_Nmax['z']; 
    return max;
   }
  } 
 else {*n=0; return const_z_;}
}


float GSTscatterData::minZ(int *n)
{
 if (zvar_)
  {
   if(_Nmin['z']>=0) { if(n!=NULL) *n=_Nmin['z'];  return (*V['z'])[_Nmin['z']];}
   else
   { 
   
    float min=(*V['z'])[0];
    _Nmin['z']=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V['z'])[i]< min){ min=(*V['z'])[i]; _Nmin['z']=i; }
    }
    if(n!=NULL) *n=_Nmin['z']; 
    return min;
   }
  } 
 else {*n=0; return const_z_;}
 }

float  GSTscatterData::minF(char key, int *n)
{
 if(V.find(key)==V.end()) {*n=0; return GST_EMPTY_VALUE ; 
   GSTwarning("Attempt to find minimum for non-existing data field");}
 if( (key=='x') && (!xvar_)) {*n=0; return const_x_;}
 if( (key=='y') && (!yvar_)) {*n=0; return const_y_;}
 if( (key=='z') && (!zvar_)) {*n=0; return const_z_;}
 
 if(_Nmin[key]>=0) { if(n!=NULL) *n=_Nmin[key];  return (*V[key])[_Nmin[key]];}
   else
   { 
   
    float min=(*V[key])[0];
    _Nmin[key]=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V[key])[i]< min){ min=(*V[key])[i]; _Nmin[key]=i; }
    }
    if(n!=NULL) *n=_Nmin[key]; 
    return min;
   }
}


float GSTscatterData::maxF(char key,int *n)
{
 if(V.find(key)==V.end()) {*n=0; return GST_EMPTY_VALUE ;
   GSTwarning("Attempt to find maximum for non-existing data field");}
 if( (key=='x') && (!xvar_)) {*n=0; return const_x_;}
 if( (key=='y') && (!yvar_)) {*n=0; return const_y_;}
 if( (key=='z') && (!zvar_)) {*n=0; return const_z_;}
 if(_Nmax[key]>=0) { if(n!=NULL) *n=_Nmax[key];  return (*V[key])[_Nmax[key]];}
   else
   { 
   
    float max=(*V[key])[0];
    _Nmax[key]=0;
    int i;
    for(i=1;i<Np_;i++)
    {
     if((*V[key])[i]> max){ max=(*V[key])[i]; _Nmax[key]=i; }
    }
    if(n!=NULL) *n=_Nmax[key]; 
    return max;
   }

}

void GSTscatterData::max(GSTpointMap &p, GSTpointNum *n)
{
  int f;
  GSTscatterMap::iterator it;
  for (it=V.begin();it!=V.end();it++)
   {
     f=0;
     if(V.find(it->first)==V.end()) f=-1;
     if(it->first=='x' && (!xvar_)) f=1;
     if(it->first=='y' && (!yvar_)) f=2;
     if(it->first=='z' && (!zvar_)) f=3;        
     switch(f)
       {       
       case 0:            
        if(_Nmax[it->first]>=0) { if(n!=NULL) (*n)[it->first]=_Nmax[it->first]; 
                               p[it->first]=(*V[it->first])[_Nmax[it->first]];}
         else
          {   
           float max=(*V[it->first])[0];
	   _Nmax[it->first]=0;
           int i;
           for(i=1;i<Np_;i++)
            {
             if((*V[it->first])[i] > max){ max=(*V[it->first])[i]; _Nmax[it->first]=i; }
            }
           if(n!=NULL) (*n)[it->first]=_Nmax[it->first];  p[it->first]=max;
          }
	 break;
       case -1: (*n)[it->first]=0; p[it->first]=GST_EMPTY_VALUE; 
                GSTwarning("Attempt to find maximum for non-existing data field");
		break; 
       case 1:(*n)['x']=0; p['x']=const_x_; break;
       case 2:(*n)['y']=0; p['y']=const_y_; break;
       case 3:(*n)['z']=0; p['z']=const_z_; break;          
      }  
     }   
 }


void GSTscatterData::min(GSTpointMap &p, GSTpointNum *n)
{ 
  int f;
  GSTscatterMap::iterator it;
  for (it=V.begin();it!=V.end();it++)
     {
     f=0;
     if(V.find(it->first)==V.end()) f=-1;
     if(it->first=='x' && (!xvar_)) f=1;
     if(it->first=='y' && (!yvar_)) f=2;
     if(it->first=='z' && (!zvar_)) f=3;        
     switch(f)
       {       
       case 0:        
        if(_Nmin[it->first]>=0) { if(n!=NULL) (*n)[it->first]=_Nmin[it->first]; 
                                p[it->first]=(*V[it->first])[_Nmin[it->first]];}
        else
         {    
          float min=(*V[it->first])[0];
	  _Nmin[it->first]=0;
          int i;
          for(i=1;i<Np_;i++)
           {
            if((*V[it->first])[i]< min){ min=(*V[it->first])[i]; _Nmin[it->first]=i; }
            }
         if(n!=NULL)(*n)[it->first]=_Nmin[it->first];
         p[it->first]=min; 
	} 
	 break;
       case -1: (*n)[it->first]=0; p[it->first]=GST_EMPTY_VALUE; 
                GSTwarning("Attempt to find minimum for non-existing data field"); 
		break;
       case 1:(*n)['x']=0; p['x']=const_x_; break;
       case 2:(*n)['y']=0; p['y']=const_y_; break;
       case 3:(*n)['z']=0; p['z']=const_z_; break;       
      }	 	 
   }    
}

void GSTscatterData:: sort(int Nx, int Ny)

{
  if(!xvar_)GSTerror("There is no field x \n");
  if(!yvar_)GSTerror("There is no field y \n");

  int i,j,k;
  float X0,Y0;
  GSTpointMap A1,A2;

  X0=minX()-GST_FLOAT_ACCURACY;  
  Y0=minY()-GST_FLOAT_ACCURACY; 

  _Br=Nx;
  _Bc=Ny;

  _dx=(maxX()-X0)/Nx;
  _dy=(maxY()-Y0)/Ny;

 
 if(!_sort) if(!fieldIns('_')) GSTerror( " name '_'  \n");
 
 if (_sort)
  {
    delete[] _Num[0];
    delete[] _Num;
  }
  
 _Num=new int*[Nx];
 _Num[0]=new int[Nx*Ny];
 for(i=1; i<Nx; i++)
  {
   _Num[i]=_Num[0]+i*Ny;
  }
  
  for(i=0;i<Nx ;i++)
   {
    for(j=0; j<Ny; j++)
    {
     _Num[i][j]=0;
    }
   }
  
 
 for(k=0;k<Np_;k++)
  {
    i=(int)floor(((*V['x'])[k]-X0)/_dx);
    if(i==Nx) i=Nx-1;
    j=(int)floor(((*V['y'])[k]-Y0)/_dy);
    if(j==Ny) j=Ny-1;
     (*V['_'])[k]=i*Ny+j;
     _Num[i][j]=_Num[i][j]+1; // сейчас тут- чило точек принадлежащих данному блоку      
   }
 
 

int fl;
int Npoint=0; //номер точки в который записывается текущее значение- выше отсортированный массив
              // ниже  не отсортированный.   
	         
 
  for(i=0; i<Nx; i++) // по номерам блоков 
   {
    for(j=0; j<Ny ;j++)
     {
      fl=0;
      for(k=Npoint; k<Np_ && fl<=_Num[i][j] ; k++) //по точкам
       { 	
	 if( (*V['_'])[k]==i*Ny+j)
	  {	  
	   if(k!=Npoint)
	    {
	      allF(Npoint,A1);
              allF(k,A2);
              setAllF(Npoint,A2);
              setAllF(k,A1);    
	    }
	   Npoint++;
	   fl++;
	  }
	 
	
       
       } 
     } // по точкам
   }
 
  
 _Num[0][0]=_Num[0][0]-1; 
    
 
  for(i=0; i< Nx; i++)
   {
   for(j=0; j<Ny; j++)
    { 
     if(i+j>0)
     {      
     _Num[i][j]=_Num[i][j-1]+_Num[i][j];  //номер последнего эллемента принадлежащего блоку
     
      }
    }     
   }
_sort=1;   
}





void GSTscatterData::sort(float Lx, float Ly, int* nx, int *ny)
{


if(!xvar_)GSTerror("The is no field x \n");
if(!yvar_)GSTerror("The is no field y \n");

_dx=Lx;
_dy=Ly;

int i,j,k;
GSTpointMap A1,A2;
float X1,X2;
float Y1,Y2;

X1=minX()-GST_FLOAT_ACCURACY;
X2=maxX()+GST_FLOAT_ACCURACY;
Y1=minY()-GST_FLOAT_ACCURACY;
Y2=maxY()+GST_FLOAT_ACCURACY;

if(((X2-X1)/Lx)- floor((X2-X1)/Lx)>0)
_Br=(int)floor((X2-X1)/Lx)+1;
else _Br=(int)floor((X2-X1)/Lx);

if(((Y2-Y1)/Ly)- floor((Y2-Y1)/Ly)>0)
_Bc= (int)floor((Y2-Y1)/Ly)+1;
else _Bc= (int)floor((Y2-Y1)/Ly);


*nx=_Br;
*ny=_Bc;

if(_Br==1 && _Bc==1)GSTwarning(" One block only in Array\n"); 
 if(!_sort) if(!fieldIns('_')) GSTerror( " name '_'  \n");
 
 if (_sort)
  {
    delete[] _Num[0];
    delete[] _Num;
  }
  
 _Num=new int*[_Br];
 _Num[0]=new int[_Br*_Bc];
 for(i=1; i<_Br; i++)
  {
   _Num[i]=_Num[0]+i*_Bc;
  }
  
  for(i=0;i<_Br ;i++)
   {
    for(j=0; j<_Bc; j++)
    {
     _Num[i][j]=0;
    }
   }
  
 
 for(k=0;k<Np_;k++)
  {
    i=(int)floor(((*V['x'])[k]-X1)/Lx);
    if(i==_Br) i=_Br-1;
    j=(int)floor(((*V['y'])[k]-Y1)/Ly);
    if(j==_Bc) j=_Bc-1;
     (*V['_'])[k]=i*_Bc+j;
     _Num[i][j]=_Num[i][j]+1; // сейчас тут- чило точек принадлежащих данному блоку
      
   } 

int fl;
int Npoint=0; //номер точки в который записывается текущее значение- выше отсортированный массив
              // ниже  не отсортированный.    

        	      
  for(i=0; i<_Br; i++) // по номерам блоков 
   {
    for(j=0; j<_Bc ;j++)
     {
      
       fl=0;
      for(k=Npoint; k<Np_ && fl<=_Num[i][j]; k++) //по точкам
       {	
	
	 if( (*V['_'])[k]==i*_Bc+j)
	  {	  
	   if(k>Npoint)
	    {	      
	      allF(Npoint,A1);
              allF(k,A2);
              setAllF(Npoint,A2);
              setAllF(k,A1);    
	    }	   
	   Npoint++;
	   fl=fl+1;
	   
	  }
	 
	
       
       } 
     } 
   } 
  
  
 _Num[0][0]=_Num[0][0]-1;    
 
  for(i=0; i< _Br; i++)
   {
   for(j=0; j<_Bc; j++)
    {   
     if(i+j>0)    
     _Num[i][j]=_Num[i][j-1]+_Num[i][j];  //номер последнего эллемента принадлежащего блоку
  }
}

_sort=1;
}

 
int GSTscatterData::np_block(int bi, int bj)
{ 
 //if(!_sort) GSTerror("Array is not sorted \n"); 
 //if(_Br==1 && _Bc==1 )   GSTwarning(" One blok ounly\n");
  
 //if( 0<= bi && bi<_Br && 0<=bj && bj< _Bc )
  //{
   if (bi==0 && bj==0) return (_Num[bi][bj]+1);
   return (_Num[bi][bj]-_Num[bi][bj-1]);
  // }  
  //else GSTerror("error number of block \n");
  //return -1; 
}

 void GSTscatterData:: N_block( float x, float y, int &bi, int &bj)
 {
  bi=(int)floor((x- minX()-GST_FLOAT_ACCURACY)/_dx);
  bj=(int)floor((y-minY()-GST_FLOAT_ACCURACY)/_dy); 
 }
  
bool  GSTscatterData::BlockExist(int bi, int bj)
{
 if( 0<= bi && bi<_Br && 0<=bj && bj< _Bc) return 1;
 else return 0;
}

bool GSTscatterData::BlockSort()
{
 return _sort;
}

void GSTscatterData::BlockNumber(int &nBx, int &nBy)
{
 nBx=_Br;
 nBy=_Bc;
}    

void GSTscatterData::Block_step(float &dx, float &dy)
{ 
 if (_dx==0 && _dy==0 )GSTwarning("Array is not sorted\n");
 dx=_dx;
 dy=_dy; 
}      

int GSTscatterData::bxyz(int bi, int bj, int i, float &xv, float &yv, float &zv)
{
 if(bi==0 && bj==0)
 {
  xv=(*V['x'])[i];
  yv=(*V['y'])[i];
  zv=(*V['z'])[i];
  return i;
 }
 else 
 {
  xv=(*V['x'])[_Num[bi][bj-1]+1+i];
  yv=(*V['y'])[_Num[bi][bj-1]+1+i];
  zv=(*V['z'])[_Num[bi][bj-1]+1+i];
  return _Num[bi][bj-1]+1+i;
 }
}


float GSTscatterData:: lbf(int bi, int bj, int i, char key)
{
  if(bi==0 && bj==0)      
     return f((i),key);     
   else      
     return f((_Num[bi][bj-1]+1+i),key);  
}


float GSTscatterData:: bf(int bi, int bj, int i, char key)
{
  if(!_sort) GSTerror("Array is not sorted \n"); 
  if(_Br==1 && _Bc==1) GSTwarning(" One blok ounlyin Array\n");
  if( 0<= bi && bi<_Br && 0<=bj && bj< _Bc )
  {
    if(bi==0 && bj==0)
     { 
      if(i<=_Num[0][0])
      return f((i),key);//???
      else GSTerror("number of point in block [0][0] is too big\n");
      }
     
    else
     { 
      if((_Num[bi][bj-1]+1+i)<=_Num[bi][bj])
       return f((_Num[bi][bj-1]+1+i),key);
       else GSTerror("number of point in block is too big\n");        
     }
   }  
  else GSTerror("error number of block\n ");
  return GST_EMPTY_VALUE;
}



void GSTscatterData::allBf(int bi, int bj, int i, GSTpointMap &pm)
{
 if(!_sort) GSTerror("Array is not sorted \n"); 
 if(_Br==1 && _Bc==1) GSTwarning(" One blok ounly in Array\n");
 if( 0<= bi && bi<_Br && 0<=bj && bj< _Bc )
  {  
   if(bi==0 && bj==0)
     { 
      if(i<=_Num[0][0])
      return allF((i),pm); //?
      else GSTerror("number of point in block [0][0] is too big\n");
      }  
   else
   {  
     if((_Num[bi][bj-1]+1+i)<=_Num[bi][bj])   
     allF( (_Num[bi][bj-1]+1+i), pm);
     else GSTerror("number of point in block is too big\n");
   }
  }
  else GSTerror("error number of block\n ");
}





   
