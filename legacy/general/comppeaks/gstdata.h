/*������ �� 29.01.06, �� ��������� �������������*/
#ifndef GSTGRIDDEDDATA_H
#define GSTGRIDDEDDATA_H
#include <string>
#include <vector>
#include <map>

#include "gstdefs.h"

enum GSTgrid_type {node, area};

typedef std::vector<float> GSTdataVector;
typedef std::map<char,float> GSTpointMap;
typedef std::map<char,GSTdataVector*> GSTscatterMap;
typedef std::map<char,int> GSTpointNum;

typedef struct
	{
	  float x,y,z;
	}
	GSTxyzType;

class GSTgriddedData
{
  public:
   
    GSTgriddedData(int N_row=1, int N_col=1, double X_min=0, double X_max=1, 
                   double Y_min=0, double Y_max=1, GSTgrid_type type=node);
                   
    /* ����������� ������.
       ������� ������ ������ GSTgriddedData. �.�. ���������� ������ ��� ������� �������� 
       N_row X N_col.
       �������� ��������� ���������� ��������� � ����������� ���� ������ � ������ � 
       ��� ��������������  � ������� ����-�������, ��������� ����.
              
       ��� ������� �� ��������� GSTgriddedData() ����� ��������  ��������� �������� ����������.
       ����� �������� �������� ������ ���������� ������ ����������, � �������� �������� ����������� 
      (�� ���������)���������� ����� ������������ �������� �� ���������.
       GSTgriddedData(5,5).
       
       �������� GSTgrid_type type ����� ��������� �������� node(��������, ������� ����� 
       ����� �������� � ���  ������� ����� ��������� � ����� �����) � area (�������� 
       ����� ��������� � ������� ��������)*/
    
    GSTgriddedData(const GSTgriddedData &rhs);
    /* �����������-����������, �������������� �������� �����������. ����� ������ ��������.
    */             
                   
    ~GSTgriddedData();   /*���������� ������*/ 
       
    /*�������� ����-������� ������ */
    
     void size(int &N_row, int &N_col) const;
    /* �� ���������� ���� GSTgriddedData ������ ���� ���� N_row � ����� �������� N_col.
     ���������� ������� �������� ����������, � ������� ������������ ��������� */
     
    void setArea(double X_min, double X_max, double Y_min, double Y_max, GSTgrid_type type);
    /*� GSTgriddedData ���������� �������� ���������� �� ���������� �����*/
    
    int getArea(double &X_min, double &X_max, double &Y_min, double &Y_max, GSTgrid_type &type) const;
    /* �� GSTgriddedData ������ �������� ���������� ����� � ���������� �� � ���������� ��
     ����������*/
     
    void setType(GSTgrid_type type);
    /*������������� ��� ����� � type � ������������� �������� ����� dX � dY*/
    
    int getIntervals(double &d_X, double &d_Y) const;
    /* �� GSTgriddedData ������ �������� � ���������� �� � ���������� �� ����������*/
    
    float value(int r, int c) const;
    /* ������ �������� ������� � r-��� ������ c -� ������� */
    
    float value(double x, double y) const;
    /* ������ �������� ������� � ����� ��������� � x,y ��� �,y ������ � ������������� �������
     ��������� */
    
    void setValue(float v);
    /* ����� ��������� v �� ��� ������ ������� */
            
    void setValue(int r, int c, float v);
    /*����� �������� v � ������ [r][c]*/
    
    void setValue(double x, double y, float v);
    /* ����� �������� v � ������ ,��������� � x,y ��� �,y ������ � ������������� ������� ���������
     */
    
    char *getComment(char *comm) const;
    /* ������ ����������� �� GSTgriddedData.*/
    
    void setComment(const char *comm);
    /* ����� ����������� � GSTgriddedData.*/
    
    float **data() const;
    /* ���������� ��������� (�����) ������ ������� �������, �� ��������� �������� �����
     ��������������� 
    ���������� � �������. ������ �������������: (*D).data()[i][j] */
    
    bool haveEqualArea(GSTgriddedData &d) const; 
    /*���������� X_min, X_max, Ymin , Y_max  ������� ���������� ���� GSTgriddedData,
     � ��������������� �����������  ���������� d. ���������� 1-���� ����������, 0- ���� ������.
     */    
    bool haveEqualGrid(GSTgriddedData &d) const;
    /*��������� ��������� ����� ( N_row,  N_col,  X_min,  X_max, Y_min,  Y_max, GSTgrid_type)
    ������� ���������� ���� GSTgriddedData � ���������� �� ���������� d.*/
    
    float qmin() const;
    /*���������� ����������� �������� � �����*/
    
    float qmin(float x0, float x1, float y0, float y1) const;
    /*���������� ����������� �������� �� ��������� �����, ������������*/
    /*������������: float x0, float x1, float y0, float y1*/     
    
    float qmax() const;
    /*���������� ������������ �������� � �����*/
    
    float qmax(float x0, float x1, float y0, float y1) const;
    /*���������� ������������ �������� �� ��������� �����, ������������*/
    /*������������: float x0, float x1, float y0, float y1*/
    
     float qmidl(float x0, float x1, float y0, float y1)const;
    /*���������� ������� �������� �� ��������� �����, ������������*/
    /*������������: float x0, float x1, float y0, float y1*/
    
    void minimize();

    GSTgriddedData &operator=(const GSTgriddedData&);
    
    friend GSTgriddedData operator+(const GSTgriddedData&, const GSTgriddedData&);  
    friend GSTgriddedData operator+(const GSTgriddedData&, float);
    friend GSTgriddedData operator+(float, const GSTgriddedData&);
    friend GSTgriddedData operator-(const GSTgriddedData&, const GSTgriddedData&);  
    friend GSTgriddedData operator-(const GSTgriddedData&, float);
    friend GSTgriddedData operator-(float, const GSTgriddedData&);
    friend GSTgriddedData operator*(const GSTgriddedData&, const GSTgriddedData&);  
    friend GSTgriddedData operator*(const GSTgriddedData&, float);
    friend GSTgriddedData operator*(float, const GSTgriddedData&);
    friend GSTgriddedData operator/(const GSTgriddedData&, const GSTgriddedData&);  
    friend GSTgriddedData operator/(const GSTgriddedData&, float);
    friend GSTgriddedData operator/(float, const GSTgriddedData&);
    
  private:
  
    void resize(int N_row, int N_col);
    /*
    ������������� ������ ����� � NrowxNcol. ��������� ������� ��� ���� ���������.
    ����������� ���������.
    ��������� ���������� ������� � ��� ����� �� ����������. ��� ���������������.
    */
    void resize(int N_row, int N_col, double X_min, double X_max, double Y_min, double Y_max, GSTgrid_type type);
    /*
    ������������� ������ ����� � NrowxNcol. ��������� ������� ��� ���� ���������.
    ����������� ���������.
    ��������������� ����� ��������� ���������� ������� � ��� �����. ��� ��������������.
    */
    void resize(const GSTgriddedData &spec);
    /*
    ��� ��������� ����� ��������������� ������ ��, ��� � ������� spec.
    ���������� ������� ���������, ���������� ������� spec �� ����������.
    */
    
    int Nr;			/*Number of rows - corresponds to X axis, numbered from Xmin to Xmax*/
    int Nc;			/*Number of columns - corresponds to Y axis, numbered from Ymin to Ymax*/
    double Xmin;		/*X minimum value (corresponds to first row for node grid type)*/
    double Xmax;		/*X maximum value (corresponds to last row for node grid type)*/
    double Ymin;		/*Y minimum value (corresponds to first column for node grid type)*/
    double Ymax;		/*Y maximum value (corresponds to last column for node grid type)*/
    GSTgrid_type grid_type;	/*Grid type: values correspond to nodes (node) or averaged over areas (area)*/
    double dX;			/*Interval between neighbouring X values*/
    double dY;			/*Interval between neighbouring Y values*/
    float **V;			/*Data array*/
    char comment[GST_STRLN];	/*Text comment*/
};

class GSTscatterData
/*������ ����� ������ ��� �������� ��������� � ���������� ������� (������� "�����") ������������ � ����� ����� ������������.
� ����� ������ ������������ ���� ����� ������������ ����� ������ ���������� GSTpointMap �.�. � ����������� ������� �� ���������� � ������ ������ ��� ����� ��� ������ ����.
������� ������ �������� ��� ����� ��������. �.�. ��� ������� ���� ���������� ��������� ������� ������. 
*/
{
  public:
  
    GSTscatterData(const char *fields="xyz",float x_const=0, float y_const=0, float z_const=0);
    ~GSTscatterData();
    /*� ������ ��������� ������������ ������ � ������������ ������� ����������� �������� �����.
    (���� '_' �� ����� �������� ������ ����!)
     ����������� ��������� ���������� ���� ����-�� �� ��������� x,y,z �������� ���������� ��� 
     ���� �����.�������� ���� z=3 GSTscatterData("x,y,t",0,0,3) ���������� �������� ���������
     ��������� � ������ ���� ����� (� ����� ������ x,y) �������������� �� �����.
     �������� ������������� ������������ �� ��������� ��� � ��������� ������ ����������.*/
     
    /*���� ����� � ���������� ���� GSTscatterData ���������� ���������� ���� GSTpointMap � ������
    ������� �����, �� ����� ��������� ������ ���� � ������ ������, � � ����������� ���� ����� 
    �������� ������ ��������.*/ 
    
    float x(int i) const;	                	// ���������� ���������� x ����� i
    float y(int i) const;				// ���������� ���������� y ����� i
    float z(int i) const;				// ���������� ���������� z ����� i
    void xyz(int i, float &xv, float &yv, float &zv) const; // � ���������� xv yv zv ����������
                                                        // �������������� ���������� ����� i
    
    void xyz(int i, GSTxyzType &p) const;		        //� ���������� ���� GSTxyzType
                                                         //���������� ���������� ����� i
							
    
    int setX(int i, float xv);				//������������� �������� xv ���������� x 
                                                        //����� i 
    int setY(int i, float yv);				//������������� �������� yv ���������� y 
                                                        //����� i 
    int setZ(int i, float zv);				//������������� �������� zv ���������� z 
                                                        //����� i 
    int setXYZ(int i, float xv, float yv, float zv);	//������������� �������� xv yv zv
                                                        // ����������� ����� i 
							
    int setXYZ(int i, GSTxyzType &p);                   //������������� �������� ���������� �
                                                        // ����������  ����  GSTxyzType  �����������
							//  ����� i  
    
    void allF(int i, GSTpointMap &pm) const;            // ������ ��� ���� ����� i � ���������� 
                                                        // ���� GSTpointMap 
    float f(int i, char key) const;                     // ���������� �������� ���� key ����� i
    
    int setAllF(int i, const GSTpointMap &pm);          // ������������� �������� ���� �����  
                                                        //  �����  i ������ ��� � ���������� ����
							// GSTpointMap.
							
    int setF(int i, char key, float v);                 // ������������� �������� v ���� key 
                                                        // �����    i.
    int setF( char key, float v);                       // ������������� �������� v ���� key
     							// �� ���� ������ GSTscatterData.
    
    int pntDel(int i);                                  // ������� ����� i. �������� ������� ��
                                                        // �������
    
    int pntAdd(GSTpointMap &pm);                        // ��������� ����� � ����� �������. ��������
                                                        // ����� �������� �� ���������� ����
							// GSTpointMap
    
    int numPnt() const;					// ���������� ����� �����
    int numFld() const;					// ���������� ����� �����
    
    void resize(int N, float val=GST_EMPTY_VALUE);      // ������������� ����� ����� � N
    						        // ���� N ������ ��������, �� ��������
						        // ����� ����� ��������������� � val
						   
    void resize(int N, GSTpointMap &pm);                // Not implemented yet
    void reserve(int N); 	                         // ���������� ������ ���  ����� (�����
                                                        // ������ ��� �����������)
    
    int fieldIns(char key, float val=GST_EMPTY_VALUE); // �������� ����  key, ������ ���������
                                                        //- ��������, ������� ����� ��������� 
							//����� ���� (��� ���� �����), ��� �������
							//�� ��������� ����� ���������� ������
							// ��������.
							    
    int fieldDel(char key);				// ������� ����  key
    
    GSTdataVector *field(char key);                     //���������� ��������� �� ������
                                                        // GSTdataVector ���� key	
						        
    float *field_plain(char key);                       // ���������� ��������� ������ ������� ����
                                                        // float  ���� key
    int  namesFLD(std::string &name) const;	        //���������� ����� ����� GSTscatterData 
                                 			// � ���������� name, ���������� �����
							// �����.
							
/*����������� ������������(�����������) �������� ���� X(Y,Z), ���� ���������� n ��������, � ��� ������������ ����� �����, �������  ����������� ������������(�����������) �������� ���� */	
							
    float maxX(int *n=NULL);                           
    float minX(int *n=NULL);
    
    float maxY(int *n=NULL);
    float minY(int *n=NULL);
    
    float maxZ(int *n=NULL);
    float minZ(int *n=NULL);
  
/*���������� ����������� (������������) �������� ���� key , ���� ���������� n ��������, � ��� ������������ ����� �����, �������  ����������� ������������(�����������) �������� ���� */		
    float minF(char key, int *n=NULL);
    float maxF(char key,int *n=NULL);
   
 /*� ���������� p ����������� ����������� (������������) �������� ����� GSTscatterData, ����
  ���������� n ��������, � ��� ������������ ������ �����, ������� ������������  �����������
  (������������)�������� � �����.
  ���  GSTpointNum, ��� �� ��� � GSTpointMap ������ �� 2-� ����� ������ ���� �har-��� ����, ������
  ���� int- ����� �����   map<char,int> GSTpointNum */  
  
    void min(GSTpointMap &p, GSTpointNum *n=NULL);
    void max(GSTpointMap &p, GSTpointNum *n=NULL);					
    
    /*���������� �� ������ � ������� ��������� � �� ��������������� �������� */    
    void sort(int Nx, int Ny);  // ��������� ������� �� ����� Nx Ny -����� ������ �� x, y. 
 
    void sort(float Lx, float Ly, int* nx, int *ny); // ��������� ������� �� ����� Lx, Ly ������ 
                                                  // ����� ��  x,y, �  nx,ny- ������������ �������
						  //  ���������� ������ �� x � �� y/ 
 
    int np_block(int bi, int bj); 
    /* ���������� ����� ����� ������������� �����  bi, bj
      ����� ����������� ������� ��� ��������!*/
                 
     void  N_block( float x, float y, int &bi, int &bj);      
     /* � bi  bj ������������ ����� ����� �� x,y �������������.
     ���� ����� ������� �� ������� ����� ������, �� ����� ����� ������������� ��� ������, ��� 
     ���� ������. �.�. ���� ����� ����� ��� ��� ����� ���������� � ���������� ���������.       
							 */ 
     bool BlockExist(int bi, int bj);
     /* ���������� 1, ����  ���� � ����� ������� ����������
       0- ���� ������ ����� ���.*/
       
     bool BlockSort();
     /*���������� 1 ���� ������ ������������, 0- ���� ���������� �� ���� �����������*/
     
     void BlockNumber(int &nBx, int &nBy);
     /*���������� ����� ������ �� x � �� y � nBx nBy �������������� */
     
     void Block_step(float &dx, float &dy);
     /* � dx dy ����� �������� ��� (������ �����) �� x y*/
     
     int bxyz(int bi, int bj, int i, float &xv, float &yv, float &zv);
     /*������ ������� void xyz(int i, float &xv, float &yv, float &zv).
     i ����� ����� ������  �����, ��������� 0� 0 �� N-1 , N-����� ����� � �����.
     bi bj- ����� �����.
     ���������� ����� ����� � ����� ���������.     
      
     ��������!!! ����������� �������. �� �������� �������� ����������������� �������,
     ������������ ������ ������ ����� � ������������ ������ ����������� ������ ����� 
     � ������ ����� */
     
     
     float bf(int bi, int bj, int i, char key); 
 
     /* ������ ������� float f(int i, char key);  ���������� �������� ���� key �����  i,
     ��� i ����� ����� ������  �����, ��������� 0� 0 �� N-1 , N-����� ����� � �����.
     bi bj- ����� �����*/
     
     float  lbf(int bi, int bj, int i, char key); 
     /*��������! ����������� (�.�. ������ ��� �������� ������������� ������������ �� 
     ����������) ������� ���������� �������  float bf(int bi, int bj, int i, char key); */
     
    
    void allBf(int bi, int bj, int i, GSTpointMap &pm);
  /* ������ �������  void allF(int i, GSTpointMap &pm) ������ ��� ���� ����� i (����������
   ����� ����� � �����) � ���������� pm, bi bj- ����� �����*/  
   						
							
    //The following function is for debug purposes only
    void chknp();
    
  protected:
    
    int Nf_;
    int Np_;
    bool xvar_, yvar_, zvar_;
    float const_x_, const_y_, const_z_;
    std::vector<char> keys_;
    GSTscatterMap V;
    
    GSTpointNum _Nmin ; // ������ ����� � ����������
    GSTpointNum _Nmax ;
    
    int **_Num; // � ������ ��������� ������� _Num[i][j] �������� ����� ���������� ���������
                // �������������� ����� [i][j]  
    //float* _X; //������ � ��������� ������
   // float* _Y;
    float _dx;
    float _dy;
    int _Br; //����� ������ �� x
    int _Bc; //����� ������ �� y
    bool  _sort;
};
 
 

class GSTxyzData
{
  public:
  
    GSTxyzData();
    GSTxyzData(int N);
    GSTxyzType pnt(int i);
    void pntSet(int i, GSTxyzType &p);
    void pntSet(int i, float x, float y, float z);
    void pntSetZ(int i, float z);
    int length();
    std::string getComment(std::string &comm);
    void setComment(std::string &comm);
    void resize(int N);
    void append(GSTxyzData *tail);
    void replace(GSTxyzData *repl);
          
  private:
  
    int Np;
    std::vector<GSTxyzType> V;
    std::string comment;
};

#endif
