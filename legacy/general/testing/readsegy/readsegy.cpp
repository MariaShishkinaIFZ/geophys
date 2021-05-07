#include <vector>
#include <cstdio>

using namespace std;

typedef vector<float> fvec;
typedef vector<fvec> ftable;
typedef struct
	{
	float x,y,z;
	}
	triple;
typedef vector<triple> coords_type;

typedef struct
	{
	  float xs, ys, zs;		// source coordinates
	  float dt;			// sample interval
	  int Ns;			// Number of samples per data trace
	  int Ntraces;			// Number of traces in file

	  ftable traces;		// matrix of traces
	  fvec t0;			// time delay between source explosion and start of recording
	  coords_type rec_coords;	// vector of reciver's coordinates
	  coords_type sou_coords;	// vector of reciver's coordinates
	}
	CSG_data_type;

	
	
int read_segy(CSG_data_type &data, const char *filename)
{
  
   
   void mcread(char*  dummy, const int nbytes, const int nitems, FILE* file);
   void msread(short* dummy, const int nbytes, const int nitems, FILE* file);
   void miread(int*   dummy, const int nbytes, const int nitems, FILE* file);
   void mfread(float* dummy, const int nbytes, const int nitems, FILE* file);

  
  int i, j;
  int chan, cdp;
  char dummy[3601];  // dummy buffer for skipping unused header information
  float rec_elev, rec_x, rec_y;
  float sou_x, sou_y;  
  float value;
  FILE* in;
  
  
    if ((in = fopen(filename,"r")) == NULL) {
      printf("\nERROR:  cannot find file %s\n\n",filename);
      exit(1);
      
    }
    
   printf("---chek point---%s\n",filename);
   
   
   mcread(dummy,     3216, 1, in);
   mfread(&data.dt,    2, 1, in);
   mcread(dummy,        2, 1, in);
   miread(&data.Ns,     2, 1, in);
   mcread(dummy,      378, 1, in);
   
   j=0;
   
   while (fread(&dummy, 12, 1, in) == 1) {

      miread(&chan,        4, 1, in); 	//13-16
      mcread(dummy,        4, 1, in);
      miread(&cdp,         4, 1, in);	//21-24
      mcread(dummy,       16, 1, in);
      mfread(&rec_elev,    4, 1, in);	//41-44
      mcread(dummy,       28, 1, in);
      mfread(&sou_x,       4, 1, in);	//73-76
      mfread(&sou_y,       4, 1, in);	//77-80
//    mcread(dummy,       36, 1, in);
      mfread(&rec_x,       4, 1, in);	//81-84
      mfread(&rec_y,       4, 1, in);	//85-88
      mcread(dummy,       20, 1, in);
      printf("---chek point---\n");
      //mfread(&data.t0[j],  2, 1, in);	//109-110
      mfread(&rec_y,       2, 1, in);	//109-110
      mcread(dummy,      130, 1, in);
      
      printf("---chek point---1\n");

      data.rec_coords[j].x=sou_x;
      data.rec_coords[j].y=sou_y;
      data.rec_coords[j].z=0;			///////
      
      data.sou_coords[j].x=rec_x;
      data.sou_coords[j].y=rec_y;
      data.sou_coords[j].z=rec_elev;
      
      printf("---chek point---\n");
      

      for (int i = 0 ; i < data.Ns; i++) {

         mfread(&value, 4, 1, in);


//            if (ieee != 1)
//               value = ibm2ieee(value);

	    data.traces[j][i]=value;
	    
//            if      (xmode == 0)
//               fprintf(out,"%4d   %8.*f %12.4f\n",chan,           prec, time, value);
//            else if (xmode == 1)
//               fprintf(out,"%4d   %8.*f %12.4f\n",cdp,            prec, time, value);
//            else if (xmode == 2)
//               fprintf(out,"%8.2f %8.*f %12.4f\n",disttable[cdp], prec, time, value);
         

      }
      
      j++;
   }
   
   data.Ntraces=j;

   fclose(in);
   

}


int main(int narg, char **argv)
{  
  int i;
  FILE* fo;
  CSG_data_type shot;

  read_segy(shot,argv[1]);

  fo=fopen(argv[2],"wt");
  
  for (i=0;i<shot.Ns;i++)
    fprintf(fo,"%g\n",shot.traces[shot.Ntraces][i]);
}


//////////////////////////////////////////////////////////////
//
// SUBROUTINES
//
//////////////////////////////////////////////////////////////

void mcread(char* dummy, const int nbytes, const int nitems, FILE* file) {
   if (fread (dummy, nbytes, nitems, file) == 0) {
      printf("\nERROR:  unexpected end of SEGY file\n");
      exit(1);
   }
}
void msread(short* dummy, const int nbytes, const int nitems, FILE* file) {
   if (fread (dummy, nbytes, nitems, file) == 0) {
      printf("\nERROR:  unexpected end of SEGY file\n");
      exit(1);
   }
}
void miread(int* dummy, const int nbytes, const int nitems, FILE* file) {
   if (fread (dummy, nbytes, nitems, file) == 0) {
      printf("\nERROR:  unexpected end of SEGY file\n");
      exit(1);
   }
}
void mfread(float* dummy, const int nbytes, const int nitems, FILE* file) {
   if (fread (dummy, nbytes, nitems, file) == 0) {
      printf("\nERROR:  unexpected end of SEGY file\n");
      exit(1);
   }
}

