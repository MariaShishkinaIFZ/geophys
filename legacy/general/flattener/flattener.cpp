#include <stdio.h>
#include <math.h>
#include <stdlib.h>
float *a, *b, *x, *y, **sums;
int N, K;
//N - number of data points
//K - polinom power
//K<=N
char filename[256];
FILE* InFile=NULL;

void count_num_lines(){
   //count number of lines in input file - number of equations
   int nelf=0;       //non empty line flag
   do{
       nelf = 0;
       while(fgetc(InFile)!='\n' && !feof(InFile)) nelf=1;
       if(nelf) N++;
   }while(!feof(InFile));
}

void freematrix(){
   //free memory for matrixes
   int i;
   for(i=0; i<K+1; i++){
       delete [] sums[i];
   }
   delete [] a;
   delete [] b;
   delete [] x;
   delete [] y;
   delete [] sums;
}

void allocmatrix(){
   //allocate memory for matrixes
   int i,j,k;
   a = new float[K+1];
   b = new float[K+1];
   x = new float[N];
   y = new float[N];
   sums = new float*[K+1];
   if(x==NULL || y==NULL || a==NULL || sums==NULL){
       printf("\nNot enough memory to allocate. N=%d, K=%d\n", N, K);
       exit(-1);
   }
   for(i=0; i<K+1; i++){
       sums[i] = new float[K+1];
       if(sums[i]==NULL){
	   printf("\nNot enough memory to allocate for %d equations.\n", K+1);
       }
   }
   for(i=0; i<K+1; i++){
       a[i]=0;
       b[i]=0;
       for(j=0; j<K+1; j++){
	   sums[i][j] = 0;
       }
   }
   for(k=0; k<N; k++){
       x[k]=0;
       y[k]=0;
   }
}

void readmatrix(){
   int i=0,j=0, k=0;
   //read x, y matrixes from input file
   for(k=0; k<N; k++){
       fscanf(InFile, "%f", &x[k]);
       fscanf(InFile, "%f", &y[k]);
   }
   //init square sums matrix
   for(i=0; i<K+1; i++){
       for(j=0; j<K+1; j++){
	   sums[i][j] = 0;
	   for(k=0; k<N; k++){
	       sums[i][j] += pow(x[k], i+j);
	   }
       }
   }
   //init free coefficients column
   for(i=0; i<K+1; i++){
       for(k=0; k<N; k++){
	   b[i] += pow(x[k], i) * y[k];
       }
   }
}

void printresult(){
   //print polynom parameters
   int i=0;
   for(i=0; i<K+1; i++){
       fprintf(stdout,"%f ", a[i]);
   }
    fprintf(stdout,"\n ");
}

/*void testsolve(){
   //test that ax=y
   FILE* OutFile;
   int i=0, j=0;
   float sum, s, begin, end, step;
   begin = x[0][1];
   end = x[N-1][1];
   step = (end-begin)/20;
   printf("\n");
   for(i=0; i<N; i++){
       sum = 0;
       for(j=0; j<N; j++){
	   sum += x[i][j]*a[j];
       }
       printf("%f\t%f\n", s, y[i]);
   }
   OutFile = fopen("test.txt", "wt");
   for(s = begin; s<=end; s+=step){
       sum = 0;
       for(j=0; j<N; j++){
	   sum += pow(s, j)*a[j];
       }
       fprintf(OutFile, "%f\t%f\n", s, sum);
   }
   fclose(OutFile);
}*/
/*void printresult(){
   int i=0;
   printf("\n");
   printf("Result\n");
   for(i=0; i<N; i++){
       printf("A%d = %f\n", i, a[i]);
   }
}*/
void diagonal(){
   int i, j, k;
   float temp=0;
   for(i=0; i<K+1; i++){
       if(sums[i][i]==0){
	   for(j=0; j<K+1; j++){
	       if(j==i) continue;
	       if(sums[j][i] !=0 && sums[i][j]!=0){
		   for(k=0; k<K+1; k++){
		       temp = sums[j][k];
		       sums[j][k] = sums[i][k];
		       sums[i][k] = temp;
		   }
		   temp = b[j];
		   b[j] = b[i];
		   b[i] = temp;
		   break;
	       }
	   }
       }
   }
}

void cls(){
   for(int i=0; i<25; i++) printf("\n");
}

void errquit(const char *s)
{
  puts(s); exit(EXIT_FAILURE);
}

int main (int narg, char **argv)
{
   int i=0,j=0, k=0;
   
  if (narg < 3) errquit("Use: curvefitter <data_file> <power of approximation polinom K<N> .\n");
  
   InFile=fopen(argv[1],"rt");
   count_num_lines();
//    fprintf(stderr,"Number of points: N=%d.\n", N);
  
   sscanf(argv[2],"%d", &K);

   

   allocmatrix();
   rewind(InFile);
   //read data from file
   readmatrix();
   //check if there are 0 on main diagonal and exchange rows in that case
   diagonal();
   fclose(InFile);
   //printmatrix();
   //process rows
   for(k=0; k<K+1; k++){
       for(i=k+1; i<K+1; i++){
	   if(sums[k][k]==0){
	       printf("\nSolution is not exist.\n");
	   }
	   float M = sums[i][k] / sums[k][k];
	   for(j=k; j<K+1; j++){
	       sums[i][j] -= M * sums[k][j];
	   }
	   b[i] -= M*b[k];
       }
   }
   //printmatrix();
   for(i=(K+1)-1; i>=0; i--){
       float s = 0;
       for(j = i; j<K+1; j++){
	   s = s + sums[i][j]*a[j];
       }
       a[i] = (b[i] - s) / sums[i][i];
   }

   //InFile = fopen(filename, "rt");
   //readmatrix();
   //fclose(InFile);
   //printmatrix();
   //testsolve();
   printresult();
   freematrix();
}
