#include <cstdio>
#include <cstdlib>
#include "submem.h"

void submem_error(const char *message)
{
  fputs(message,stderr);
  exit(EXIT_FAILURE);
}

void *fmemalloc(size_t nunits, size_t unitsz)
{
  void *ptr;
  if (!(ptr=calloc(nunits,unitsz)))
  {
	 printf("Insufficient memory. Program terminated.\n");
	 perror("System error ");
	 exit(EXIT_FAILURE);
  }
  return ptr;
}

void fmemfree(void *ptr)
{
  free(ptr);
}

void alloc_float_matrix(float ***Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  *Buf=(float**)fmemalloc(nstr,sizeof(float*));
  for (i=0;i<nstr;i++) (*Buf)[i]=(float*)fmemalloc(ncol,sizeof(float));
}

void free_float_matrix(float **Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  for (i=0;i<nstr;i++) free(Buf[i]);
  free(Buf);
}

void alloc_double_matrix(double ***Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  *Buf=(double**)fmemalloc(nstr,sizeof(double*));
  if (!Buf) submem_error("Memory allocation error at alloc_double_matrix [1]\n");
  for (i=0;i<nstr;i++) 
    if (!((*Buf)[i]=(double*)fmemalloc(ncol,sizeof(double))))
      submem_error("Memory allocation error at alloc_double_matrix [2]\n");
}

void free_double_matrix(double **Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  for (i=0;i<nstr;i++) free(Buf[i]);
  free(Buf);
}

void alloc_uns_matrix(unsigned ***Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  *Buf=(unsigned**)fmemalloc(nstr,sizeof(unsigned*));
  for (i=0;i<nstr;i++) (*Buf)[i]=(unsigned*)fmemalloc(ncol,sizeof(unsigned));
}

void free_uns_matrix(unsigned **Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  for (i=0;i<nstr;i++) free(Buf[i]);
  free(Buf);
}

void alloc_int_matrix(int ***Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  *Buf=(int**)fmemalloc(nstr,sizeof(int*));
  for (i=0;i<nstr;i++) (*Buf)[i]=(int*)fmemalloc(ncol,sizeof(int));
}

void free_int_matrix(int **Buf, size_t nstr, size_t ncol)
{
  unsigned i;
  for (i=0;i<nstr;i++) free(Buf[i]);
  free(Buf);
}

void alloc_float_3d(float ****Buf, size_t nlayers, size_t n2dstr, size_t n2dcol)
{
  unsigned k,j;
  *Buf=(float***)fmemalloc(nlayers,sizeof(float**));
  for (k=0;k<nlayers;k++) 
  {
    (*Buf)[k]=(float**)fmemalloc(n2dstr,sizeof(float*));
    for (j=0;j<n2dstr;j++)
      (*Buf)[k][j]=(float*)fmemalloc(n2dcol,sizeof(float));
  }
}

void free_float_3d(float ***Buf, size_t nlayers, size_t n2dstr, size_t n2dcol)
{
  unsigned k,j;
  for (k=0;k<nlayers;k++)
  {
    for (j=0;j<n2dstr;j++)
      fmemfree(Buf[k][j]);
    fmemfree(Buf[k]);
  }
  fmemfree(Buf);
}
