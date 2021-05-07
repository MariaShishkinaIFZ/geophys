/*Subroutines for memory allocation*/

void *fmemalloc(size_t nunits, size_t unitsz);
void fmemfree(void *ptr);
void alloc_float_matrix(float ***Buf, size_t nstr, size_t ncol);
void free_float_matrix(float **Buf, size_t nstr, size_t ncol);
void alloc_uns_matrix(unsigned ***Buf, size_t nstr, size_t ncol);
void free_uns_matrix(unsigned **Buf, size_t nstr, size_t ncol);
void alloc_int_matrix(int ***Buf, size_t nstr, size_t ncol);
void free_int_matrix(int **Buf, size_t nstr, size_t ncol);
void alloc_float_3d(float ****Buf, size_t nlayers, size_t n2dstr, size_t n2dcol);
void free_float_3d(float ***Buf, size_t nlayers, size_t n2dstr, size_t n2dcol);
