#ifndef _CSC_H_
#define _CSC_H_

#include "triple.h"


template <class T, class ITYPE>
struct Triple;

template <class T, class ITYPE>
class Csc
{
public:
	Csc ():nz(0), m(0), n(0) {}				// default constructor
	Csc (ITYPE size,ITYPE rows, ITYPE cols);
	Csc (const Csc<T, ITYPE> & rhs);		// copy constructor
	~Csc();
	Csc<T, ITYPE> & operator=(const Csc<T, ITYPE> & rhs);	// assignment operator
	Csc (Triple<T, ITYPE> * triples, ITYPE size, ITYPE rows, ITYPE cols);
	Csc (ITYPE * ri, ITYPE * ci, T * val, ITYPE size, ITYPE rows, ITYPE cols);

	ITYPE colsize() const { return n;} 
	ITYPE rowsize() const { return m;} 
	ITYPE * getjc() const { return jc;} 
	ITYPE * getir() const { return ir;} 
	T * getnum() const { return num;} 

private:
	void Resize(ITYPE nsize);

	ITYPE * jc ;	// col pointers, size n+1
	ITYPE * ir;		// row indices, size nnz 
	T * num;	// numerical values, size nnz 
	
	ITYPE nz;
	ITYPE m;		//  number of rows
	ITYPE n;		//  number of columns

	template <class U, class UTYPE>
	friend class Sym;
	template <class U, class UTYPE>
	friend class BiCsb;

	template <typename U, typename UTYPE>
	friend void csc_gaxpy (const Csc<U, UTYPE> & A, U * x, U * y);

	template <typename U, typename UTYPE>
	friend void csc_gaxpy_trans (const Csc<U, UTYPE> & A, U * x, U * y);
};

/* y = A*x+y */
template <typename T, typename ITYPE>
void csc_gaxpy (const Csc<T, ITYPE> & A, T * x, T * y)
{
    for (ITYPE j = 0 ; j < A.n ; ++j)	// for all columns of A
    {
	for (ITYPE k = A.jc [j] ; k < A.jc [j+1] ; ++k)	// scale jth column with x[j]
	{
		y [ A.ir[k] ] += A.num[k]  * x [j] ;
	}
    }
}

/* y = A' x + y */
template <typename T, typename ITYPE>
void csc_gaxpy_trans(const Csc<T,ITYPE> & A, T * x, T * y)
{
	for (ITYPE j = 0; j< A.n; ++j)
	{
		for(ITYPE k= A.jc[j]; k < A.jc[j+1]; ++k)
		{
			y[j] += A.num[k] * x [ A.ir[k] ]; 
		}
	}
}

#include "csc.cpp"	// Template member function definitions need to be known to the compiler
#endif

