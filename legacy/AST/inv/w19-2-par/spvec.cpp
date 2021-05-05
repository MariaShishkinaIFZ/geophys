#include "csc.h"
#include "spvec.h"
#include "utility.h"
#include "randgen.h"
#include <cassert>

// constructor that generates a junk dense vector 
template <class T, class ITYPE>
Spvec<T,ITYPE>::Spvec (ITYPE dim): n(dim)
{
	assert(n != 0);
	arr = new T[dim];
}

template <class T, class ITYPE>
Spvec<T,ITYPE>::Spvec (T * darr, ITYPE dim): n(dim)
{
	assert(n != 0);
	arr = new T[dim];

	for(ITYPE i=0; i< n; ++i)
	{
		arr[i] = darr[i];
	}
}


// copy constructor
template <class T, class ITYPE>
Spvec<T,ITYPE>::Spvec (const Spvec<T, ITYPE> & rhs): n(rhs.n) 
{
	if(n > 0)
	{
		arr = new T[n];

		for(ITYPE i=0; i< n; ++i)		
			arr[i]= rhs.arr[i];
	}
}

template <class T, class ITYPE>
Spvec<T,ITYPE> & Spvec<T,ITYPE>::operator= (const Spvec<T,ITYPE> & rhs)
{
	if(this != &rhs)		
	{
		if(n > 0)
		{
			delete [] arr;
		}

		n	= rhs.n;
		if(n > 0)
		{
			arr = new T[n];
			for(ITYPE i=0; i< n; ++i)		
				arr[i]= rhs.arr[i];
		}
	}
	return *this;
}


template <class T, class ITYPE>
Spvec<T,ITYPE>::~Spvec()
{
	if ( n > 0)
	{
		delete [] arr;
	}
}

template <class T, class ITYPE>
Spvec<T,ITYPE> & Spvec<T,ITYPE>::operator+=(const Matmul< Csc<T, ITYPE>, Spvec<T,ITYPE> > & matmul)
{
	if((n == matmul.op1.rowsize()) && (matmul.op1.colsize() == matmul.op2.size()))		// check compliance
	{
		csc_gaxpy(matmul.op1, const_cast< T * >(matmul.op2.arr), arr); 
	}
	else
	{
		cout<< "Detected noncompliant matvec..." << endl;
	}
	return *this;
}

template <class T, class ITYPE>
Spvec<T,ITYPE> & Spvec<T,ITYPE>::operator+=(const Matmul< BiCsb<T, ITYPE>, Spvec<T,ITYPE> > & matmul)
{
	if((n == matmul.op1.rowsize()) && (matmul.op1.colsize() == matmul.op2.size()))		// check compliance
	{
		bicsb_gaxpy(matmul.op1, matmul.op2.arr, arr); 
	}
	else
	{
		cout<< "Detected noncompliant matvec..." << endl;
	}
	return *this;
}

// populate the vector with random entries
// currently, only works for T "double" 
template <class T, class ITYPE>
void Spvec<T,ITYPE>::fillrandom()
{
	RandGen G;

	for(ITYPE i=0; i< n; ++i)
	{
		arr[i] = G.RandReal();
	}
}

// populate the vector with zeros
template <class T, class ITYPE>
void Spvec<T,ITYPE>::fillzero()
{
	for(ITYPE i=0; i< n; ++i)
	{
		arr[i] = 0;
	}
}

