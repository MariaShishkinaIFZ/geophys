#ifndef _SPVEC_H_
#define _SPVEC_H_

#include "csc.h"
#include "bicsb.h"
#include "matmul.h"

template <class T, class ITYPE>
class Spvec
{
public:
	Spvec (): n(0) {};
	Spvec (ITYPE dim);				
	Spvec (T * darr, ITYPE dim);
	Spvec (const Spvec<T,ITYPE> & rhs);		
	~Spvec();
	Spvec<T,ITYPE> & operator=(const Spvec<T, ITYPE> & rhs);	

	T& operator[] (const ITYPE nIndex)
	{
		return arr[nIndex];
	}

	//! Delayed evaluations using compositors for SpMV operation...  y <- y + Ax
	Spvec<T,ITYPE> & operator+=(const Matmul< Csc<T, ITYPE>, Spvec<T,ITYPE> > & matmul);	
	Spvec<T,ITYPE> & operator+=(const Matmul< BiCsb<T, ITYPE>, Spvec<T,ITYPE> > & matmul);

	void fillzero();
	void fillrandom();

	ITYPE size() const { return n;}
	T * getarr(){ return arr;} 

private:
	T * arr;
	ITYPE n;
};

#include "spvec.cpp"
#endif

