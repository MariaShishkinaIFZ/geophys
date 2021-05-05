#ifndef _BICSB_H
#define _BICSB_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>		// for std:accumulate()
#include <limits>		// C++ style numeric_limits<T>
#include "csc.h"
#include "mortoncompare.h"

using namespace std;

// CSB variant where nonzeros "within each block" are sorted w.r.t. the bit-interleaved order
// Implementer's (Aydin) notes:
//	- to ensure correctness in BlockPar, we use square blocks (lowcolmask = highcolmask)
template <class T, class ITYPE>
class BiCsb
{
public:
	BiCsb ():nz(0), m(0), n(0), nbc(0), nbr(0) {}	// default constructor (dummy)

	BiCsb (ITYPE size,ITYPE rows, ITYPE cols, int workers);
	BiCsb (ITYPE size,ITYPE rows, ITYPE cols, ITYPE * ri, ITYPE * ci, T * val, int workers, ITYPE forcelogbeta = 0);

	BiCsb (const BiCsb<T, ITYPE> & rhs);			// copy constructor
	~BiCsb();
	BiCsb<T,ITYPE> & operator=(const BiCsb<T,ITYPE> & rhs);	// assignment operator
	BiCsb (Csc<T, ITYPE> & csc, int workers);
	
	ofstream & PrintStats(ofstream & outfile) const;
	ITYPE colsize() const { return n;} 
	ITYPE rowsize() const { return m;} 
	bool isPar() const { return ispar; }
	float RowImbalance() const;
        float ColImbalance() const;

private:
	void Init(int workers, ITYPE forcelogbeta = 0);
	void SubSpMV(ITYPE * btop, ITYPE bstart, ITYPE bend, const T * x, T * suby) const;
	void SubSpMVTrans(ITYPE col, ITYPE rowstart, ITYPE rowend, const T * __restrict x, T * __restrict suby) const;
	void BMult(ITYPE ** chunks, ITYPE start, ITYPE end, const T * x, T * y, ITYPE ysize) const;
	void BTransMult(ITYPE col, ITYPE rowstart, ITYPE rowend, const T * x, T * y, ITYPE ysize) const;
	void BlockPar(ITYPE start, ITYPE end, const T * __restrict subx, T * __restrict suby, 
					ITYPE rangebeg, ITYPE rangeend, ITYPE cutoff) const;
	void BlockParT(ITYPE start, ITYPE end, const T * __restrict subx, T * __restrict suby, 
					ITYPE rangebeg, ITYPE rangeend, ITYPE cutoff) const;
	void SortBlocks(pair<ITYPE, pair<ITYPE,ITYPE> > * pairarray, T * val);

	ITYPE ** top ;	// pointers array (indexed by higher-order bits of the coordinate index), size ~= ntop+1
	ITYPE * bot;	// contains lower-order bits of the coordinate index, size nnz 
	T * num;	// contains numerical values, size nnz

	bool ispar;
	ITYPE nz;		// # nonzeros
	ITYPE m;		// # rows
	ITYPE n;		// # columns
	ITYPE blcrange;		// range indexed by one block

	ITYPE nbc;		// #{column blocks} = #{blocks in any block row}
	ITYPE nbr; 		// #{block rows)
	
	ITYPE rowlowbits;	// # lower order bits for rows
	ITYPE rowhighbits;
	ITYPE highrowmask;  	// mask with the first log(m)/2 bits = 1 and the other bits = 0  
	ITYPE lowrowmask;

	ITYPE collowbits;	// # lower order bits for columns
	ITYPE colhighbits;
	ITYPE highcolmask;  // mask with the first log(n)/2 bits = 1 and the other bits = 0  
	ITYPE lowcolmask;

	MortonCompare<ITYPE> mortoncmp;	// comparison operator w.r.t. the N-morton layout

	template <typename U, typename UTYPE>
	friend void bicsb_gaxpy (const BiCsb<U, UTYPE> & A, const U * x, U * y);
	
	template <typename U, typename UTYPE>
	friend void bicsb_gaxpy_trans (const BiCsb<U, UTYPE> & A, const U * x, U * y);
};

// Operation y = A*x+y 
template <typename T, typename ITYPE>
void bicsb_gaxpy (const BiCsb<T, ITYPE> & A, const T * x, T * y)
{
	ITYPE ysize = A.lowrowmask + 1;			// size of the output subarray (per block row - except the last)

	if( A.isPar() && ( A.RowImbalance() > 2) )
	{	
		cilk_for (ITYPE i = 0 ; i < A.nbr ; ++i)	// for all unredundant block rows of A 
		{
			ITYPE *  btop = A.top [i];		// get the pointer to this block row
			ITYPE rhi = ((i << A.rowlowbits) & A.highrowmask);
			T * suby = &y[rhi];

			ITYPE thsh = BREAKEVEN * ysize;
			vector<ITYPE*> chunks;
			chunks.push_back(btop);
			for(ITYPE j =0; j < A.nbc; )
			{
				ITYPE count = btop[j+1] - btop[j];
				if(count < thsh && j < A.nbc)
				{
					while(count < thsh && j < A.nbc)
					{
						count += btop[(++j)+1] - btop[j]; 
					}
					chunks.push_back(btop+j);	// push, but exclude the block that caused the overflow
				}
				else
				{
					chunks.push_back(btop+(++j));	// don't exclude the overflow block if it is the only block in that chunk
				}
			}
			// In std:vector, the elements are stored contiguously so that we can 
			// treat &chunks[0] as an array of pointers to ITYPE w/out literally copying it to ITYPE**
			if(i==(A.nbr-1))	// last iteration
			{
				A.BMult(&chunks[0], 0, chunks.size()-1, x, suby,  A.rowsize() - ysize*i);
			}
			else
			{
				A.BMult(&chunks[0], 0, chunks.size()-1, x, suby, ysize);
			}
		}
	}
	else
	{
		cilk_for (ITYPE i = 0 ; i < A.nbr ; ++i)    // for all unredundant block rows of A 
                {
			ITYPE * btop = A.top [i];                       // get the pointer to this block row
                        ITYPE rhi = ((i << A.rowlowbits) & A.highrowmask);
                        T * suby = &y[rhi];

			A.SubSpMV(btop, 0, A.nbc, x, suby);
		}
	}
}

// Operation y = (A^t)*x+y 
template <typename T, typename ITYPE>
void bicsb_gaxpy_trans (const BiCsb<T, ITYPE> & A, const T * x, T * y)
{
	ITYPE ysize = A.lowcolmask + 1;			// size of the output subarray (per block column - except the last)

	if(A.isPar() && ( A.ColImbalance() > 2 ))
        {
		cilk_for (ITYPE j = 0 ; j < A.nbc; ++j)	// for all unredundant block columns of A 
		{
			ITYPE rhi = ((j << A.collowbits) & A.highcolmask);
			T * suby = &y[rhi];

			if(j==(A.nbc-1))	// last iteration
			{
				A.BTransMult(j, 0, A.nbr, x, suby, A.colsize() - ysize*j);
			}
			else
			{
				A.BTransMult(j, 0, A.nbr, x, suby, ysize);
			}
		}
	}
	else
	{
		cilk_for (ITYPE j =0; j< A.nbc; ++j)
		{
			ITYPE rhi = ((j << A.collowbits) & A.highcolmask);
                        T * suby = &y[rhi];

			A.SubSpMVTrans(j, 0, A.nbr, x, suby);
		}
    	}
}

#include "bicsb.cpp"	
#endif
