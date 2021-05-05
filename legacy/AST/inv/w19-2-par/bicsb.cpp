#include "bicsb.h"
#include "utility.h"
#include <iostream>
#include <cassert>

// Choose block size as big as possible given the following constraints
// 1) The bot array is addressible by ITYPE
// 2) The parts of x & y vectors that a block touches fits into L2 cache
// 3) There's enough parallel slackness for block rows (at least SLACKNESS * CILK_NPROC)

template <class T, class ITYPE>
void BiCsb<T, ITYPE>::Init(int workers, ITYPE forcelogbeta)
{
	ispar = (workers > 1);
	ITYPE roundrowup = nextpoweroftwo(m);
	ITYPE roundcolup = nextpoweroftwo(n);

	// if indices are negative, highestbitset returns -1, 
	// but that will be caught by the sizereq below
	ITYPE rowbits = highestbitset(roundrowup);
	ITYPE colbits = highestbitset(roundcolup);
	bool sizereq;
	if (ispar)
	{
		sizereq = ((IntPower<ITYPE>(2,rowbits) > SLACKNESS * workers) 
			&& (IntPower<ITYPE>(2,colbits) > SLACKNESS * workers));
	}
	else
	{
		sizereq = ((rowbits > 1) && (colbits > 1));
	}

	if(!sizereq)
	{
		cerr << "Matrix too small for this library" << endl;
		return;
	}

	rowlowbits = rowbits-1;	
	collowbits = colbits-1;	
	ITYPE inf = numeric_limits<ITYPE>::max();
	ITYPE maxbits = highestbitset(inf);

	rowhighbits = rowbits-rowlowbits;	// # higher order bits for rows (has at least one bit)
	colhighbits = colbits-collowbits;	// # higher order bits for cols (has at least one bit)
	if(ispar)
	{
		while(IntPower<ITYPE>(2,rowhighbits) < SLACKNESS * workers)
		{
			rowhighbits++;
			rowlowbits--;
		}
	}

	// calculate the space that suby occupies in L2 cache
	ITYPE yL2 = IntPower<ITYPE>(2,rowlowbits) * sizeof(T);
	while(yL2 > L2SIZE)
	{
		yL2 /= 2;
		rowhighbits++;
		rowlowbits--;
	}

	// calculate the space that subx occupies in L2 cache
	ITYPE xL2 = IntPower<ITYPE>(2,collowbits) * sizeof(T);
	while(xL2 > L2SIZE)
	{
		xL2 /= 2;
		colhighbits++;
		collowbits--;
	}
	
	// blocks need to be square for correctness (maybe generalize this later?) 
	while(rowlowbits+collowbits > maxbits)
	{
		if(rowlowbits > collowbits)
		{
			rowhighbits++;
			rowlowbits--;
		}
		else
		{
			colhighbits++;
			collowbits--;
		}
	}
	while(rowlowbits > collowbits)
	{
		rowhighbits++;
		rowlowbits--;
	}
	while(rowlowbits < collowbits)
	{
		colhighbits++;
		collowbits--;
	}
	assert (collowbits == rowlowbits);

	lowrowmask = IntPower<ITYPE>(2,rowlowbits) - 1;
	lowcolmask = IntPower<ITYPE>(2,collowbits) - 1;
	if(forcelogbeta != 0)
	{
		ITYPE candlowmask  = IntPower<ITYPE>(2, forcelogbeta) -1;
		cout << "Forcing beta to "<< (candlowmask+1) << " instead of the chosen " << (lowrowmask+1) << endl;
		cout << "Warning : No checks are performed on the beta you have forced, anything can happen !" << endl;
		lowrowmask = lowcolmask = candlowmask;
		rowlowbits = collowbits = forcelogbeta;
		rowhighbits = rowbits-rowlowbits; 
		colhighbits = colbits-collowbits; 
	}
	else 
	{
		double sqrtn = sqrt(sqrt(static_cast<double>(m) * static_cast<double>(n)));
		ITYPE logbeta = static_cast<ITYPE>(ceil(log2(sqrtn))) + 2;
		if(rowlowbits > logbeta)
		{
			rowlowbits = collowbits = logbeta;
			lowrowmask = lowcolmask = IntPower<ITYPE>(2, logbeta) -1;
			rowhighbits = rowbits-rowlowbits;
	                colhighbits = colbits-collowbits;
		}
		cout << "Beta chosen to be "<< (lowrowmask+1) << endl;
	}
	highrowmask = ((roundrowup - 1) ^ lowrowmask);
	highcolmask = ((roundcolup - 1) ^ lowcolmask);
	
	// nbc = #{block columns} = #{blocks in any block row},  nbr = #{block rows)
	ITYPE blcdimrow = lowrowmask + 1;
        ITYPE blcdimcol = lowcolmask + 1;
        nbr = static_cast<ITYPE>(ceil(static_cast<double>(m) / static_cast<double>(blcdimrow)));
        nbc = static_cast<ITYPE>(ceil(static_cast<double>(n) / static_cast<double>(blcdimcol)));
	
	blcrange = (lowrowmask+1) * (lowcolmask+1);	// range indexed by one block
	mortoncmp = MortonCompare<ITYPE>(rowlowbits, collowbits, lowrowmask, lowcolmask);
}

// Constructing empty BiCsb objects (size = 0) are not allowed.
template <class T, class ITYPE>
BiCsb<T, ITYPE>::BiCsb (ITYPE size, ITYPE rows, ITYPE cols, int workers): nz(size),m(rows),n(cols)
{
	assert(nz != 0 && n != 0 && m != 0);
	Init(workers);

	num = new T[nz];
	bot = new ITYPE[nz];
	top = new ITYPE* [nbr];	

	for(ITYPE i=0; i<nbr; ++i)
		top[i] = new ITYPE[nbc+1]; 
}

// copy constructor
template <class T, class ITYPE>
BiCsb<T, ITYPE>::BiCsb (const BiCsb<T,ITYPE> & rhs)
: nz(rhs.nz), m(rhs.m), n(rhs.n), blcrange(rhs.blcrange), nbr(rhs.nbr), nbc(rhs.nbc), 
rowhighbits(rhs.rowhighbits), rowlowbits(rhs.rowlowbits), highrowmask(rhs.highrowmask), lowrowmask(rhs.lowrowmask), 
colhighbits(rhs.colhighbits), collowbits(rhs.collowbits), highcolmask(rhs.highcolmask), lowcolmask(rhs.xlowcolmask),
mortoncmp(rhs.mortoncmp), ispar(rhs.ispar)
{
	if(nz > 0)
	{
		num = new T[nz];
		bot = new ITYPE[nz];

		for(ITYPE i=0; i< nz; ++i)	
			num[i]= rhs.num[i];
		for(ITYPE i=0; i< nz; ++i)	
			bot[i]= rhs.bot[i];
	}
	if ( nbr > 0)
	{
		top = new ITYPE* [nbr];

		for(ITYPE i=0; i<nbr; ++i)
			top[i] = new ITYPE[nbc+1]; 

		for(ITYPE i=0; i<nbr; ++i)
			for(ITYPE j=0; j <= nbc; ++j) 
				top[i][j] = rhs.top[i][j];
	}
}

template <class T, class ITYPE>
BiCsb<T, ITYPE> & BiCsb<T, ITYPE>::operator= (const BiCsb<T, ITYPE> & rhs)
{
	if(this != &rhs)		
	{
		if(nz > 0)	// if the existing object is not empty
		{
			// make it empty
			delete [] bot;
			delete [] num;
		}
		if(nbr > 0)
		{
			for(ITYPE i=0; i<nbr; ++i)
				delete [] top[i];
			delete [] top;
		}

		ispar 	= rhs.ispar;
		nz	= rhs.nz;
		n	= rhs.n;
		m   	= rhs.m;
		nbr 	= rhs.nbr;
		nbc 	= rhs.nbc;
		blcrange = rhs.blcrange;

		rowhighbits = rhs.rowhighbits;
		rowlowbits = rhs.rowlowbits;
		highrowmask = rhs.highrowmask;
		lowrowmask = rhs.lowrowmask;

		colhighbits = rhs.colhighbits;
		collowbits = rhs.collowbits;
		highcolmask = rhs.highcolmask;
		lowcolmask= rhs.lowcolmask;
		mortoncmp = rhs.mortoncmp;

		if(nz > 0)	// if the copied object is not empty
		{
			num = new T[nz];
			bot = new ITYPE[nz];

			for(ITYPE i=0; i< nz; ++i)	
				num[i]= rhs.num[i];
			for(ITYPE i=0; i< nz; ++i)	
				bot[i]= rhs.bot[i];
		}
		if(nbr > 0)
		{
			top = new ITYPE* [nbr];

			for(ITYPE i=0; i<nbr; ++i)
				top[i] = new ITYPE[nbc+1]; 

			for(ITYPE i=0; i<nbr; ++i)
				for(ITYPE j=0; j <= nbc; ++j) 
					top[i][j] = rhs.top[i][j];
		}
	}
	return *this;
}

template <class T, class ITYPE>
BiCsb<T, ITYPE>::~BiCsb()
{
	if( nz > 0)
	{
		delete [] bot;
		delete [] num;
	}
	if ( nbr > 0)
	{
		for(ITYPE i=0; i<nbr; ++i)
			delete [] top[i];
		delete [] top;
	}
}

template <class T, class ITYPE>
BiCsb<T, ITYPE>::BiCsb (Csc<T, ITYPE> & csc, int workers):nz(csc.nz),m(csc.m),n(csc.n)
{
	typedef std::pair<ITYPE, ITYPE> ipair;
	typedef std::pair<ITYPE, ipair> mypair;

	assert(nz != 0 && n != 0 && m != 0);
	Init(workers);

	num = new T[nz];
	bot = new ITYPE[nz];
	top = new ITYPE* [nbr];
	for(ITYPE i=0; i<nbr; ++i)
		top[i] = new ITYPE[nbc+1]; 

	mypair * pairarray = new mypair[nz];
	ITYPE k = 0;
	for(ITYPE j = 0; j < n; ++j)
	{
		for (ITYPE i = csc.jc [j] ; i < csc.jc[j+1] ; ++i)	// scan the jth column
		{
			// concatenate the higher/lower order half of both row (first) index and col (second) index bits 
			ITYPE hindex = (((highrowmask &  csc.ir[i] ) >> rowlowbits) << colhighbits)
										| ((highcolmask & j) >> collowbits);
			ITYPE lindex = ((lowrowmask &  csc.ir[i]) << collowbits) | (lowcolmask & j) ;

			// i => location of that nonzero in csc.ir and csc.num arrays
			pairarray[k++] = mypair(hindex, ipair(lindex,i));
		}
	}
	sort(pairarray, pairarray+nz);	// sort according to hindex

	SortBlocks(pairarray, csc.num);

	delete [] pairarray;
}

// Assumption: rowindices (ri) and colindices(ci) are "parallel arrays" sorted w.r.t. column index values
template <class T, class ITYPE>
BiCsb<T, ITYPE>::BiCsb (ITYPE size, ITYPE rows, ITYPE cols, ITYPE * ri, ITYPE * ci, T * val, int workers, ITYPE forcelogbeta)
:nz(size),m(rows),n(cols)
{
	typedef std::pair<ITYPE, ITYPE> ipair;
	typedef std::pair<ITYPE, ipair> mypair;

	assert(nz != 0 && n != 0 && m != 0);
	Init(workers, forcelogbeta);

	num = new T[nz];
	bot = new ITYPE[nz];
	top = new ITYPE* [nbr];
	for(ITYPE i=0; i<nbr; ++i)
		top[i] = new ITYPE[nbc+1]; 

	mypair * pairarray = new mypair[nz];
	for(ITYPE k = 0; k < nz; ++k)
	{
		// concatenate the higher/lower order half of both row (first) index and col (second) index bits 
		ITYPE hindex = (((highrowmask &  ri[k] ) >> rowlowbits) << colhighbits)	| ((highcolmask & ci[k]) >> collowbits);	
		ITYPE lindex = ((lowrowmask &  ri[k]) << collowbits) | (lowcolmask & ci[k]) ;

		// k is stored in order to retrieve the location of this nonzero in val array
		pairarray[k] = mypair(hindex, ipair(lindex, k));
	}
	sort(pairarray, pairarray+nz);	// sort according to hindex
	
	SortBlocks(pairarray, val);
	
	delete [] pairarray;
}

template <class T, class ITYPE>
void BiCsb<T, ITYPE>::SortBlocks(pair<ITYPE, pair<ITYPE,ITYPE> > * pairarray, T * val)
{
 	typedef pair<ITYPE, pair<ITYPE, ITYPE> > mypair;	
	ITYPE cnz = 0;
	ITYPE ldim = IntPower<ITYPE>(2,colhighbits);	// leading dimension (not always equal to nbc)
	for(ITYPE i = 0; i < nbr; ++i)
	{
		for(ITYPE j = 0; j < nbc; ++j)
		{
			top[i][j] = cnz;
			ITYPE prevcnz = cnz; 
			std::vector<mypair> blocknz;

			while(cnz < nz && pairarray[cnz].first == ((i*ldim)+j) )	// as long as we're in this block
			{
				ITYPE lowbits = pairarray[cnz].second.first;
				ITYPE rlowbits = ((lowbits >> collowbits) & lowrowmask);
				ITYPE clowbits = (lowbits & lowcolmask);
				ITYPE bikey = BitInterleaveLow(rlowbits, clowbits);
				
				blocknz.push_back(mypair(bikey, pairarray[cnz++].second));
			}
			// sort the block into bitinterleaved order
			sort(blocknz.begin(), blocknz.end());

			for(ITYPE k=prevcnz; k<cnz ; ++k)
			{
				bot[k] = blocknz[k-prevcnz].second.first;
				num[k] = val[blocknz[k-prevcnz].second.second];
			}
		}
		top[i][nbc] = cnz;
	}
	assert(cnz == nz);
}

template <class T, class ITYPE>
float BiCsb<T, ITYPE>::RowImbalance() const
{
	assert(top[nbr-1][nbc] == nz);
	// get the average without the last left-over blockrow
	float rowave = static_cast<float>(*(top[nbr-1])) / (nbr-1);
	ITYPE rowmax = 0;
	for(ITYPE i=1; i< nbr; ++i)
	{
		rowmax = std::max(rowmax, *(top[i]) - *(top[i-1]));
	}
	return static_cast<float>(rowmax) / rowave;
}

template <class T, class ITYPE>
float BiCsb<T, ITYPE>::ColImbalance() const
{
	vector<float> sum(nbc-1);
	cilk_for(ITYPE j=1; j< nbc; ++j) 	// ignore the last block column
	{
		ITYPE * blocknnz = new ITYPE[nbr];  // nnz per block responsible
        	for(ITYPE i=0; i<nbr; ++i)
        	{
                	blocknnz[i] = top[i][j] - top[i][j-1];
        	}
        	sum[j-1] = std::accumulate(blocknnz, blocknnz + (nbr-1), 0);         // ignore the last block row
        	delete [] blocknnz;
	}
	float colave = std::accumulate(sum.begin(), sum.end(), 0.0) / static_cast<float>(nbc-1);
        vector<float>::iterator colmax = std::max_element(sum.begin(), sum.end());
	return (*colmax) / colave;
}



/**
  * @param[ITYPE**] chunks {an array of pointers, ith entry is an address pointing to the top array }
  * 	That address belongs to the the first block in that chunk
  * 	chunks[i] is valid for i = {start,start+1,...,end} 
  *	chunks[0] = btop
  **/ 
template <class T, class ITYPE>
void BiCsb<T, ITYPE>::BMult(ITYPE** chunks, ITYPE start, ITYPE end, const T * x, T * y, ITYPE ysize) const
{
	assert(end-start > 0);	// there should be at least one chunk
	if (end-start == 1) 	// single chunk
	{
		if((chunks[end] - chunks[start]) == 1)	// chunk consists of a single (normally dense) block 
		{
			ITYPE chi = ( (chunks[start] - chunks[0])  << collowbits);

			// m-chi > lowcolmask for all blocks except the last skinny tall one.
			// if the last one is regular too, then it has m-chi = lowcolmask+1
			if(ysize == (lowrowmask+1) && (m-chi) > lowcolmask )	// parallelize if it is a regular/complete block 	
			{
				const T * __restrict subx = &x[chi];
				BlockPar( *(chunks[start]) , *(chunks[end]), subx, y, 0, blcrange, BREAKEVEN * ysize);
			}
			else 		// otherwise block parallelization will fail 
			{
				SubSpMV(chunks[0], chunks[start]-chunks[0], chunks[end]-chunks[0], x, y);
			}
		}
		else 	// a number of sparse blocks with a total of at most O(\beta) nonzeros
		{
			SubSpMV(chunks[0], chunks[start]-chunks[0], chunks[end]-chunks[0], x, y);
		}  
	}
	else
	{
		// divide chunks into half 
		ITYPE mid = (start+end)/2;

		cilk_spawn BMult(chunks, start, mid, x, y, ysize);
		if(SYNCHED)
		{ 
			BMult(chunks, mid, end, x, y, ysize);
		}
		else
		{
			T * temp = new T[ysize];
			std::fill_n(temp, ysize, 0.0);

			BMult(chunks, mid, end, x, temp, ysize);
			cilk_sync;
			
			for(ITYPE i=0; i<ysize; ++i)
				y[i] += temp[i];

			delete [] temp;
		}
	}
}

template <class T, class ITYPE>
void BiCsb<T, ITYPE>::BTransMult(ITYPE col, ITYPE rowstart, ITYPE rowend, const T * x, T * y, ITYPE ysize) const
{
	ITYPE * blocknnz = new ITYPE[rowend-rowstart];	// nnz per block responsible
	for(ITYPE i=rowstart; i<rowend; ++i)
	{
		blocknnz[i-rowstart] = top[i][col+1] - top[i][col];
	}
	ITYPE sum = std::accumulate(blocknnz, blocknnz + (rowend-rowstart), 0);		// initial value = 0
	delete [] blocknnz;

	if(sum < BREAKEVEN * ysize)		// not enough nonzeros to amortize new array formation
	{
		// NOTE: even if a single block remains, it's sparse enough that we don't have to parallelize
		if(sum != 0)	// avoid function call and loop set up
			SubSpMVTrans(col, rowstart, rowend, x, y);

	}
	else if(rowend-rowstart == 1)	// only one block remains and it is relatively dense
	{
		if(ysize == (lowcolmask+1))     // parallelize if it is a regular (complete) block
		{
			ITYPE chi = (rowstart << rowlowbits);
			const T * __restrict subx = &x[chi];

			BlockParT( top[rowstart][col] , top[rowstart][col+1],  subx, y, 0, blcrange, BREAKEVEN * ysize);
		}
		else
		{
			SubSpMVTrans(col, rowstart, rowend, x, y);
		}
	}
	else
	{
		ITYPE rowmid = (rowend+rowstart)/2;	

		cilk_spawn BTransMult(col, rowstart, rowmid, x, y, ysize);
		if(SYNCHED)
		{
			BTransMult(col, rowmid, rowend, x, y, ysize);
		}
		else
		{
			T * temp = new T[ysize];
			std::fill_n(temp, ysize, 0.0);
		
			BTransMult(col, rowmid, rowend, x, temp, ysize);
			cilk_sync;
			
			for(ITYPE i=0; i<ysize; ++i)
				y[i] += temp[i];

			delete [] temp;
		}
	}
}

// double* restrict a; --> No aliases for a[0], a[1], ...
// bstart/bend: block start/end index (to the top array)
template <class T, class ITYPE>
void BiCsb<T, ITYPE>::SubSpMV(ITYPE * __restrict btop, ITYPE bstart, ITYPE bend, const T * __restrict x, T * __restrict suby) const
{
#ifdef STATS
        subspmvcalls() += 1;
#endif
	ITYPE * __restrict r_bot = bot;
	T * __restrict r_num = num;

	for (ITYPE j = bstart ; j < bend ; ++j)		// for all blocks inside that block row
	{
		// get higher order bits for column indices
		ITYPE chi = (j << collowbits);
		const T * __restrict subx = &x[chi];

		for (ITYPE k = btop[j] ; k < btop[j+1] ; ++k)	// for all nonzeros within ith block (expected =~ nnz/n = c)
		{
			ITYPE rli = ((r_bot[k] >> collowbits) & lowrowmask);
			ITYPE cli = (r_bot[k] & lowcolmask);
#ifdef BWTEST
			MultAdd<UNROLL> (suby[rli], r_num[k], subx[cli]);
#else
			suby [rli] += r_num[k] * subx [cli] ;
#endif
		}
	}
}

template <class T, class ITYPE>
void BiCsb<T, ITYPE>::SubSpMVTrans(ITYPE col, ITYPE rowstart, ITYPE rowend, 
							const T * __restrict x, T * __restrict suby) const
{
	ITYPE * __restrict r_bot = bot;
	T * __restrict r_num = num;

	for(ITYPE i= rowstart; i < rowend; ++i)
	{
		// get the starting point for accessing x
		ITYPE chi = (i << rowlowbits);
		const T * __restrict subx = &x[chi];

		for (ITYPE k = top[i][col] ; k < top[i][col+1] ; ++k)
		{
			// Note the swap in cli/rli
			ITYPE cli = ((r_bot[k] >> collowbits) & lowrowmask);
			ITYPE rli = (r_bot[k] & lowcolmask);
			suby [rli] += r_num[k] * subx [cli] ;
		}
	}	
}

// Parallelize the block itself (A*x version)
// start/end: element start/end positions (indices to the bot array)
// bot[start...end] always fall in the same block
// PRECONDITION: rangeend-rangebeg is a power of two 
// TODO: we rely on the particular implementation of lower_bound for correctness, which is dangerous !
//		 what if lhs (instead of rhs) parameter to the comparison object is the splitter?
template <class T, class ITYPE>
void BiCsb<T, ITYPE>::BlockPar(ITYPE start, ITYPE end, const T * __restrict subx, T * __restrict suby, 
							   ITYPE rangebeg, ITYPE rangeend, ITYPE cutoff) const
{
	assert(IsPower2(rangeend-rangebeg));
#ifdef STATS
	blockparcalls() += 1;
#endif
	if(end - start < cutoff)
	{
		ITYPE * __restrict r_bot = bot;
		T * __restrict r_num = num;

		for (ITYPE k = start ; k < end ; ++k)	
		{
			ITYPE rli = ((r_bot[k] >> collowbits) & lowrowmask);
			ITYPE cli = (r_bot[k] & lowcolmask);
#ifdef BWTEST
                        MultAdd<UNROLL> (suby[rli], r_num[k], subx[cli]);
#else                   
                        suby [rli] += r_num[k] * subx [cli] ;
#endif 			
		}
	}
	else
	{
		// Lower_bound is a version of binary search: it attempts to find the element value in an ordered range [first, last) 
		// Specifically, it returns the first position where value could be inserted without violating the ordering
		ITYPE halfrange = (rangebeg+rangeend)/2;
		ITYPE qrt1range = (rangebeg+halfrange)/2;
		ITYPE qrt3range = (halfrange+rangeend)/2;

		ITYPE * mid = std::lower_bound(&bot[start], &bot[end], halfrange, mortoncmp);
		ITYPE * left = std::lower_bound(&bot[start], mid, qrt1range, mortoncmp);
		ITYPE * right = std::lower_bound(mid, &bot[end], qrt3range, mortoncmp);

		/* -------
		   | 0 2 |
		   | 1 3 |
		   ------- */
		// subtracting two pointers pointing to the same array gives you the # of elements separating them
		// we're *sure* that the differences are 1) non-negative, 2) small enough to be indexed by an ITYPE
		ITYPE size0 = static_cast<ITYPE> (left - &bot[start]);
		ITYPE size1 = static_cast<ITYPE> (mid - left);
		ITYPE size2 = static_cast<ITYPE> (right - mid);
		ITYPE size3 = static_cast<ITYPE> (&bot[end] - right);

		ITYPE ncutoff = std::max<ITYPE>(cutoff/2, MINNNZTOPAR);
	    
		// We can choose to perform [0,3] in parallel and then [1,2] in parallel
		// or perform [0,1] in parallel and then [2,3] in parallel
		// Decision is based on the balance, i.e. we pick the more balanced parallelism
		if( ( absdiff(size0,size3) + absdiff(size1,size2) ) < ( absdiff(size0,size1) + absdiff(size2,size3) ) )
		{	
			cilk_spawn BlockPar(start, start+size0, subx, suby, rangebeg, qrt1range, ncutoff);	// multiply subblock_0
			BlockPar(end-size3, end, subx, suby, qrt3range, rangeend, ncutoff);					// multiply subblock_3
			cilk_sync;

			cilk_spawn BlockPar(start+size0, start+size0+size1, subx, suby, qrt1range, halfrange, ncutoff);	// multiply subblock_1
			BlockPar(start+size0+size1, end-size3, subx, suby, halfrange, qrt3range, ncutoff);				// multiply subblock_2
			cilk_sync;
		}
		else
		{
			cilk_spawn BlockPar(start, start+size0, subx, suby, rangebeg, qrt1range, ncutoff);		// multiply subblock_0
			BlockPar(start+size0, start+size0+size1, subx, suby, qrt1range, halfrange, ncutoff);	// multiply subblock_1
			cilk_sync;

			cilk_spawn BlockPar(start+size0+size1, end-size3, subx, suby, halfrange, qrt3range, ncutoff);	// multiply subblock_2
			BlockPar(end-size3, end, subx, suby, qrt3range, rangeend, ncutoff);								// multiply subblock_3
			cilk_sync;
		}
	}
}

// Parallelize the block itself (A'*x version)
// start/end: element start/end positions (indices to the bot array)
// bot[start...end] always fall in the same block
template <class T, class ITYPE>
void BiCsb<T, ITYPE>::BlockParT(ITYPE start, ITYPE end, const T * __restrict subx, T * __restrict suby, 
								   ITYPE rangebeg, ITYPE rangeend, ITYPE cutoff) const
{
	if(end - start < cutoff)
	{
		ITYPE * __restrict r_bot = bot;
		T * __restrict r_num = num;

		for (ITYPE k = start ; k < end ; ++k)	
		{
			// Note the swap in cli/rli
			ITYPE cli = ((r_bot[k] >> collowbits) & lowrowmask);
			ITYPE rli = (r_bot[k] & lowcolmask);
			suby [rli] += r_num[k] * subx [cli] ;
		}
	}
	else
	{
		ITYPE halfrange = (rangebeg+rangeend)/2;
		ITYPE qrt1range = (rangebeg+halfrange)/2;
		ITYPE qrt3range = (halfrange+rangeend)/2;

		// Lower_bound is a version of binary search: it attempts to find the element value in an ordered range [first, last) 
		// Specifically, it returns the first position where value could be inserted without violating the ordering
		ITYPE * mid = std::lower_bound(&bot[start], &bot[end], halfrange, mortoncmp);
		ITYPE * left = std::lower_bound(&bot[start], mid, qrt1range, mortoncmp);
		ITYPE * right = std::lower_bound(mid, &bot[end], qrt3range, mortoncmp);

		/* -------
		   | 0 1 |
		   | 2 3 |
		   ------- */
		// subtracting two pointers pointing to the same array gives you the # of elements separating them
		// we're *sure* that the differences are 1) non-negative, 2) small enough to be indexed by an ITYPE
		ITYPE size0 = static_cast<ITYPE> (left - &bot[start]);
		ITYPE size1 = static_cast<ITYPE> (mid - left);
		ITYPE size2 = static_cast<ITYPE> (right - mid);
		ITYPE size3 = static_cast<ITYPE> (&bot[end] - right);

		ITYPE ncutoff = std::max<ITYPE>(cutoff/2, MINNNZTOPAR);
	    
		// We can choose to perform [0,3] in parallel and then [1,2] in parallel
		// or perform [0,2] in parallel and then [1,3] in parallel
		// Decision is based on the balance, i.e. we pick the more balanced parallelism
		if( ( absdiff(size0,size3) + absdiff(size1,size2) ) < ( absdiff(size0,size2) + absdiff(size1,size3) ) )
		{	
			cilk_spawn BlockParT(start, start+size0, subx, suby, rangebeg, qrt1range, ncutoff);	// multiply subblock_0
			BlockParT(end-size3, end, subx, suby, qrt3range, rangeend, ncutoff);				// multiply subblock_3
			cilk_sync;

			cilk_spawn BlockParT(start+size0, start+size0+size1, subx, suby, qrt1range, halfrange, ncutoff);// multiply subblock_1
			BlockParT(start+size0+size1, end-size3, subx, suby, halfrange, qrt3range, ncutoff);				// multiply subblock_2
			cilk_sync;
		}
		else
		{
			cilk_spawn BlockParT(start, start+size0, subx, suby, rangebeg, qrt1range, ncutoff);	// multiply subblock_0
			BlockParT(start+size0+size1, end-size3, subx, suby, halfrange, qrt3range, ncutoff);	// multiply subblock_2
			cilk_sync;

			cilk_spawn BlockParT(start+size0, start+size0+size1, subx, suby, qrt1range, halfrange, ncutoff);// multiply subblock_1
			BlockParT(end-size3, end, subx, suby, qrt3range, rangeend, ncutoff);				// multiply subblock_3
			cilk_sync;
		}
	}
}

// Print stats to an ofstream object
template <class T, class ITYPE>
ofstream & BiCsb<T, ITYPE>::PrintStats(ofstream & outfile) const 
{
	if(nz == 0)
	{
		outfile << "## Matrix Doesn't have any nonzeros" <<endl;
		return outfile;
	}
	const ITYPE ntop = nbr * nbc; 	

	outfile << "## Average block is of dimensions "<< lowrowmask+1 << "-by-" << lowcolmask+1 << endl;
	outfile << "## Number of real blocks is "<< ntop << endl;
	outfile << "## Row imbalance is " << RowImbalance() << endl;
	outfile << "## Col imbalance is " << ColImbalance() << endl;
	
	std::vector<int> blocksizes(ntop);
	for(ITYPE i=0; i<nbr; ++i)
	{
		for(ITYPE j=0; j < nbc; ++j) 
		{
			blocksizes[i*nbc+j] = static_cast<int> (top[i][j+1]-top[i][j]);
		}
	}	
	sort(blocksizes.begin(), blocksizes.end());
	outfile<< "## Total nonzeros: "<< accumulate(blocksizes.begin(), blocksizes.end(), 0) << endl;

	outfile << "## Nonzero distribution (sorted) of blocks follows: \n" ;
	for(ITYPE i=0; i< ntop; ++i)
	{	
		outfile << blocksizes[i] << "\n";
	}
	outfile << endl;
	return outfile;
}


