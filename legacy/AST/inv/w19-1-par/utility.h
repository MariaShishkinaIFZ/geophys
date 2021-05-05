#ifndef _UTILITY_H
#define _UTILITY_H

#define __int64 long long
#include <climits>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

#if defined STATS || defined BWTEST
        #define CILKPARALLEL
#endif

#ifdef BWTEST
	#define UNROLL 100
#else
	#define UNROLL 1
#endif

#ifdef CILKPARALLEL
	#include <cilk.h>
	#define SYNCHED __cilkrts_synched()
	#define DETECT __cilkscreen_enable_checking()
	#define ENDDETECT __cilkscreen_disable_checking()
#else
	#define cilk_for for
	#define cilk_main main
	#define cilk_spawn
	#define cilk_sync
	#define SYNCHED true
#endif

#ifdef STATS
	#include <reducer_opadd.h>
	cilk::hyperobject< cilk::reducer_opadd<__int64> > blockparcalls;
	cilk::hyperobject< cilk::reducer_opadd<__int64> > subspmvcalls;
#endif

using namespace std;

#define SLACKNESS 8
#define KBYTE 1024
#define L2SIZE (256*KBYTE)	// less than half of the L2 Cache (L2 should hold x & y at the same time) - back to 256
#define CLSIZE 64
#define REPEAT 100
#define BREAKEVEN 3		// A block (or subblock) with less than (BREAKEVEN * dimension) nonzeros won't be parallelized
#define MINNNZTOPAR 128		// A block (or subblock) with less than MINNNZTOPAR nonzeros won't be parallelized
#define EPSILON 0.0001

// "absolute" difference macro that has no possibility of unsigned wrap
#define absdiff(x,y)   ( (x) > (y) ? (x-y) : (y-x))

template <typename ITYPE>
ITYPE CumulativeSum (ITYPE * arr, ITYPE size)
{
    ITYPE prev;
    ITYPE tempnz = 0 ;
    for (ITYPE i = 0 ; i < size ; ++i)
    {
		prev = arr[i];
		arr[i] = tempnz;
		tempnz += prev ;	
    }
    return (tempnz) ;		    // return sum
}


template <typename T>
T machineEpsilon()
{	
	T machEps = 1.0;
 	do {
       		machEps /= static_cast<T>(2.0);
       		// If next epsilon yields 1, then break, because current
       		// epsilon is the machine epsilon.
    	}
    	while ((T)(static_cast<T>(1.0) + (machEps/static_cast<T>(2.0))) != 1.0);
 
    	return machEps;
}


template < typename T >
struct absdiff : binary_function<T, T, T>
{
        T operator () ( T const &arg1, T const &arg2 ) const
        {
                using std::abs;
                return abs( arg1 - arg2 );
        }
};


template <typename VEC, typename ITYPE>
void Verify(VEC & control, VEC & test, string name, ITYPE m)
{
        vector<double>error(m);
        transform(&control[0], (&control[0])+m, &test[0], error.begin(), absdiff<double>());
        vector<double>::iterator max = max_element(error.begin(), error.end());
        cout << "Max error is: " << *max << " on y[" << max-error.begin()<<"]=" << test[max-error.begin()] << endl;
        double machEps = machineEpsilon<double>();
        cout << "Absolute machine epsilon is: " << machEps <<" and y[" << max-error.begin() << "]*EPSILON becomes "
                        << machEps * test[max-error.begin()] << endl;

	double sqrtm = sqrt(static_cast<double>(m));
	cout << "sqrt(n) * relative error is: " << abs(machEps * test[max-error.begin()]) * sqrtm << endl;
        if ( (abs(machEps * test[max-error.begin()]) * sqrtm) < abs(*max))
        {
                cout << "*** ATTENTION ***: error is more than sqrt(n) times the relative machine epsilon" << endl;
        }

#ifdef DEBUG
        cout << "\n Errors: ";  // Print any error larger than EPSILON
        remove_copy_if(error.begin(), error.end(), ostream_iterator<double>(cout, " " ), bind2nd(less<double>(), EPSILON));
        cout << "\n On locations: ";
        for(ITYPE i=0; i<m; ++i)
        {
                if(error[i] > EPSILON)
                {
                        cout << i << " ";
                }
        }
#endif
}


template <int D>
void MultAdd(double & a, const double & b, const double & c)
{
	for(int i=0; i<D; i++)
	{
		a += b * c;
	}	
	
}

// bit interleave x and y, and return result
// only the lower order bits of x and y are assumed valid
template <typename ITYPE>
ITYPE BitInterleaveLow(ITYPE x, ITYPE y)
{
	ITYPE z = 0; // z gets the resulting Morton Number.
	int ite = sizeof(z) * CHAR_BIT / 2;

	for (int i = 0; i < ite; ++i) 
	{
		// bitwise shift operations have precedence over bitwise OR and AND
  		z |= (x & (1 << i)) << i | (y & (1 << i)) << (i + 1);
	}
	return z;
}

// bit interleave x and y, and return result z (which is twice in size)
template <typename ITYPE, typename OTYPE>
OTYPE BitInterleave(ITYPE x, ITYPE y)
{
	OTYPE z = 0; // z gets the resulting Morton Number.
	int ite = sizeof(x) * CHAR_BIT;

	for (int i = 0; i < ite; ++i) 
	{
		// bitwise shift operations have precedence over bitwise OR and AND
  		z |= (x & (1 << i)) << i | (y & (1 << i)) << (i + 1);
	}
	return z;
}

template <typename ITYPE>
ITYPE IntPower(ITYPE base, ITYPE exponent)
{
	ITYPE i = 1; 
	ITYPE power = 1;

	while ( i <= exponent ) 
	{
		power *= base;
		i++;
	}
	return power;
}

// T should be uint32, uint64, int32 or int64; force concept requirement
template <typename T>
bool IsPower2(T x)
{
	return ( (x>0) && ((x & (x-1)) == 0));
}

unsigned int nextpoweroftwo(unsigned int v)
{
	// compute the next highest power of 2 of 32(or 64)-bit n
	// essentially does 1 << (lg(n - 1)+1).

	unsigned int n = v-1;

	// any "0" that is immediately right to a "1" becomes "1" (post: any zero has at least two "1"s to its left) 
	n |= n >> 1;

	// turn two more adjacent "0" to "1" (post: any zero has at least four "1"s to its left)
	n |= n >> 2;
	n |= n >> 4;	// post: any zero has at least 8 "1"s to its left
	n |= n >> 8;	// post: any zero has at least 16 "1"s to its left
	n |= n >> 16;	// post: any zero has at least 32 "1"s to its left

	return ++n;
}

// 64-bit version
// note: least significant bit is the "zeroth" bit
// pre: v > 0
unsigned int highestbitset(unsigned __int64 v)
{
	// b in binary is {10,1100, 11110000, 1111111100000000 ...}  
	const unsigned __int64 b[] = {0x2ULL, 0xCULL, 0xF0ULL, 0xFF00ULL, 0xFFFF0000ULL, 0xFFFFFFFF00000000ULL};
	const unsigned int S[] = {1, 2, 4, 8, 16, 32};
	int i;

	unsigned int r = 0; // result of log2(v) will go here
	for (i = 5; i >= 0; i--) 
	{
		if (v & b[i])	// highestbitset is on the left half (i.e. v > S[i] for sure)
		{
			v >>= S[i];
			r |= S[i];
		} 
	}
	return r;
}

__int64 highestbitset(__int64 v)
{
	if(v < 0)
	{
		cerr << "Indices can not be negative, aborting..." << endl;
		return -1;
	}
	else
	{
		unsigned __int64 uv = static_cast< unsigned __int64 >(v);
		unsigned __int64 ur = highestbitset(uv);
		return static_cast< __int64 > (ur);
	}
}
		
// 32-bit version 
// note: least significant bit is the "zeroth" bit
// pre: v > 0
unsigned int highestbitset(unsigned int v)
{
	// b in binary is {10,1100, 11110000, 1111111100000000 ...}  
	const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
	const unsigned int S[] = {1, 2, 4, 8, 16};
	int i;

	unsigned int r = 0; 
	for (i = 4; i >= 0; i--) 
	{
		if (v & b[i])	// highestbitset is on the left half (i.e. v > S[i] for sure)
		{
			v >>= S[i];
			r |= S[i];
		} 
	}
	return r;
}

int highestbitset(int v)
{
	if(v < 0)
	{
		cerr << "Indices can not be negative, aborting..." << endl;
		return -1;
	}
	else
	{	
		unsigned int uv = static_cast< unsigned int> (v);
		unsigned int ur = highestbitset(uv);
		return static_cast< int > (ur); 
	}
}

#endif

