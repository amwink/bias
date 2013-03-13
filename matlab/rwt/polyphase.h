/* 

   File Name: polyphase.h
   Author: Alle Meije Wink <a.wink@vumc.nl>
   
   polyphase decomposition in the frequency domain
   helper routines for other functions such as fsidwt, fisidwt

   - macros for complex numbers and operations
   - macros for complex and imaginary exponentials
   - routines for allocating matrices of (complex) doubles
   - routines for providing complex Fourier basis functions
   - routines for splitting signals into different phases
   - routines for merging phase signals back into the original

   In itself, polyphase decomposition is not more efficient in the frequency domain
   than in the time domain for filters much smaller than the input signal.
   However, in case the input signal is in the frequency domain, it is beneficial 
   to do the decomposition without going to the time domain.
   
   works only on column vectors (of the input matrix)
   does not work on [>2]D signals

*/

/***********************************************************************************/
/* include necesary header files                                                   */
/***********************************************************************************/

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/***********************************************************************************/
/* definitions                                                                     */
/***********************************************************************************/

/* the elusive number pi */
#define pi 3.141592653589793116

/* return 1 if x is an integer, 0 otherwise */
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)

/***********************************************************************************/
/* macros for handling complex numbers                                             */
/***********************************************************************************/

/* assign a dcomplex value z1 to a dcomplex variable z2 */
#define cassign(z2,z1) (z2.r = z1.r, z2.i = z1.i)

/* assign a dcomplex value (r1,r2) to a dcomplex variable z2 */
#define cassignrr(z,r1,r2) (z.r = r1, z.i = r2)

/* assign a dcomplex value z2 to two reals (r1,r2) */
#define rrassignc(r1,r2,z) (r1 = z.r, r2 = z.i)

/* add a dcomplex value z2 to a complex variable z1 */
#define addcc(z2,z1) ( z2.r += z1.r, z2.i += z1.i )

/* real and imaginary part of a complex exponential */
#define cexp(y,x) ( y.r = exp(x.r)*cos(x.i), y.i = exp(x.r)*sin(x.i)  )
#define r_cexp(x,y) (exp(x)*cos(y))
#define i_cexp(x,y) (exp(x)*sin(y))

/* real and imaginary part of an imaginary exponential */
#define iexp(y,x) ( y.r = cos(x), y.i = sin(x)  )
#define r_iexp(y) (cos(y))
#define i_iexp(y) (sin(y))

/* real and imaginary part of a complex product */
#define mulcc(c,a,b) ( c.r = a.r*b.r - a.i*b.i, c.i = a.i*b.r + a.r*b.i )
#define r_mulcc(a,b) (a.r*b.r-a.i*b.i)
#define i_mulcc(a,b) (a.i*b.r+a.r*b.i)

/* complex product output where 1 input is given as a <dComplex>, and 
                                1 input is given as two doubles: (re,im) */
#define mulcrr(d,a,b,c) (d.r = a.r*b - a.i*c, d.i = a.i*b + a.r*c)
#define r_mulcrr(a,b,c) (a.r*b-a.i*c)
#define i_mulcrr(a,b,c) (a.i*b+a.r*c)

/* add complex product to complex output (r1,r2) */
#define addmulrcc(r1,r2,z1,z2) (r1 += z1.r*z2.r - z1.i*z2.i, r2 += z1.i*z2.r + z1.r*z2.i)

/* add complex product to output where 1 input is given as a <dComplex>, and 
                                       1 input is given as two doubles: (re,im) */
#define addmulcrr(d,a,b,c) (d.r += a.r*b - a.i*c, d.i += a.i*b + a.r*c)

/***********************************************************************************/
/* type definitions for numbers, vectors and matrices                              */
/***********************************************************************************/

typedef struct {double r,i;} dComplex;

typedef dComplex *dComplexVec;
typedef dComplex **dComplexMat;

typedef double *doubleVec;
typedef double **doubleMat;

/***********************************************************************************/
/* allocation                                                                      */
/***********************************************************************************/

/* allocate a matrix of complex doubles */
dComplexMat makedComplexMat(int width, int height);

/* deallocate a matrix of complex doubles */
void freedComplexMat(dComplexMat theMatrix);

/* turn a 1D array of doubles into a matrix */ 
doubleMat DoubleMake2D( doubleVec array1D,
		        int width, int height);

/* turn a 1D array of complex doubles into a matrix */ 
dComplexMat dComplexMake2D( dComplexVec array1D,
			    int width, int height);

/***********************************************************************************/
/* shifts                                                                          */
/***********************************************************************************/

/* create shift exponentials for within each phase */
dComplexMat kExponentials(int Q, int MoverQ);

/* create shift exponentials for different frequency blocks */
dComplexMat lExponentials(int Q);

/* compute the required offsets for polyphase transforms at each level */
dComplexMat PolyOffsetExponentials(int len, int levels);

/* compute the required offsets for monophase transforms at each level */
dComplexMat MonoOffsetExponentials(int len, int levels);    

/***********************************************************************************/
/* transforms                                                                      */
/***********************************************************************************/

/* The column vectors of all matrices M, M1, and M2
   respectively, mentioned below, must be Fourier 
   transforms (frequency representations) of time signals. 
   The filters must also be in the freq. domain. */
   
/* convert column vectors of a matrix M into the polyphase representation */
dComplexMat PolyPhase( doubleVec xreal, doubleVec ximag, 
		       dComplexMat shifts, 
		       int MoverQ, int Q);

/* convert column vectors of a matrix M
   into the mono phase representation */
int MonoPhase( dComplexMat input, 
	       doubleVec outreal, doubleVec outimag, 
	       dComplexMat shifts, 
	       int MoverQ, int Q);

/* multiply the column vectors of matrix M by 
   frequency representations of two QM filters 
   and store both filtered versions of M */
dComplexMat FilterSplit( dComplexMat workspacec,
			 doubleMat hrMat, doubleMat hiMat,
			 int dim, int Q);

/* multiply the column vectors of two matrices M1 and M2
   by frequency representations of two QM filters and 
   store the sum of the filtered versions in M1 */
int FilterMerge( dComplexMat workspacec,dComplexMat workspaced,
		 doubleMat hrMat, doubleMat hiMat,
		 int MoverQ, int Q);
