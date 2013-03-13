/*

File Name: multi1Drdwt.c

(C) Alle Meije Wink <a.wink@vumc.nl>

redundant wavelet transform on the column vectors of a matrix

applies rdwt only along the first dimension of a 2D matrix
(compare calling mean, sum or fft to a 2D matrix)

does not work on [>2]D signals

SYNTAX: [a b] = multi1Drdwt(x,h,l);
inputs          x [m*n]  = 2D matrix of column vectors
	        h        = orthonormal scaling funtion
		l        = levels of decomposition
outputs         a [m*n]  = approximation signals
                b [m*ln] = detail signals

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

#define mat(a, i, j) (*(a + (dim1*(j)+i))) 
#define max(A,B) (A > B ? A : B)
#define min(A,B) (A < B ? A : B)
#define even(x)  ((x & 1) ? 0 : 1)
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])

{
  double *x, *h, *approx, *detail;
  int xdim1, xdim2, 
      hdim1, hdim2, 
      dims, hlen, level;
  double mtest;

  /* check for correct # of input variables */
  if (nrhs!=3){
    mexErrMsgTxt("3 parameters required: signal X, filter H, level L");
    return;
  }

  /* check for correct # of dimensions */
  dims=mxGetNumberOfDimensions(prhs[0]);
  if (dims>2) {
    mexErrMsgTxt("multi1Drdwt(x,h,l): dimensionality of x is too high");
    return;
  }

  /* get input signals and dimensions */
  x = mxGetPr(prhs[0]);
  xdim1 = mxGetM(prhs[0]); 
  xdim2 = mxGetN(prhs[0]); 

  /* get filter and dimensions */
  h = mxGetPr(prhs[1]);
  hdim1 = mxGetM(prhs[1]); 
  hdim2 = mxGetN(prhs[1]); 

  if ((hdim1*hdim2)<=1) {
    mexErrMsgTxt("Filter of 0 or 1 coefficients is useless");
    return;
  }
    
  if (min(hdim1,hdim2)>1) {
    mexErrMsgTxt("1D filter expected");
    return;
  }

  hlen = max(hdim1,hdim2);

  /* get levels of decomposition */
  level = (int) *mxGetPr(prhs[2]);
  if (level < 0) {
    mexErrMsgTxt("Max. level of decomposition must be non-negative");
    return;
  }

  /* Check the ROW dimension of input */
  if (xdim1 > 1) {
    mtest = (double) xdim1/pow(2.0, (double) level);
    if (!isint(mtest)) {
      mexErrMsgTxt("The matrix row dimension must be of size m*2^(L)");
      return;
    } 
  } else {
    mexErrMsgTxt("Cannot decompose signals of length 1");
    return;
  }

  /* Create matrix for approximation part */
  plhs[0] = mxCreateDoubleMatrix(xdim1,xdim2,mxREAL);
  approx  = mxGetPr(plhs[0]);

  /* Create matrix for detail part */
  plhs[1] = mxCreateDoubleMatrix(xdim1,level*xdim2,mxREAL);
  detail  = mxGetPr(plhs[1]);

  multiMRDWT1D(x, xdim1, xdim2, h, hlen, level, approx, detail);
}

multiMRDWT1D(double *input, 
	     int dim1, int dim2, 
	     double *Phi, 
	     int lPhi, 
	     int lev,
	     double *yl, double *yh)
{
  double  *phi, *psi, *xdummyl, *ydummyll, *ydummyhh;
  long i;
  int levIter, newDim1, dOffset, col, row;
  int rownum, rbCol, sample_f;

  /* allocate working vectors */
  xdummyl = (double *)mxCalloc(max(dim1,dim2)+lPhi-1,sizeof(double));
  ydummyll = (double *)mxCalloc(max(dim1,dim2),sizeof(double));
  ydummyhh = (double *)mxCalloc(max(dim1,dim2),sizeof(double));

  /* allocate space for filters */
  phi = (double *)mxCalloc(lPhi,sizeof(double));
  psi = (double *)mxCalloc(lPhi,sizeof(double));

  /* analysis lowpass and highpass */
  for (i=0; i<lPhi; i++){
    phi[i] = Phi[lPhi-i-1];
    psi[i] =Phi[i];
  }
  for (i=0; i<lPhi; i+=2)
    psi[i] = -psi[i];
  
  /* keep track of the dimensions */
  newDim1 = dim1;

  /* copy original into approximation */
  for (i=0; i<dim1*dim2; i++)
    yl[i] = input[i];
  
  /* main loop */
  sample_f = 1;
  for (levIter=0; levIter < lev; levIter++){
    /* actual (level dependent) column offset */
    /* only one detail channel per resolution */
    dOffset = dim2*(levIter);
      
    /* go by columns in case of a 2D signal*/
    rbCol = dim1/newDim1;               /* # of row blocks per column */

    for (col=0; col<dim2; col++){          /* loop over column           */

      for (rownum=0; rownum<rbCol; rownum++){  /* loop within one column     */

	/* store in dummy variables */
	row = -sample_f + rownum;
	for (i=0; i<newDim1; i++){    
	  row = row + sample_f;
	  xdummyl[i] = mat(yl, row, col);  
	}

	/* perform filtering: frowst LL/LH, then HL/HH */
	fpconv(xdummyl, newDim1, phi, psi, lPhi, ydummyll, ydummyhh); 

	/* restore dummy variables in matrices */
	row = -sample_f + rownum;
	for (i=0; i<newDim1; i++){    
	  row = row + sample_f;
	  mat(yl, row, col) = ydummyll[i];  
	  mat(yh, row, dOffset+col) = ydummyhh[i];  
	}
      }
    }
    sample_f = sample_f*2;
    newDim1 = newDim1/2;
  }
}

fpconv(double *x_in, int lx, double *h0, double *h1, int lh,
       double *x_outl, double *x_outh)
{
  int i, j;
  double x0, x1;

  for (i=lx; i < lx+lh-1; i++)
    x_in[i] = x_in[i-lx];
  for (i=0; i<lx; i++){
    x0 = 0;
    x1 = 0;
    for (j=0; j<lh; j++){
      x0 = x0 + x_in[j+i]*h0[lh-1-j];
      x1 = x1 + x_in[j+i]*h1[lh-1-j];
    }
    x_outl[i] = x0;
    x_outh[i] = x1;
  }
}
