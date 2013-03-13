/*

File Name: multi1Dirdwt.c

(C) Alle Meije Wink <a.m.wink@vumc.nl>

applies irdwt only along the first dimension of a 2D 
matrix of approximation signals

does not work on [>2]D approximation signals

SYNTAX: r = multi1Drdwt(a,b,h,l);
inputs      a [m*n]  = approximation signals
            b [m*ln] = detail signals
	    h        = orthonormal scaling funtion
            l        = levels of decomposition
outputs     r [m*n]  = 2D matrix of reconstructed column vectors

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
  double *x, *h,  *approx, *detail;
  int hlen, level, dims,
      adim1, adim2, 
      ddim1, ddim2, 
      hdim1, hdim2;      
  double mtest;

  /* check for correct # of input variables */
  if (nrhs!=4){
    mexErrMsgTxt("4 parameters required: approximation C, detail D, filter H, level L");
    return;
  }

  /* check for correct # of dimensions of the approximation*/
  dims=mxGetNumberOfDimensions(prhs[0]);
  if (dims>2) {
    mexErrMsgTxt("multi1Drdwt(c,d,h,l): dimensionality of c is too high");
    return;
  }
  /* get approximation signal and dimensions */
  approx = mxGetPr(prhs[0]);
  adim1 = mxGetM(prhs[0]); 
  adim2 = mxGetN(prhs[0]); 

  /* get detail signal and dimensions */
  detail = mxGetPr(prhs[1]);
  ddim1 = mxGetM(prhs[1]); 
  ddim2 = mxGetN(prhs[1]); 

  /* get filter and dimensions */  
  h = mxGetPr(prhs[2]);
  hdim1 = mxGetM(prhs[2]); 
  hdim2 = mxGetN(prhs[2]); 

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
  level = (int) *mxGetPr(prhs[3]);
  if (level < 0) {
    mexErrMsgTxt("The number of levels, L, must be a non-negative integer");
    return;
  }  

  /* check for consistency of rows and columns of approx, detail */
  if (adim1>1) {
    if((adim1 != ddim1) | ((adim2*level) != ddim2)){
      mexErrMsgTxt("Dimensions of first two input matrices not consistent!");
      return;
    }
    mtest = (double) adim1/pow(2.0, (double) level);
    if (!isint(mtest)) {
      mexErrMsgTxt("The matrix row dimension must be of size m*2^(L)");
      return;
    }
  } else {
    mexErrMsgTxt("Cannot reconstruct signals of length 1");
    return;
  }

  /* Create matrix for the reconstructed signal */
  plhs[0] = mxCreateDoubleMatrix(adim1,adim2,mxREAL);
  x = mxGetPr(plhs[0]);

  multiMIRDWT1D(x, adim1, adim2, h, hlen, level, approx, detail);
}

multiMIRDWT1D(double *output, 
       int dim1, int dim2, 
       double *h, 
       int lPhi, 
       int lev,
       double *yl, double *yh)
{
  double  *phi, *psi, *xdummyl, *ydummyll, *ydummyhh;
  double *xdummyh, *xh;
  long i, j;
  int levIter, newDim1, dOffset, row, col, lPhimin1;
  int rownum, rbCol, sample_f;

  /* allocate working vectors */
  xdummyl = (double *)mxCalloc(max(dim1,dim2),sizeof(double));
  ydummyll = (double *)mxCalloc(max(dim1,dim2)+lPhi-1,sizeof(double));
  ydummyhh = (double *)mxCalloc(max(dim1,dim2)+lPhi-1,sizeof(double));

  /* allocate space for filters */
  phi = (double *)mxCalloc(lPhi,sizeof(double));
  psi = (double *)mxCalloc(lPhi,sizeof(double));

  /* analysis lowpass and highpass */
  for (i=0; i<lPhi; i++){
    phi[i] = h[i]/2;
    psi[i] = h[lPhi-i-1]/2;
  }
  for (i=1; i<=lPhi; i+=2)
    psi[i] = -psi[i];
  
  lPhimin1 = lPhi - 1;

  /* 2^L */
  sample_f = 1;
  for (i=1; i<lev; i++)
    sample_f = sample_f*2;
  newDim1 = dim1/sample_f;

  /* restore yl in x */
  for (i=0;i<dim1*dim2;i++)
    output[i] = yl[i];
  
  /* main loop */
  for (levIter=lev;levIter>0; levIter--){
    /* actual (level dependent) column offset */
    dOffset = dim2*(levIter-1);
    
    /* go by columns in case of a 2D signal*/
    rbCol = dim1/newDim1;                 /* # of row blocks per column */
    for (col=0; col<dim2; col++){            /* loop over column */
      for (rownum=0; rownum<rbCol; rownum++){    /* loop within one column */
	/* store in dummy variables */
	row = -sample_f + rownum;
	for (i=0; i<newDim1; i++){    
	  row = row + sample_f;
	  ydummyll[i+lPhimin1] = mat(output, row, col);  
	  ydummyhh[i+lPhimin1] = mat(yh, row, dOffset+col);  
	}
	/* perform filtering and adding: first LL/LH, then HL/HH */
	bpconv(xdummyl, newDim1, phi, psi, lPhi, ydummyll, ydummyhh); 
	/* store dummy variables in matrices */
	row = -sample_f + rownum;
	for (i=0; i<newDim1; i++){    
	  row = row + sample_f;
	  mat(output, row, col) = xdummyl[i];  
	}
      }
    }
    sample_f = sample_f/2;
    newDim1 = newDim1*2;
  }
}

bpconv(double *x_out, int lx, double *g0, double *g1, int lh,
       double *x_inl, double *x_inh)
{
  int i, j;
  double x0;
 
  for (i=lh-2; i > -1; i--){
    x_inl[i] = x_inl[lx+i];
    x_inh[i] = x_inh[lx+i];
  }
  for (i=0; i<lx; i++){
    x0 = 0;
    for (j=0; j<lh; j++)
      x0 = x0 + x_inl[j+i]*g0[lh-1-j] +
	x_inh[j+i]*g1[lh-1-j];
    x_out[i] = x0;
  }
}


