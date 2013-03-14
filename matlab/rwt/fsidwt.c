/* File Name: fsidwt.c

   (C) Alle Meije Wink <a.wink@vumc.nl>
   
   redundant wavelet transform in the frequency domain on the column vectors of a matrix

   -- this is the direct route in http://dx.doi.org/10.1016/j.sigpro.2009.11.022
   Alle Meije Winka, Jos B.T.M. Roerdink (2010), Signal Processing 90(6): 1779 - 1787
   "Polyphase decompositions and shift-invariant discrete wavelet transforms in the frequency domain"
  
   applies fsidwt only along the first dimension of a 2D matrix (does not work on a row vector!)
   (calling mean, sum or fft on a 2D matrix also works on the columns, but does work on a row vector)
   
   does not work on [>2]D signals
   
   SYNTAX: [a b] = fisidwt(x,h,l)
   
   inputs          x [m*n]  = 2D matrix of column vectors
                   h        = freq. repr. of scaling (column 1) and wavelet (column 2)
                   l        = levels of decomposition
   outputs         a [m*n]  = approximation signals
                   b [ml*n] = detail signals
*/

#include "polyphase.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{
  double *inputr, *inputi,       /* real and imaginary parts of input vector */
         *hr, *hi,               /* real and imaginary parts of filter */
         *ycr, *yci, *ydr, *ydi; /* real and imaginary parts of output signals */
  int len, nsig,                 /* dim1 and dim2 of the input isignal */
      hlen, nh,                  /* dim1 and dim2 of the filter */
      dims, level;               /* input dimensionality (for check), decomposition level */
  int i;
  double mtest;                  /* test the level of decomposition */

  /* check for correct # of input variables */
  if (nrhs<3) {
      mexErrMsgTxt("rFWD(X,H,L): 3 parameters required: signal X, filter H, level L");
      return;
  } /* if (nrhs) */

  /* check for correct # of dimensions */
  dims=mxGetNumberOfDimensions(prhs[0]);
  if (dims>2) {
    mexErrMsgTxt("rFWD(X,H,L): dimensionality of X is too high");
    return;
  } /* if (dims) */

  /* get input signals and dimensions */
  inputr = mxGetPr(prhs[0]);
  inputi = mxGetPi(prhs[0]);
  len = mxGetM(prhs[0]); 
  nsig = mxGetN(prhs[0]); 

  /* check if input has an imaginary part */
  if (!mxGetPi(prhs[0])) {
    /* mexWarnMsgTxt("Input signal has no complex part. Creating..."); */
    inputi=mxCalloc(len*nsig,sizeof(double));   
  } /* if (mxGetPi) */
    
  /* get filter and dimensions */
  hr = mxGetPr(prhs[1]);
  hi = mxGetPi(prhs[1]);
  hlen = mxGetM(prhs[1]); 
  nh = mxGetN(prhs[1]); 
    
  if (nh!=2) {
    mexErrMsgTxt("frequency filter needs 2 columns: highpass and lowpass");
    return;
  } /* if (nh) */

  /* check if filter has an imaginary part */
  if (!mxGetPi(prhs[1])) {
    /* mexWarnMsgTxt("Input signal has no complex part. Creating..."); */
    hi=mxCalloc(len*nsig,sizeof(double));   
  } /* if (mxGetPi) */

  /* get levels of decomposition */
  level = (int) *mxGetPr(prhs[2]);
  if (level < 0) {
    mexErrMsgTxt("Max. level of decomposition must be non-negative");
    return;
  } /* if (level) */

  /* Check the ROW dimension of input */
  if (len > 1) {
    mtest = (double) len/pow(2.0, (double) level);
    if (!isint(mtest)) {
      mexErrMsgTxt("Column vectors must have length m*2^(L)");
      return;
    } /* if (isint) */
  } else {
    mexErrMsgTxt("Cannot decompose signals of length 1");
    return;
  }  /* else */

  /* Create matrix for approximation part */
  plhs[0] = mxCreateDoubleMatrix(len,nsig,mxCOMPLEX);
  ycr     = mxGetPr(plhs[0]);
  yci     = mxGetPi(plhs[0]);

  /* Create matrix for detail part */
  plhs[1] = mxCreateDoubleMatrix(len,level*nsig,mxCOMPLEX);
  ydr     = mxGetPr(plhs[1]);
  ydi     = mxGetPi(plhs[1]);

  multiMRFWD1D(inputr, inputi,            /* real parts and imaginary parts of signals */
	       len, nsig,		  /* signal length and number of input signals */
	       hr, hi,                    /* real and imaginary parts of filters */
	       level,                     /* level of decomposition */
	       ycr, yci, ydr, ydi);       /* real and imag parts of approx and detail signals */

}  /* mexFunction() */

/* computation of the redundant DWT in the frequency domain */
int multiMRFWD1D( doubleVec xreal, doubleVec ximag, 
		  int siglen, int nsig,
		  doubleVec hreal, doubleVec himag,
		  int lev,
		  doubleVec creal, doubleVec cimag,
		  doubleVec dreal, doubleVec dimag) {

  long s,d,L;
  long i,j,k,q,Q,NQ;
  long T=nsig*siglen;

  dComplexVec
    workspacec=(dComplexVec)mxCalloc(siglen,sizeof(dComplex)),     /* working memory for c(vec) */
    workspaced=(dComplexVec)mxCalloc(siglen,sizeof(dComplex));     /* working memory for d(vec) */

  dComplexMat
    Hfilter2d,Gfilter2d,Detail2d,Approx2d, 
    hcomp=makedComplexMat(2,siglen);                           /* working memory for filter */  
  
  doubleMat
    xrMat=DoubleMake2D(xreal,nsig,siglen),
    xiMat=DoubleMake2D(ximag,nsig,siglen),                         /*          memory for input */
    crMat=DoubleMake2D(creal,nsig,siglen),
    ciMat=DoubleMake2D(cimag,nsig,siglen),                         /*         memory for approx */
    drMat=DoubleMake2D(dreal,lev*nsig,siglen),
    diMat=DoubleMake2D(dimag,lev*nsig,siglen),                     /*         memory for detail */
    hrMat=DoubleMake2D(hreal,2,siglen),
    hiMat=DoubleMake2D(himag,2,siglen);                            /*         memory for filter */

  /* copy filter into complex array                                                             */

  for (i=0; i<siglen; i++) {	
    cassignrr(hcomp[0][i],hrMat[0][i],hiMat[0][i]);
    cassignrr(hcomp[1][i],hrMat[1][i],hiMat[1][i]);
  } /* for(i) */

  /* loop over signals and keep track of signal offset                                          */
  /* (the n-level transform is applied toe each separate signal)                                */

  for (s=0; s<nsig; s++) {
        
    /* loop over levels, keep track of detail offset                                            */
    /* (each level produces another detail signal)                                              */

    for (L=0; L<lev; L++) {      	          
      
      /* if (!level) make c and d by multiplying input with (original) filter                   */
      /* otherwise make c,d by multiplying prev approx with (next-level) filter                 */ 
      
      if (!L) {
	
	/*  no filter subsampling for first level                                               */	
	Q=1;
	NQ=siglen;

	/* pointer to first detail signal                                                       */	
 	d=s;
  	
	for (i=0; i<siglen; i++) {
	  mulcrr(workspaced[i],hcomp[1][i],xrMat[s][i],xiMat[s][i]);
	  mulcrr(workspacec[i],hcomp[0][i],xrMat[s][i],xiMat[s][i]);
	} /* for(i) */
	
      } else {
	
	/*  prepare filter subsampling for next level */
	d+=nsig;
	Q*=2;
	NQ/=2;

	/* initialise next level: pointer to next detail signal                                 */
	
	/* facilitate skipping through filter and data */
	Hfilter2d=(dComplexMat)dComplexMake2D(hcomp[0],NQ,Q);
	Gfilter2d=(dComplexMat)dComplexMake2D(hcomp[1],NQ,Q);
	Detail2d=(dComplexMat)dComplexMake2D(workspaced,Q,NQ);
	Approx2d=(dComplexMat)dComplexMake2D(workspacec,Q,NQ);

	for (q=0; q<Q; q++)         
	  for (i=0; i<NQ; i++) {     
	    mulcc(Detail2d[q][i],Gfilter2d[i][0],Approx2d[q][i]);
	    mulcc(Approx2d[q][i],Hfilter2d[i][0],Approx2d[q][i]);

	  } /* for(i) */	

	/* undo redimensionalisation */
	mxFree(Hfilter2d);
	mxFree(Gfilter2d);
	mxFree(Detail2d);
	mxFree(Approx2d);

      } /* if (!L) */
      
      /* copy detail into output                                                                */
      for (i=0; i<siglen; i++) 
	rrassignc(drMat[d][i],diMat[d][i],workspaced[i]);
      
    } /* for (L) */        

    /* if final level reached: copy approx into output                                          */
    for (i=0; i<siglen; i++)
      rrassignc(crMat[s][i],ciMat[s][i],workspacec[i]);
    
  } /* for (s) */
  
  mxFree(workspacec);
  mxFree(workspaced);
  mxFree(hcomp);
  mxFree(xrMat);
  mxFree(xiMat);
  mxFree(crMat);
  mxFree(ciMat);
  mxFree(drMat);
  mxFree(diMat);
  mxFree(hrMat);
  mxFree(hiMat);
  
} /* multiMRFWD1D */
