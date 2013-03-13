/* File Name: fisidwt.c
   
   (C) Alle Meije Wink <a.wink@vumc.nl>
   
   inverse redundant wavelet transform in the frequency domain on the column vectors of a matrix
   
   -- this is the direct route in http://dx.doi.org/10.1016/j.sigpro.2009.11.022
   Alle Meije Winka, Jos B.T.M. Roerdink (2010), Signal Processing 90(6): 1779 - 1787
   "Polyphase decompositions and shift-invariant discrete wavelet transforms in the frequency domain"
   
   applies fisidwt only along the first dimension of a 2D matrix (does not work on a row vector!)
   (calling mean, sum or fft on a 2D matrix also works on the columns, but does work on a row vector)
   
   does not work on [>2]D signals
   
   SYNTAX: y = fisidwt(a,b,h,l)
   
   outputs         a [m*n]  = approximation signals
                   b [ml*n] = detail signals
                   h        = freq. representation of scaling function (column 1) 
                                                      and wavelet (column 2)
                   l        = levels of decomposition
   outputs         y [m*n]  = 2D matrix of column vectors

   the format of the output is in accordance with that of mirdwt of the Rice Wavelet Toolbox (rwt)
   for 1D signals, and with multi1Dirdwt (modification of rdwt to do 1D transform on a matrix) for
   2D signals.
*/

#include "polyphase.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) 
{
  double *outputr, *outputi,       /* real and imaginary parts of output vector */
         *hr, *hi,                 /* real and imaginary parts of filter */
         *ycr, *yci, *ydr, *ydi;   /* real and imaginary parts of output signals */
  int clen, cwid,                  /* dim1 and dim2 of the input isignal */
      dlen, dwid,                  /* dim1 and dim2 of the input isignal */
      hlen, nh,                    /* dim1 and dim2 of the filter */
      dims, level;                 /* input dimensionality (for check), decomposition level */
  int i;
  double mtest;                    /* test the level of decomposition */

  /* check for correct # of input variables */
  if (nrhs<4) {
      mexErrMsgTxt("rFWD(A,D,H,L): 4 parameters required: signal X, filter H, level L");
      return;
  } 

  /* check for correct # of dimensions */
  dims=mxGetNumberOfDimensions(prhs[0]);
  if (dims>2) {
    mexErrMsgTxt("rFWD(A,D,H,L): dimensionality of A is too high");
    return;
  } 

  dims=mxGetNumberOfDimensions(prhs[1]);
  if (dims>2) {
    mexErrMsgTxt("rFWD(A,D,H,L): dimensionality of D is too high");
    return;
  } 

  /* get input signals and dimensions: yc*/
  ycr  = mxGetPr(prhs[0]);
  yci  = mxGetPi(prhs[0]);
  clen = mxGetM(prhs[0]); 
  cwid = mxGetN(prhs[0]); 

  /* get input signals and dimensions: yd */
  ydr  = mxGetPr(prhs[1]);
  ydi  = mxGetPi(prhs[1]);
  dlen = mxGetM(prhs[1]); 
  dwid = mxGetN(prhs[1]); 

  /* check if approximation has an imaginary part */
  if (!mxGetPi(prhs[0])) {
    /* mexWarnMsgTxt("Approximation signal has no complex part. Creating..."); */
    yci=mxCalloc(clen*cwid,sizeof(double));   
  }
  /* check if detail has an imaginary part */
  if (!mxGetPi(prhs[1])) {
    /* mexWarnMsgTxt("Detail signal has no complex part. Creating..."); */
    ydi=mxCalloc(dlen*dwid,sizeof(double));   
  }

  /* get filter and dimensions */
  hr = mxGetPr(prhs[2]);
  hi = mxGetPi(prhs[2]);
  hlen = mxGetM(prhs[2]); 
  nh = mxGetN(prhs[2]); 
  
  if (nh!=2) {
    mexErrMsgTxt("frequency filter needs 2 columns: highpass and lowpass");
    return;
  }

  /* check if filter has an imaginary part */
  if (!mxGetPi(prhs[2])) {
    /* mexWarnMsgTxt("Input signal has no complex part. Creating..."); */
    hi=mxCalloc(hlen*nh,sizeof(double));   
  }

  /* get levels of decomposition */
  level = (int) *mxGetPr(prhs[3]);
  if (level < 0) {
    mexErrMsgTxt("Max. level of decomposition must be non-negative");
    return;
  }

  /* Check the ROW dimension of input */
  if (clen > 1) {
    mtest = (double) clen/pow(2.0, (double) level);
    if (!isint(mtest)) {
      mexErrMsgTxt("Column vectors must have length m*2^(L)");
      return;
    }
  } else {
    mexErrMsgTxt("Cannot decompose signals of length 1");
    return;
  } 

  if ( (clen != dlen) || (dwid != (level*cwid) ) ) {
    mexErrMsgTxt("Dimensions of input signals and level of decomposition do not match");
    return;
  }

  /* Create matrix for approximation part */
  plhs[0] = mxCreateDoubleMatrix(clen,cwid,mxCOMPLEX);
  outputr     = mxGetPr(plhs[0]);
  outputi     = mxGetPi(plhs[0]);

  multiMRFWD1D(outputr, outputi, /* real parts and imaginary parts of signals */
	       clen, cwid,	 /* signal length and number of input signals */
	       hr, hi,           /* real and imaginary parts of filters */
	       level,            /* level of decomposition */
	       ycr, yci, 
	       ydr, ydi);        /* real and imaginary parts of approximation and detail signals */

} 

/* reconstruction of the redundant DWT in the frequency domain */
int multiMRFWD1D( doubleVec yreal, doubleVec yimag, 
		  int siglen, int nsig,
		  doubleVec hreal, doubleVec himag,
		  int lev,
		  doubleVec creal, doubleVec cimag,
		  doubleVec dreal, doubleVec dimag) {

  long s,d,L;
  long i,j,k,kk;
  long H=siglen/2,
       T=nsig*siglen;

  dComplexVec
    workspacec=(dComplexVec)mxCalloc(siglen,sizeof(dComplex)),     /* working memory for c(vec) */
    workspaced=(dComplexVec)mxCalloc(siglen,sizeof(dComplex));     /* working memory for d(vec) */

  dComplexMat
    hcomp=makedComplexMat(2*lev,siglen);                           /* working memory for filter */  

  doubleMat
    yrMat=DoubleMake2D(yreal,nsig,siglen),                        /*          memory for output */
    yiMat=DoubleMake2D(yimag,nsig,siglen),
    crMat=DoubleMake2D(creal,nsig,siglen),                         /*         memory for approx */
    ciMat=DoubleMake2D(cimag,nsig,siglen),
    drMat=DoubleMake2D(dreal,lev*nsig,siglen),                     /*         memory for detail */
    diMat=DoubleMake2D(dimag,lev*nsig,siglen);

  /* copy filter into complex array                                                             */
  /* make filters for every level: efficient if >1 signal                                       */

  for (k=0; k<(2*lev); k+=2) {

    kk=k+1;
    
    if (!k)

      for (i=0, j=siglen; i<siglen; i++, j++) {	
	cassignrr(hcomp[0][i],hreal[i]/2,himag[i]/2);
	cassignrr(hcomp[1][i],hreal[j]/2,himag[j]/2);
      } /* for(i) */
    
    else {

      /* copy filter from previous level                                                        */
      for (i=0, j=siglen; i<siglen; i++, j++) {	
	cassign(hcomp[k] [i],hcomp[k-2] [i]);
	cassign(hcomp[kk][i],hcomp[kk-2][i]);
      } /* for(i) */

      /* subsample filter and put in first half                                                 */
      for (i=0, j=0; i<H; i++, j+=2) {
	cassign(hcomp[k] [i],hcomp[k] [j]);
	cassign(hcomp[kk][i],hcomp[kk][j]);
      } /* for(i) */
      
      /* copy first half to second half                                                         */      
      for (i=0, j=H; i<H; i++, j++) {
	cassign(hcomp[k] [j],hcomp[k] [i]);
	cassign(hcomp[kk][j],hcomp[kk][i]);
	
      } /* for(i) */      
    
    } /* else */

  } /* for k */

  /* loop over signals and keep track of signal offset                                          */
  /* (the n-level transform is applied toe each separate signal)                                */

  for (s=0; s<nsig; s++) {
    
    /* loop over levels, keep track of detail offset                                            */
    /* (each level produces another detail signal)                                              */
    
    for (L=lev; L>0; L--) {
      
      /* if (L==lev) make c and d by multiplying input with (final) filter                      */
      /* otherwise make c,d by multiplying prev signals with (prev-level) filter                */ 

      if (L==lev) {
	
	/* start with original filters                                                          */
	k=2*(lev-1);
	kk=k+1;

	/* pointer to last detail signal                                                        */	
 	d=(lev-1)*T;
  	
	for (i=0, j=d; i<siglen; i++, j++) {
	  mulcrr(workspaced[i],hcomp[kk][i],drMat[s][j],diMat[s][j]);
	  mulcrr(workspacec[i],hcomp[k ][i],crMat[s][i],ciMat[s][i]);
	} /* for(i) */
	
      } else {
	
	/*  prepare filter for next level                                                       */
	k-=2;
	kk-=2;
	
	/* initialise next level: pointer to next detail signal                                 */
	d-=T;
	
	for (i=0, j=d; i<siglen; i++, j++) {
	  mulcrr(workspaced[i],hcomp[kk][i],drMat[s][j],diMat[s][j]);
	  mulcc (workspacec[i],hcomp[k] [i],workspacec[i]);
	} /* for(i) */
	
      } /* if (L==lev) */
      
      /* add current detail to current approximation                                            */
      for (i=0; i<siglen; i++) 
	addcc(workspacec[i],workspaced[i]);
      
    } /* for (L) */

    /* if final level reached: copy approx into output                                        */
    for (i=0; i<siglen; i++)
      rrassignc(yrMat[s][i],yiMat[s][i],workspacec[i]);
 
  } /* for (s) */

  mxFree(hcomp);
  mxFree(yrMat);
  mxFree(yiMat);
  mxFree(crMat);
  mxFree(ciMat);
  mxFree(drMat);
  mxFree(diMat);

} /* multiMRFWD1D */

