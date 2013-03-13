/* 
   File Name: rfwd.c
   (C) Alle Meije Wink <a.wink@vumc.nl>
   
   redundant wavelet transform in the frequency domain on the column vectors of a matrix

   -- this is the polyphase decomposition route in http://dx.doi.org/10.1016/j.sigpro.2009.11.022
   Alle Meije Winka, Jos B.T.M. Roerdink (2010), Signal Processing 90(6): 1779 - 1787
    "Polyphase decompositions and shift-invariant discrete wavelet transforms in the frequency domain"
   
   applies rfwd only along the first dimension of a 2D matrix (does not work on a row vector!)
   (calling mean/sum/fft on a matrix also works on the columns, but does work on a row vector)
   
   does not work on [>2]D signals
   
   SYNTAX: [a b] = multi1Drfwd(x,h,l)
   
   inputs          x [m*n]  = 2D matrix of column vectors
                   h        = freq. representation of scaling function (column 1) and wavelet (column 2)
                   l        = levels of decomposition
   outputs         a [m*n]  = approximation signals
                   b [m*ln] = detail signals
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
  } 

  /* check for correct # of dimensions */
  dims=mxGetNumberOfDimensions(prhs[0]);
  if (dims>2) {
    mexErrMsgTxt("rFWD(X,H,L): dimensionality of X is too high");
    return;
  } 

  /* get input signals and dimensions */
  inputr = mxGetPr(prhs[0]);
  inputi = mxGetPi(prhs[0]);
  len = mxGetM(prhs[0]); 
  nsig = mxGetN(prhs[0]); 

  /* check if input has an imaginary part */
  if (!mxGetPi(prhs[0])) {
    /* mexWarnMsgTxt("Input signal has no complex part. Creating..."); */
    inputi=mxCalloc(len*nsig,sizeof(double));   
  }
    
  /* get filter and dimensions */
  hr = mxGetPr(prhs[1]);
  hi = mxGetPi(prhs[1]);
  hlen = mxGetM(prhs[1]); 
  nh = mxGetN(prhs[1]); 
    
  if (nh!=2) {
    mexErrMsgTxt("frequency filter needs 2 columns: highpass and lowpass");
    return;
  }

  /* check if filter has an imaginary part */
  if (!mxGetPi(prhs[1])) {
    /* mexWarnMsgTxt("Input signal has no complex part. Creating..."); */
    hi=mxCalloc(len*nsig,sizeof(double));   
  }

  /* get levels of decomposition */
  level = (int) *mxGetPr(prhs[2]);
  if (level < 0) {
    mexErrMsgTxt("Max. level of decomposition must be non-negative");
    return;
  }

  /* Check the ROW dimension of input */
  if (len > 1) {
    mtest = (double) len/pow(2.0, (double) level);
    if (!isint(mtest)) {
      mexErrMsgTxt("Column vectors must have length m*2^(L)");
      return;
    }
  } else {
    mexErrMsgTxt("Cannot decompose signals of length 1");
    return;
  } 

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
	       ycr, yci, ydr, ydi);       /* real and imaginary parts of approximation and detail signals */

} 

/* computation of the redundant DWT in the frequency domain */
int multiMRFWD1D( doubleVec xreal, doubleVec ximag, 
		  int siglen, int nsig,
		  doubleVec hreal, doubleVec himag,
		  int lev,
		  doubleVec creal, doubleVec cimag,
		  doubleVec dreal, doubleVec dimag) {

  int s,d,L;
  int Q,MoverQ;

  /* working memory */
  dComplexMat workspacec,workspaced;

  /* shifts */
  dComplexMat 
    PolyShifts=PolyOffsetExponentials(siglen,lev),
    MonoShifts=MonoOffsetExponentials(siglen,lev);

  doubleMat
    xrMat=DoubleMake2D(xreal,nsig,siglen),
    xiMat=DoubleMake2D(ximag,nsig,siglen),
    crMat=DoubleMake2D(creal,nsig,siglen),
    ciMat=DoubleMake2D(cimag,nsig,siglen),
    drMat=DoubleMake2D(dreal,lev*nsig,siglen),
    diMat=DoubleMake2D(dimag,lev*nsig,siglen),
    hrMat=DoubleMake2D(hreal,2,siglen),
    hiMat=DoubleMake2D(himag,2,siglen);


  /* for (d=0; d<siglen; d++) { 
     for (s=0; s<nsig; s++)
     printf("s[%d][%d] = (%0.2f,%0.2f) ",s,d,xrMat[s][d],xiMat[s][d]);
     printf("\n");
     }
  */
  
  /* loop over signals and keep track of signal offset           */
  /* (the n-level transform is applied toe each separate signal) */
  for (s=0; s<nsig; s++) {

    /* initialise level 0 */
    Q=1;
    MoverQ=siglen;
    d=s;
    
    /* loop over levels, keep track of detail offset */
    /* (each level requires another detail signal)   */
    for (L=0; L<lev; L++) {
      
      /* go to the multiphase representation */ 
      if (L>0)
	workspacec=PolyPhase(crMat[s],ciMat[s],PolyShifts,MoverQ,Q);
      else 
	workspacec=PolyPhase(xrMat[s],xiMat[s],PolyShifts,MoverQ,Q);

      /* multiply the phase signals with the filters */
      workspaced=FilterSplit(workspacec,hrMat,hiMat,MoverQ,Q); 
      
      /* go back to the monophase representation */    
      MonoPhase( workspacec,crMat[s],ciMat[s],MonoShifts,MoverQ,Q);
      MonoPhase( workspaced,drMat[d],diMat[d],MonoShifts,MoverQ,Q);
           
      /* initialise next level */
      Q*=2;     
      MoverQ/=2;
      d+=nsig;
    }

  }

  mxFree(xrMat);
  mxFree(xiMat);
  mxFree(crMat);
  mxFree(ciMat);
  mxFree(drMat);
  mxFree(diMat);
  mxFree(hrMat);
  mxFree(hiMat);

}
