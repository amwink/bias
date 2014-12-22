/*********************************************************************************************
Function init_mask.c
Returns a binary edge mask Y={y1...yn} for the input detail image X={x1...xn}
based on thresholding the input magnitudes as follows:
yi=1 if |xi|>=T 
yi=0 if |xi|<T

Last modified: February 2001.
Author: Aleksandra Pizurica [Aleksandra.Pizurica@Telin.UGent.be] Ghent University, Dept. TELIN

Copyright: All software, documentation, and related files in this distribution
           are Copyright (c) 2001 Aleksandra Pizurica. All rights reserved.
**********************************************************************************************/ 


#include "mex.h"


void workFcn(double **X, double *Y, int M, int N, double Threshold)
{  
  int i,j,m,n;

  for(i=0; i<M; i++)
    for(j=0; j<N; j++)
       Y[i+j*M]=X[i][j];


  for(i=0; i<M; i++){
    for(j=0; j<N; j++){ 
       if(abs(X[i][j])<Threshold)
          Y[i+j*M]=0;    
       else
           Y[i+j*M]=1;             
    }
  }
             
}  /* End of workFcn */




/*Defining the gateway function, i.e. the mex function*/

#define ARRAY_IN  prhs[0]
#define THRESHOLD  prhs[1]
#define ARRAY_OUT plhs[0]

 
void mexFunction(
               int nlhs,
               mxArray  *plhs[],
               int     nrhs,
               const mxArray  *prhs[]
               )
{
  int row, col, i, M, N;
  double Threshold;
  double **ArrayIn, *ArrayOut;

  if (nlhs!=1)
      mexErrMsgTxt("init_mask requires one output argument [ARRAY_OUT].");
  if (nrhs!=2)
      mexErrMsgTxt("mexample requires two input arguments [ARRAY_IN THRESHOLD]");

  M=mxGetM(ARRAY_IN);
  N=mxGetN(ARRAY_IN);

  Threshold=mxGetScalar(THRESHOLD);
    
  /* Allocate memory for return matrix */
  ARRAY_OUT = mxCreateDoubleMatrix(M, N, mxREAL);  /* create a M x N
                                                            full, numeric,
                                                            real-valued
                                                            array */
  ArrayOut = mxGetPr(ARRAY_OUT);


  /* Dynamic allocation of memory for ArrayIn */

  ArrayIn = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    ArrayIn[i]=malloc(N*sizeof(double));

  /* Convert ARRAY_IN to a 2x2 C array (MATLAB stores a two-dimensional matrix 
     in memory as a one-dimensional array) */

         for (col=0; col < mxGetN(ARRAY_IN); col++)
            for (row=0; row < mxGetM(ARRAY_IN); row++)
               ArrayIn[row][col] = (mxGetPr(ARRAY_IN))[row+col*mxGetM(ARRAY_IN)];



/*Call workFcn function*/ 
workFcn(ArrayIn,ArrayOut,M,N,Threshold);

for(i=0;i<M;i++)
 free(ArrayIn[i]);
free(ArrayIn);
 
}/*mexFcn*/







