
/****************************************************************************                  
Magnitudes of the wavelet coefficients in the specified sector
(selected by input mask)

The output is a row vector (for further histogram analysis).

Author: Aleksandra Pizurica <Aleksandra.Pizurica@telin.UGent.be> 
        Dept.TELIN, Ghent University
Last modified: February 2001.
*****************************************************************************/ 
   
  
#include "mex.h"
#include <math.h>


double absol(double a)
{
 double b;
 if(a<0)
   b=-a;
else
   b=a;
return b;
}

void workFcn(double *Y, double **B, double **W, int M, int N, int Length)
{  
 
int i,j,k=-1;
double a;

for (i=0; i<M; i++)
{
  for(j=0; j<N; j++)
  {
     if(B[i][j]==1)
	  {
	   k++;
      a=W[i][j];
      Y[k]=absol(a);
      }
   }
}


} /* End of workFcn */



/*----------------------------------------------------------------------------------*/

/*Defining the gateway function, i.e. the mex function*/

#define SECTOR_BIN  prhs[0]
#define DETAIL prhs[1]
#define ARRAY_OUT plhs[0]


void mexFunction(
               int nlhs,
               mxArray  *plhs[],
               int     nrhs,
               const mxArray  *prhs[]
               )
{
  int row, col, i, j, M, N, Length;
  double *ArrayOut, **SectorBin, **Detail;

  if (nlhs!=1)
      mexErrMsgTxt("program requires one  output argument [ARRAY_OUT].");
  if (nrhs!=2)
      mexErrMsgTxt("program requires two input arguments [SECTOR_BIN DETAIL]");
  

  M=mxGetM(DETAIL);
  N=mxGetN(DETAIL);
    
 
  /* Dynamic allocation of memory for input matrices*/

  SectorBin = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    SectorBin[i]=malloc(N*sizeof(double));

  Detail = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    Detail[i]=malloc(N*sizeof(double));


  /* Convert input matrices to 2x2 C array (MATLAB stores a two-dimensional matrix 
     in memory as a one-dimensional array) */

  for (col=0; col<N; col++){
    for (row=0; row<M; row++){
		SectorBin[row][col] = (mxGetPr(SECTOR_BIN))[row+col*M];
      Detail[row][col] = (mxGetPr(DETAIL))[row+col*M];
      }
  }

  Length=0;
  for (i=0;i<M;i++)
	  for (j=0;j<N;j++)
		  if(SectorBin[i][j]==1)
			  Length++;

ARRAY_OUT = mxCreateDoubleMatrix(Length, 1, mxREAL);  
  ArrayOut = mxGetPr(ARRAY_OUT);



/*Call workFcn function*/ 
workFcn(ArrayOut,SectorBin,Detail,M,N,Length);

for(i=0;i<M;i++)
 free(Detail[i]);
free(Detail);

for(i=0;i<M;i++)
 free(SectorBin[i]);
free(SectorBin);
 
}/*mexFcn*/





