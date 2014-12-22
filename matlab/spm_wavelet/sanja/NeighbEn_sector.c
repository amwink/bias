
/****************************************************************************                  
Normalized energy of the neighboring wavelet coefficients in the specified sector
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

void workFcn(double *Y, double **B, double **W, double **Wcoarser, int H, int M, int N, int Length)
{  
 
int i,j,k=-1,u,v;
double E, Nr_neighbors=(2*H+1)*(2*H+1)-1;

for (i=H; i<M-H; i++)
{
  for(j=H; j<N-H; j++)
  {
     if(B[i][j]==1)
	  {
	   k++;
      E=0;
      for (u=i-H; u<i+1+H; u++)
        for (v=j-H; v<j+1+H; v++)
          E+=W[u][v]*W[u][v];
      E=E-W[i][j]*W[i][j]+Wcoarser[i][j]*Wcoarser[i][j];
      E=E/(Nr_neighbors+1);
      E=pow(E,0.5);
      Y[k]=E;
      }
   }
}


} /* End of workFcn */



/*----------------------------------------------------------------------------------*/

/*Defining the gateway function, i.e. the mex function*/

#define SECTOR_BIN  prhs[0]
#define DETAIL prhs[1]
#define DETAIL_COARSER prhs[2]
#define WINDOW_HALFSIZE prhs[3]
#define ARRAY_OUT plhs[0]


void mexFunction(
               int nlhs,
               mxArray  *plhs[],
               int     nrhs,
               const mxArray  *prhs[]
               )
{
  int row, col, i, j, M, N, Length, H;
  double *ArrayOut, **SectorBin, **Detail, **Detail_coarser;

  if (nlhs!=1)
      mexErrMsgTxt("NeighbEn_sector requires one  output argument [ARRAY_OUT].");
  if (nrhs!=4)
      mexErrMsgTxt("NeighbEn_sector requires four input arguments [SECTOR_BIN DETAIL DETAIL_COARSER WINDOW_HALFSIZE]");
  

  M=mxGetM(DETAIL);
  N=mxGetN(DETAIL);

  H=mxGetScalar(WINDOW_HALFSIZE);

    
 
  /* Dynamic allocation of memory for input matrices*/

  SectorBin = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    SectorBin[i]=malloc(N*sizeof(double));

  Detail = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    Detail[i]=malloc(N*sizeof(double));

  Detail_coarser = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    Detail_coarser[i]=malloc(N*sizeof(double));



  /* Convert input matrices to 2x2 C array (MATLAB stores a two-dimensional matrix 
     in memory as a one-dimensional array) */

  for (col=0; col<N; col++){
    for (row=0; row<M; row++){
		SectorBin[row][col] = (mxGetPr(SECTOR_BIN))[row+col*M];
      Detail[row][col] = (mxGetPr(DETAIL))[row+col*M];
      Detail_coarser[row][col] = (mxGetPr(DETAIL_COARSER))[row+col*M];
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
workFcn(ArrayOut,SectorBin,Detail,Detail_coarser,H,M,N,Length);

for(i=0;i<M;i++)
 free(Detail[i]);
free(Detail);

for(i=0;i<M;i++)
 free(SectorBin[i]);
free(SectorBin);
 
}/*mexFcn*/





