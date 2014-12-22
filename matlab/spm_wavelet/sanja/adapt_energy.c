/***********************************************************************************************
Function adapt_energy.c

Estimate the probability of presence of a "signal of interest" for each wavelet coefficient
in the given subband (i.e., estimate the "edge probability map" for the input subband) 
using the method of [1].

Input: 
      DETAIL...........- input wavelet subband 
      DETAIL_COARSER...- the corresponding subband from the coarser scale
      LIKELIHOOD_RATIO - ratio of the conditional densities of coefficient magnitudes 
	                     p(m|H1)/p(m|H0)
      PRIOR_RATIO......- ratio of the conditional densities of the local spatial activity indicator z:
	                     p(z|H1)/p(z|H0)
      WINDOW_HALFSIZE..- parameter that describes Window size as (Window size - 1)/2
	                     e.g., 3x3 window WINDOW_HALFSIZE=1
						 e.g., 5x5 window WINDOW_HALFSIZE=2...
      PRIOR_FACTOR.....- prior ratio P(H1)/P(H0)

Output:
      PROB.............- edge probability map for the input subband DETAIL 

[1] A. Pizurica, W. Philips, I. Lemahieu and M. Acheroy, 
   "A Versatile Wavelet Domain Noise Filtration Technique for Medical Imaging," 
    IEEE, Trans. on Medical Imaging, vol. 22, no. 3, pp. 323--331, March 2003. 

Author: Aleksandra Pizurica <Aleksandra.Pizurica@telin.UGent.be>
        Ghent University, Dept. Telecommunications and Information Processing.
Last modified: February 2001.
**********************************************************************************************/

#include "mex.h"
#include <math.h>


void workFcn(double **X, double **Xcoarser, double *P, int M, int N, 
             double *Likelihood_ratio, double *Prior_ratio, int len_l, int len_p, double r, int W)
{ 
int i,j,k,u,v,L,m,H1,H2,V1,V2,e;
double Energy,Eta,Ksi,Count;
double *Measure;

L=M*N;
Measure=malloc(L*sizeof(double)); 

for(i=0; i<M; i++)
 for(j=0; j<N; j++)
   Measure[i+j*M] = abs(X[i][j]);


 for(i=0; i<M; i++){
	for(j=0; j<N; j++){
      
      k=i+j*M;
      m=Measure[k];
      if (m>len_l-1)
        m=len_l-1;
      Ksi=Likelihood_ratio[m];

   
      /* Determine the local window*/      
      if(i>=W)
         H1=W;
      else
         H1=i;

      if(i<M-W)
         H2=W;
      else
         H2=M-i-1;

      if(j>=W)
         V1=W;
      else
         V1=j;

      if(j<N-W)
         V2=W;
      else
         V2=N-j-1;
     /* The local window specified*/      

     Energy=0;
     Count=-1;

      for(u=i-H1; u<i+H2+1; u++){
	      for(v=j-V1; v<j+V2+1; v++){
                Energy+=X[u][v]*X[u][v];
                Count+=1;
         }
      }

      Energy=(Energy-X[i][j]*X[i][j]+Xcoarser[i][j]*Xcoarser[i][j])/(Count+1);
      Energy=pow(Energy,0.5);
      e=Energy;
      if (e>len_p-1)
        e=len_p-1;
      Eta=r*Prior_ratio[e];

	   P[i+j*M]=Ksi*Eta/(1+Ksi*Eta);
     	  
	}
}

free(Measure);
             
}  /* End of workFcn */




/*Defining the gateway function, i.e. the mex function*/

#define DETAIL  prhs[0]
#define DETAIL_COARSER  prhs[1]
#define LIKELIHOOD_RATIO  prhs[2]
#define PRIOR_RATIO  prhs[3]
#define WINDOW_HALFSIZE  prhs[4]
#define PRIOR_FACTOR  prhs[5]

#define PROB plhs[0]



 
void mexFunction(
               int nlhs,
               mxArray  *plhs[],
               int     nrhs,
               const mxArray  *prhs[]
               )
{
  int row, col, i, M, N, ITERA, len_l, len_p, W;
  double Prior_factor;
  double **Detail, **Detail_coarser, *Prob, *Likelihood_ratio, *Prior_ratio;

  if (nlhs!=1)
      mexErrMsgTxt("adapt_energy requires one output arguments [PROB].");
  if (nrhs!=6)
      mexErrMsgTxt("adapt_energy requires six input arguments [DETAIL DETAIL_COARSER LIKELIHOOD_RATIO PRIOR_RATIO WINDOW_HALFSIZE PRIOR_FACTOR]");

  M=mxGetM(DETAIL);
  N=mxGetN(DETAIL);
  len_l=mxGetN(LIKELIHOOD_RATIO);
  len_p=mxGetN(PRIOR_RATIO);

  W=mxGetScalar(WINDOW_HALFSIZE);
  Prior_factor=mxGetScalar(PRIOR_FACTOR);
    
  /* Allocate memory for return matrices */
  PROB = mxCreateDoubleMatrix(M, N, mxREAL);  
  Prob = mxGetPr(PROB);

  /* Dynamic allocation of memory for ArrayIn and Mask*/

  Detail = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    Detail[i]=malloc(N*sizeof(double));

  Detail_coarser = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    Detail_coarser[i]=malloc(N*sizeof(double));

  Likelihood_ratio = malloc(len_l*sizeof(double));
  Prior_ratio = malloc(len_p*sizeof(double));




  /* Convert ARRAY_IN and INPUT_MASK to 2x2 C arrays (MATLAB stores a two-dimensional matrix 
     in memory as a one-dimensional array) */

         for (col=0; col < N; col++)
            for (row=0; row < M; row++)
	      {
               Detail[row][col] = (mxGetPr(DETAIL))[row+col*M];
               Detail_coarser[row][col] = (mxGetPr(DETAIL_COARSER))[row+col*M];
	      }

         for (col=0; col < len_l; col++)
               Likelihood_ratio[col] = (mxGetPr(LIKELIHOOD_RATIO))[col];

         for (col=0; col < len_p; col++)
               Prior_ratio[col] = (mxGetPr(PRIOR_RATIO))[col];
    

/*Call workFcn function*/ 
workFcn(Detail,Detail_coarser,Prob,M,N,Likelihood_ratio,Prior_ratio,len_l,len_p,Prior_factor,W);

for(i=0;i<M;i++)
 free(Detail[i]);
free(Detail);

for(i=0;i<M;i++)
 free(Detail_coarser[i]);
free(Detail_coarser);
 
 
}/*mexFcn*/







