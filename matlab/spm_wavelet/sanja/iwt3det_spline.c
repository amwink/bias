
/*******************************************************************                   
                 2-D INVERSE REDUNDANT WAVELET TRANSFORM

Calculated with quadratic spline. 
Uses 3 detail images and a redundant transform.

Author: Aleksandra Pizurica TELIN/RUG 1998.

********************************************************************/ 
   
  
#include "mex.h"
#include <math.h>

void wt1_mz(double *a, double *h, double *c, int p, int N, int D);


void workFcn(double *Y, double **LL, double **LH, double **HL, double **HH, int Scale, int M, int N)
{  
 
int i,j;
double *Row_In, *Row_Out, *Column_In, *Column_Out;

double k_Filt[7]={0,-0.03125, -0.21875, -0.6875, 0.6875, 0.21875, 0.03125};
double hi_Filt[7]={0,0,0.125,0.375,0.375,0.125,0};
double k_Filt_Length=7, l_Filt_Length=7, hi_Filt_Length=7;

Row_In=malloc(N*sizeof(double)); 
Column_In=malloc(M*sizeof(double));
Row_Out=malloc(N*sizeof(double)); 
Column_Out=malloc(M*sizeof(double));  


for (i=0; i<M; i++){

   for(j=0; j<N; j++)
     Row_In[j]=LL[i][j];
   wt1_mz(Row_In, hi_Filt, Row_Out, Scale, N, hi_Filt_Length);
   for(j=0; j<N; j++)
      LL[i][j]=Row_Out[j];

   for(j=0; j<N; j++)
     Row_In[j]=LH[i][j];
   wt1_mz(Row_In, hi_Filt, Row_Out, Scale, N, hi_Filt_Length);
   for(j=0; j<N; j++)
      LH[i][j]=Row_Out[j];

   for(j=0; j<N; j++)
     Row_In[j]=HL[i][j];
   wt1_mz(Row_In, k_Filt, Row_Out, Scale, N, k_Filt_Length);
   for(j=0; j<N; j++)
      HL[i][j]=Row_Out[j];

   for(j=0; j<N; j++)
     Row_In[j]=HH[i][j];
   wt1_mz(Row_In, k_Filt, Row_Out, Scale, N, k_Filt_Length);
   for(j=0; j<N; j++)
      HH[i][j]=Row_Out[j];
}

for (j=0; j<N; j++){
   for(i=0; i<M; i++)
     Column_In[i]=LL[i][j];
   wt1_mz(Column_In, hi_Filt, Column_Out, Scale, M, hi_Filt_Length);
   for(i=0; i<M; i++)
      LL[i][j]=Column_Out[i];

   for(i=0; i<M; i++)
     Column_In[i]=LH[i][j];
   wt1_mz(Column_In, k_Filt, Column_Out, Scale, M, k_Filt_Length);
   for(i=0; i<M; i++)
      LH[i][j]=Column_Out[i];

   for(i=0; i<M; i++)
     Column_In[i]=HL[i][j];
   wt1_mz(Column_In, hi_Filt, Column_Out, Scale, M, hi_Filt_Length);
   for(i=0; i<M; i++)
      HL[i][j]=Column_Out[i];

   for(i=0; i<M; i++)
     Column_In[i]=HH[i][j];
   wt1_mz(Column_In, k_Filt, Column_Out, Scale, M, k_Filt_Length);
   for(i=0; i<M; i++)
      HH[i][j]=Column_Out[i];

}
    
for(i=0; i<M; i++)
  for(j=0; j<N; j++)
    Y[i+j*M]=LL[i][j]+LH[i][j]+HL[i][j]+HH[i][j]; 

free(Row_In);
free(Row_Out);
free(Column_In);
free(Column_Out);

} /* End of workFcn */


/*---------------------------------------------------------------------------------*/

void wt1_mz(double *a, double *h, double *c, int p, int N, int D)
{  
int i,j,Dext,I1,I2,len,t;
double *ap,*he;

t=pow(2,p);

ap=malloc(3*N*sizeof(double));
Dext=(D-1)*t+1;
he=malloc(Dext*sizeof(double));

/* Solving the border problems: Input vector should be made periodical and mirrored*/

 for (i=0; i<N; i++)
   {
     ap[i+N]=a[i];
     ap[3*N-1-i]=a[i];
     ap[N-1-i]=a[i];
   }

/* Extending the impulse response h: Putting 2^p-1 zeros between the coefficients of h*/

 for(i=0; i<Dext; i++)
    he[i]=0;

 for(i=0; i<D; i++)
    he[i*t]=h[i];

/* Convolution: c=conv(ap,he) */

 I2=((D-1)/2)*t;  /*index of the element h(0) in the vector he*/

 for(j=0; j<N; j++){
   c[j]=0;
   I1=N+j;          /*index of the element a(0+j) in the vector ap*/
   if ((I1+I2+1)>=Dext)
      len=Dext;
   else
      len=I1+I2+1;
   for(i=0; i<len; i++){
      c[j]+=he[i]*ap[I1+I2-i];
   }
 }
   

free(ap);
free(he);
       
}  /* End of wt1_mz*/


/*----------------------------------------------------------------------------------*/

/*Defining the gateway function, i.e. the mex function*/

#define LOWPASS prhs[0]
#define HOR_DETAIL prhs[1]
#define VER_DETAIL prhs[2]
#define DIAG_DETAIL prhs[3]
#define SCALE prhs[4]
#define ARRAY_OUT plhs[0]


void mexFunction(
               int nlhs,
               mxArray  *plhs[],
               int     nrhs,
               const mxArray  *prhs[]
               )
{
  int row, col, i, M, N;
  int Scale;
  double *ArrayOut, **LL, **LH, **HL, **HH;

  if (nlhs!=1)
      mexErrMsgTxt("iwt3det_spline requires one  output argument [ARRAY_OUT].");
  if (nrhs!=5)
      mexErrMsgTxt("iwt3det_spline requires five input arguments [LOWPASS HOR_DETAIL VER_DETAIL DIAG_DETAIL SCALE]");
  

  M=mxGetM(LOWPASS);
  N=mxGetN(LOWPASS);

  Scale=mxGetScalar(SCALE); 
  
  /* Allocate memory for return matrices */
  ARRAY_OUT = mxCreateDoubleMatrix(M, N, mxREAL);  
  ArrayOut = mxGetPr(ARRAY_OUT);
  
 
  /* Dynamic allocation of memory for input matrices*/

  LL = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    LL[i]=malloc(N*sizeof(double));

  LH = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    LH[i]=malloc(N*sizeof(double));

  HL = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    HL[i]=malloc(N*sizeof(double));

  HH = malloc(M*sizeof(int));
  for(i=0;i<M;i++)
    HH[i]=malloc(N*sizeof(double));



  /* Convert ARRAY_IN to a 2x2 C array (MATLAB stores a two-dimensional matrix 
     in memory as a one-dimensional array) */

  for (col=0; col<N; col++){
    for (row=0; row<M; row++){
        LL[row][col] = (mxGetPr(LOWPASS))[row+col*M];
        LH[row][col] = (mxGetPr(HOR_DETAIL))[row+col*M];
        HL[row][col] = (mxGetPr(VER_DETAIL))[row+col*M];
        HH[row][col] = (mxGetPr(DIAG_DETAIL))[row+col*M];

    }
  }


/*Call workFcn function*/ 
workFcn(ArrayOut,LL,LH,HL,HH,Scale,M,N);

for(i=0;i<M;i++)
 free(LL[i]);
free(LL);

for(i=0;i<M;i++)
 free(LH[i]);
free(LH);

for(i=0;i<M;i++)
 free(HL[i]);
free(HL);

for(i=0;i<M;i++)
 free(HH[i]);
free(HH);

 
}/*mexFcn*/





