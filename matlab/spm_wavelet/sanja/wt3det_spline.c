
/*******************************************************************                   
                      2-D REDUNDANT WAVELET TRANSFORM

Calculated with quadratic spline.
It is a redundant (nondownsampled) transform, which produces a low pass
image and three detail images at each resolution scale.

Author: Aleksandra Pizurica TELIN/RUG 1998.

********************************************************************/ 
   
  
#include "mex.h"
#include <math.h>

void wt1_mz(double *a, double *h, double *c, int p, int N, int D);


void workFcn(double **X, double *LL, double *LH, double *HL, double *HH, int Scale, int M, int N)
{  
 
int i,j;
double *Row_In, *Row_Out, *Column_In, *Column_Out;
double Detail_Filt[5]={0,0.5,-0.5,0,0};
double Lowpass_Filt[5]={0.125,0.375,0.375,0.125,0};
double Detail_Filt_Length=5, Lowpass_Filt_Length=5;

Row_In=malloc(N*sizeof(double)); 
Column_In=malloc(M*sizeof(double));
Row_Out=malloc(N*sizeof(double)); 
Column_Out=malloc(M*sizeof(double));  

for (i=0; i<M; i++){
   for(j=0; j<N; j++)
     Row_In[j]=X[i][j];

   wt1_mz(Row_In, Detail_Filt, Row_Out, Scale, N, Detail_Filt_Length);
   for(j=0; j<N; j++)
      {
         HL[i+j*M]=Row_Out[j];
         HH[i+j*M]=Row_Out[j];
      }
      
   wt1_mz(Row_In, Lowpass_Filt, Row_Out, Scale, N, Lowpass_Filt_Length);
   for(j=0; j<N; j++)
      {
         LH[i+j*M]=Row_Out[j];
         LL[i+j*M]=Row_Out[j];
      } 
}
    
for (j=0; j<N; j++){
      
   for(i=0; i<M; i++)
      Column_In[i]=LL[i+j*M];
   wt1_mz(Column_In, Lowpass_Filt, Column_Out, Scale, M, Lowpass_Filt_Length);
   for(i=0; i<M; i++)
      LL[i+j*M]=Column_Out[i];
      
   for(i=0; i<M; i++)
      Column_In[i]=LH[i+j*M];
   wt1_mz(Column_In, Detail_Filt, Column_Out, Scale, M, Detail_Filt_Length);
   for(i=0; i<M; i++)
      LH[i+j*M]=Column_Out[i];
   
   for(i=0; i<M; i++)
      Column_In[i]=HL[i+j*M];
   wt1_mz(Column_In, Lowpass_Filt, Column_Out, Scale, M, Lowpass_Filt_Length);
   for(i=0; i<M; i++)
      HL[i+j*M]=Column_Out[i];
      
   for(i=0; i<M; i++)
      Column_In[i]=HH[i+j*M];
   wt1_mz(Column_In, Detail_Filt, Column_Out, Scale, M, Detail_Filt_Length);
   for(i=0; i<M; i++)
      HH[i+j*M]=Column_Out[i];
   
}

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

#define ARRAY_IN prhs[0]
#define SCALE prhs[1]
#define LOWPASS plhs[0]
#define HOR_DETAIL plhs[1]
#define VER_DETAIL plhs[2]
#define DIAG_DETAIL plhs[3]

 
void mexFunction(
               int nlhs,
               mxArray  *plhs[],
               int     nrhs,
               const mxArray  *prhs[]
               )
{
  int row, col, i, M, N;
  int Scale;
  double **ArrayIn, *LL, *LH, *HL, *HH;

  if (nlhs!=4)
      mexErrMsgTxt("wt3det_spline requires four output arguments [LOWPASS HOR_DETAIL VER_DETAIL DIAG_DETAIL].");
  if (nrhs!=2)
      mexErrMsgTxt("wt3det_spline requires two input arguments [ARRAY_IN SCALE]");
  

  M=mxGetM(ARRAY_IN);
  N=mxGetN(ARRAY_IN);

  Scale=mxGetScalar(SCALE); 
  
  /* Allocate memory for return matrices */
  LOWPASS = mxCreateDoubleMatrix(M, N, mxREAL);  
  LL = mxGetPr(LOWPASS);
  
  HOR_DETAIL = mxCreateDoubleMatrix(M, N, mxREAL);  
  LH = mxGetPr(HOR_DETAIL);

  VER_DETAIL = mxCreateDoubleMatrix(M, N, mxREAL);  
  HL = mxGetPr(VER_DETAIL);
  
  DIAG_DETAIL = mxCreateDoubleMatrix(M, N, mxREAL);  
  HH = mxGetPr(DIAG_DETAIL);


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
workFcn(ArrayIn,LL,LH,HL,HH,Scale,M,N);

for(i=0;i<M;i++)
 free(ArrayIn[i]);
free(ArrayIn);
 
}/*mexFcn*/