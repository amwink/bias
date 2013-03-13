/* 

   implementation of helper routines for 
   polyphase decomposition in the frequency domain

   see polyphase.h

   (C) Alle Meije Wink (a.wink@vumc.nl)

*/

#include "polyphase.h"
#include <string.h>

/***********************************************************************************/
/* allocation                                                                      */
/***********************************************************************************/

/* allocate a matrix of complex doubles */
dComplexMat makedComplexMat(int width, int height) {
  register int i;

  dComplexMat theMatrix = (dComplexMat)mxCalloc(width,sizeof(dComplexVec));
  theMatrix[0] = (dComplexVec)mxCalloc(width*height, sizeof(dComplex));

  for(i=1;i<width;i++)
    theMatrix[i] = theMatrix[i-1] + height;
  
  return theMatrix;
}

/* deallocate a matrix of dComplex doubles */
void freedComplexMat(dComplexMat theMatrix) {
  mxFree(theMatrix[0]);
  mxFree(theMatrix);
}

/* make an extra index for 1D arrays */ 
doubleMat DoubleMake2D(doubleVec array1D,
		       int width, int height) {

  register int i;
  doubleMat theMatrix=(doubleMat)mxCalloc(width,sizeof(doubleVec));

  theMatrix[0]=(doubleVec)array1D;

  for(i=1;i<width;i++)
    theMatrix[i] = theMatrix[i-1] + height;  

  return theMatrix;
}

/* make an extra index for 1D arrays */ 
dComplexMat dComplexMake2D(dComplexVec array1D,
			int width, int height) {

  register int i;
  dComplexMat theMatrix=(dComplexMat)mxCalloc(width,sizeof(dComplexVec));

  theMatrix[0]=(dComplexVec)array1D;

  for(i=1;i<width;i++)
    theMatrix[i] = theMatrix[i-1] + height;  

  return theMatrix;
}

/***********************************************************************************/
/* shifts                                                                          */
/***********************************************************************************/

/* create shift exponentials for within each phase */
dComplexMat kExponentials(int Q, int MoverQ) {

  int M=Q*MoverQ;
  register int j,k;

  double twopi=2*pi;
  double twopiM=twopi/M;
  double twopijM=0,twopijkM;

  dComplexMat theMatrix=makedComplexMat(Q,MoverQ);

  for (j=0; j<Q; j++) {
    twopijkM=0;
    for (k=0; k<MoverQ; k++,twopijkM+=twopijM)
      iexp(theMatrix[j][k],twopijkM);     
    twopijM+=twopiM;
  }

  return theMatrix;
  
}

/* create shift exponentials for different frequency blocks */
dComplexMat lExponentials(int Q) {

  register int j,l;

  double twopi=2*pi;
  double twopiQ=twopi/Q;
  double twopijQ=0,twopijlQ;
  dComplexMat theMatrix=makedComplexMat(Q,Q);

  for (j=0; j<Q; j++) {
    twopijlQ=0;
    for (l=0; l<Q; l++, twopijlQ+=twopijQ) 
      iexp(theMatrix[j][l],twopijlQ); 
    twopijQ+=twopiQ;
  }

  return theMatrix;

}

/* compute the required offsets for polyphase transforms at each level */
dComplexMat PolyOffsetExponentials(int len, int levels) {    

  register int L,j,k,l;
  register int ColIndex=0;
  register int SigIndex;
  int wid=(int)pow(2,levels)-1;
  int Q=1,MoverQ=len,M=len;

  dComplexMat kshifts=makedComplexMat(wid,len);
  dComplexMat kpexpo,lpexpo;

  for (L=0;L<levels;L++) {
    
    kpexpo=kExponentials(Q,MoverQ);
    lpexpo=lExponentials(Q);

    for (j=0; j<Q; j++, ColIndex++)
      for (k=0; k<MoverQ; k++) {
	SigIndex=k;	
	for (l=0; l<Q; l++,SigIndex+=MoverQ)
	  mulcc(kshifts[ColIndex][SigIndex],kpexpo[j][k],lpexpo[j][l]);	         
      }
            
    freedComplexMat(kpexpo);
    freedComplexMat(lpexpo);
    
    Q*=2;
    MoverQ/=2;

  }
 
  return kshifts;
   
}

/* compute the required offsets for monophase transforms at each level */
dComplexMat MonoOffsetExponentials(int len, int levels) {    

  int Q=(int)pow(2,levels-1),M=len;      /* Q takes the value of Q for the highest level */
  double twopi=2*pi;
  double twopiM=twopi/M;
  double twopijM=0,twopijkM=0;
  register int j,k;

  dComplexMat kshifts=makedComplexMat(Q,M);
 
  for (j=0;j<Q;j++,twopijM+=twopiM) 
    for (k=0, twopijkM=0; k<M; k++, twopijkM+=twopijM) 
      iexp(kshifts[j][k],-twopijkM);

  return kshifts;

}

/***********************************************************************************/
/* transforms                                                                      */
/***********************************************************************************/

/* convert a 1-phase signal of length M*Q into a Q-phase signal of length M */
dComplexMat PolyPhase( doubleVec xreal, doubleVec ximag, 
		       dComplexMat shifts,
		       int MoverQ, int Q) {
  
  register int j,k,l;
  int M=MoverQ*Q;            
  int level=(int)(log(Q)/log(2));
  int kShiftsCol=(int)(pow(2,level))-1;

  dComplexMat output=makedComplexMat(Q,MoverQ);

  doubleMat xrMat,xiMat;
  dComplexMat shMat;

  if (Q==1) {
    
    /* Q==1 -> only one phase, just copy the signal */
    for (k=0;k<M;k++) 
      cassignrr(output[0][k],xreal[k],ximag[k]);
   
  }  else {
    
    xrMat=DoubleMake2D(xreal,Q,MoverQ);
    xiMat=DoubleMake2D(ximag,Q,MoverQ);

    /* put the contribution of each downsampled phase in a different signal */
    for (j=0; j<Q; j++) {
      
      shMat=dComplexMake2D(shifts[kShiftsCol++],Q,MoverQ);        
      
      for (k=0; k<MoverQ; k++) 
	for (l=0; l<Q; l++)   	 	  
	  addmulcrr(output[j][k],shMat[l][k],xrMat[l][k]/Q,xiMat[l][k]/Q);       

      mxFree(shMat);        
    }
    
    mxFree(xrMat);
    mxFree(xiMat);

  } 
  
  return output;

}

/* convert a Q-phase signal of length M into a 1-phase signal of length M*Q */
int MonoPhase( dComplexMat input, 
	       doubleVec yreal, doubleVec yimag, 
	       dComplexMat shifts, 
	       int MoverQ, int Q)      
{
  register int j,k,l;                   
  int M=MoverQ*Q;                        
    
  doubleMat yrMat,yiMat;
  dComplexMat shMat;

  if (Q==1) {
    
    /* Q==1 -> only one phase, just copy the signal */
    for (k=0;k<M;k++) 
      rrassignc(yreal[k],yimag[k],input[0][k]);
    
  } else {

    /* initialise the output signal */
    memset(yreal,0,M*sizeof(double));
    memset(yimag,0,M*sizeof(double));
    
    yrMat=DoubleMake2D(yreal,Q,MoverQ);
    yiMat=DoubleMake2D(yimag,Q,MoverQ);

    /* upsample each phase signal and interleave the signals */
    for (j=0; j<Q; j++) {
      
      shMat=dComplexMake2D(shifts[j],Q,MoverQ);        
      
      for (k=0; k<MoverQ; k++) 
	for (l=0; l<Q; l++) 
	  addmulrcc(yrMat[l][k],yiMat[l][k],shMat[l][k],input[j][k]);    
      
      mxFree(shMat);
      
    }
    
    mxFree(yrMat);
    mxFree(yiMat);    
  
  }
  freedComplexMat(input);
  
}
  
/* multiply the Q-phase signal of length M with 
   the filter coefficients  (0, Q-1, 2Q-1, .., MQ-1) 
   and split it into approximation and detail 
*/
dComplexMat FilterSplit(dComplexMat workspacec, 
			doubleMat hrMat, doubleMat hiMat,
			int MoverQ, int Q) {
  
  register int j,k;
  int xlen=MoverQ*Q;
  dComplex tmp;

  /* new workspace containing the detail signal */
  dComplexMat 
    workspaced=makedComplexMat(Q,MoverQ);
  
  /* little trick to skip through the filter */
  doubleMat 
    firstHr=DoubleMake2D(hrMat[0],MoverQ,Q),
    firstHi=DoubleMake2D(hiMat[0],MoverQ,Q),
    secndHr=DoubleMake2D(hrMat[1],MoverQ,Q),
    secndHi=DoubleMake2D(hiMat[1],MoverQ,Q);

  for (j=0; j<Q; j++) 
    
    for (k=0;k<MoverQ;k++) {      
      cassign(tmp,workspacec[j][k]);

      mulcrr(workspaced[j][k],tmp,secndHr[k][0],secndHi[k][0]);
      mulcrr(workspacec[j][k],tmp,firstHr[k][0],firstHi[k][0]);
    }

  mxFree(firstHr);
  mxFree(firstHi);
  mxFree(secndHr);
  mxFree(secndHi);

  return workspaced;

}

/* multiply the Q-phase signals of length M with 
   the filter coefficients  (0, Q-1, 2Q-1, .., MQ-1) 
   and merge them into a new approximation 
   clean up the detail signal
*/
int FilterMerge(dComplexMat workspacec,dComplexMat workspaced,
		doubleMat hrMat, doubleMat hiMat,
		int MoverQ, int Q) {
  
  register int j,k;
  int xlen=MoverQ*Q;
  dComplex tmp;

  /* little trick to skip through the filter */
  doubleMat 
    firstHr=DoubleMake2D(hrMat[0],MoverQ,Q),
    firstHi=DoubleMake2D(hiMat[0],MoverQ,Q),
    secndHr=DoubleMake2D(hrMat[1],MoverQ,Q),
    secndHi=DoubleMake2D(hiMat[1],MoverQ,Q);

  for (j=0; j<Q; j++)

    for (k=0; k<MoverQ; k++) {
      cassign(tmp,workspacec[j][k]);
      mulcrr(workspacec[j][k],tmp,firstHr[k][0]/2,firstHi[k][0]/2);

      cassign(tmp,workspaced[j][k]);      
      addmulcrr(workspacec[j][k],tmp,secndHr[k][0]/2,secndHi[k][0]/2);
    }

  mxFree(firstHr);
  mxFree(firstHi);
  mxFree(secndHr);
  mxFree(secndHi);  

}

