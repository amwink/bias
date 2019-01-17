/**
 * \file fastecm.h
 * \brief compute fast eigenvector centrality map (fECM) of a 4D fMRI data set
 
 Author: AM Wink (C)
 Date: 13/04/13
 
*/

#ifndef CBA_FASTECM
#define CBA_FASTECM

#include <iostream>
#include <iomanip>

#include "nifti/nifti1_io.h"

typedef double num;

using namespace std;

/* BLOb import -- memory buffer to typename                                    */
template<typename BufElem, typename Data>

void blobImport ( void* niftibuffer, unsigned bufsize, Data *data ) {
  
  typedef typename Data::value_type value_type;

  BufElem*
    buffer = ( BufElem* ) niftibuffer;
  Data
    &dataref = *data;
  
  for ( unsigned i = 0; i < bufsize; i++ ) 
    dataref[i] = ( value_type ) ( buffer[i] );	
	
};

/* BLOb export -- typename to memory buffer (in nifti record)                  */
template<typename BufElem, typename Data>

void blobExport ( nifti_image* nim, Data *data ) {
  
  const unsigned
    bufsize = data->size();
  BufElem
    *buffer = new BufElem[bufsize];
  Data
    &dataref = *data;

  for ( unsigned i = 0; i < bufsize; i++ )
    buffer[i] = ( BufElem ) dataref[i];
  
  nim->data = ( void* ) buffer;
  
};

/* import bricks from nifti file                                               */
template<typename DataVec>

void getNiftiBricks ( nifti_image *nim, void *nifti_blob, unsigned bufsize, DataVec *vec ) {
  
  switch ( nim->datatype ) {

  case ( NIFTI_TYPE_UINT8 ) :
    blobImport<unsigned char> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_INT16 ) :
    blobImport<signed short> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_INT32 ) :
    blobImport<signed int> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_FLOAT32 ) :
    blobImport<float> ( nifti_blob, bufsize, vec );
    break;
  /*case(NIFTI_TYPE_COMPLEX64):
      blobImport<complex_float> ( nifti_blob, bufsize, vec);
      break;*/
  case ( NIFTI_TYPE_FLOAT64 ) :
    blobImport<double> ( nifti_blob, bufsize, vec );
    break;
  /*case(NIFTI_TYPE_RGB24):
      blobImport<rgb_byte> ( nifti_blob, bufsize, vec);
      break;*/
  case ( NIFTI_TYPE_INT8 ) :
    blobImport<signed char> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_UINT16 ) :
    blobImport<unsigned short> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_UINT32 ) :
    blobImport<unsigned int> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_INT64 ) :
    blobImport<signed long long> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_UINT64 ) :
    blobImport<unsigned long long> ( nifti_blob, bufsize, vec );
    break;
  case ( NIFTI_TYPE_FLOAT128 ) :
    blobImport<long double> ( nifti_blob, bufsize, vec );
    break;
  /*case(NIFTI_TYPE_COMPLEX128):
    // blobImport<complex_double> ( nifti_blob, bufsize, vec);
    // break;*/
  /*case(NIFTI_TYPE_COMPLEX256):
    // blobImport<complex_longdouble> ( nifti_blob, bufsize, vec);
    // break;*/
  /*case(NIFTI_TYPE_RGBA32):
    // blobImport<> ( nifti_blob, bufsize, vec);
    // break;*/
  default:
    cout << nim->fname << " has unsupported data type " << nim->datatype << endl;
    break;
  } /* switch ( nim->datatype )                                          */

} /* getNiftiBricks()                                                          */

/* export bricks to nifti file                                                 */
template<typename DataVec>

void setNiftiBricks ( nifti_image *nim, DataVec *vec ) {

  switch ( nim->datatype ) {

  case ( NIFTI_TYPE_UINT8 ) :
    blobExport<unsigned char> ( nim, vec );
    break;
  case ( NIFTI_TYPE_INT16 ) :
    blobExport<signed short> ( nim, vec );
    break;
  case ( NIFTI_TYPE_INT32 ) :
    blobExport<signed int> ( nim, vec );
    break;
  case ( NIFTI_TYPE_FLOAT32 ) :
    blobExport<float> ( nim, vec );
    break;
  /*case(NIFTI_TYPE_COMPLEX64):
    blobExport<complex_float> ( nim, vec);
    break;*/
  case ( NIFTI_TYPE_FLOAT64 ) :
    blobExport<double> ( nim, vec );
    break;
    /*case(NIFTI_TYPE_RGB24):
      blobExport<rgb_byte> ( nim, vec);
      break;*/
  case ( NIFTI_TYPE_INT8 ) :
    blobExport<signed char> ( nim, vec );
    break;
  case ( NIFTI_TYPE_UINT16 ) :
    blobExport<unsigned short> ( nim, vec );
    break;
  case ( NIFTI_TYPE_UINT32 ) :
    blobExport<unsigned int> ( nim, vec );
    break;
  case ( NIFTI_TYPE_INT64 ) :
    blobExport<signed long long> ( nim, vec );
    break;
  case ( NIFTI_TYPE_UINT64 ) :
    blobExport<unsigned long long> ( nim, vec );
    break;
  case ( NIFTI_TYPE_FLOAT128 ) :
    blobExport<long double> ( nim, vec );
    break;
  /*case(NIFTI_TYPE_COMPLEX128):
    blobExport<complex_double> ( nim, vec);
    break;*/
  /*case(NIFTI_TYPE_COMPLEX256):
    blobExport<complex_longdouble> ( nim, vec);
    break;*/
  /*case(NIFTI_TYPE_RGBA32):
    blobExport<> ( nim, vec);
    break;*/
  default:
    cout << nim->fname << " has unsupported data type " << nim->datatype << endl;
    break;
  } /* switch ( nim->datatype )                                          */

} /* setNiftiBricks()                                                          */

#endif /* end of CBA_FASTECM                                                   */
