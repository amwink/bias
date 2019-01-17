/**
 * \file fastecm.cpp
 * \brief compute fast eigenvector centrality map (fECM) of a 4D fMRI data set
 */

/*
  
  fastecm.cpp
  
  Author: AM Wink (C)
  Date: 13/04/13
  
  usage  :                 fastecm [-v] <filename>
      -v :                 increased verbosity
  input  : <filename>      should be 4D
  output : <prefix>fastECM is a 3D map
  return :                 0 if everything OK

  to compile:
  1. install nifti I/O libraries ( here: headers in /usr/include/nifti, shared library /usr/lib/libniftiio.so )
     -> in debian: sudo apt-get install libnifti-dev
  2. install a good c++ compiler ( here: a version of gcc that supports c++11                                 )
     -> in debian: sudo apt-get install gcc-snapshot
  3. call the c++ compiler with the right options and libraries
     -> in debian: g++ -Ofast -march=native -std=c++11 -I/usr/include/nifti -o fastecm fastecm.cpp /usr/lib/libniftiio.so

  If you use this method/program in your research, please remember to cite this paper
   in the journal "Brain Connectivity":
  Alle Meije Wink, Jan C de Munck, Ysbrand D van der Werf, Odile A van den heuvel, Frederik Barkhof
   "Fast eigenvector centrality mapping of voxel-wise connectivity in functional MRI: implementation, validation and interpretation."
    URL: http://online.liebertpub.com/doi/abs/10.1089/brain.2012.0087

*/


#include<math.h>
#include<unistd.h>

#include<string>
#include<vector>
#include<limits>
#include<numeric>
#include<algorithm>

#ifndef CBA_FASTECM
#include "fastecm.h"
#endif

enum readdata {
  NO_READ_DATA,
  DO_READ_DATA
};

using namespace std;

int main ( int   argc,
           char* argv[]
	   ) {
  
  /* Parse command line                                                        */
  int
    copt    = 0,
    errflg  = 0,
    verbose = 0;
  extern int
    optind;
  
  for ( ;; ) {
    
    if ( ( copt = getopt ( argc, argv, "v" ) ) == EOF )
      break;
    
    switch ( copt ) {
      
    case 'v':
      verbose++;
      break;
      
    default:
      errflg++;
      break;
      
    } /* switch ( copt ) */
    
  } /* for ( ;; ) */
  
    /* if the command line could not be parsed, print usage instructions         */
  if ( errflg || argc < ( optind + 1 ) ) {
    
    printf ( "usage  :                 fastecm [-v] <filename>\n" );
    printf ( "    -v :                 increased verbosity\n" );
    printf ( "input  : <filename>      should be 4D\n" );
    printf ( "output : <prefix>fastECM is a 3D map\n" );
    printf ( "return :                 0 if everything OK\n" );
    printf ( "                         1 if cannot overwite output file\n" );
    exit ( 1 );
    
  }
  
  /* no internal warnings about the file format                                */
  nifti_set_debug_level ( 0 );
  
  /* parse filename                                                            */
  string
    inputfilename = argv[optind],
    outputfilename = nifti_makebasename ( argv[optind] );
  
  outputfilename = outputfilename + "fastECM" + nifti_find_file_extension ( argv[optind] );
  
  /* check if the input file is a proper nifti file                            */
  if ( !is_nifti_file ( inputfilename.c_str() ) ) {
    
    cout << "Error during nifti integrity check of input image " << inputfilename << endl;
    exit ( 1 );
    
  }
  
  bool input_compressed = nifti_is_gzfile ( inputfilename.c_str() );
  
  /* for reading in images, distinguish between compressed and uncompressed    */
  /* files:                                                                    */
  /* 1) uncompressed files can be read, masked and stored per volume, using    */
  /*    less memory than when all data is read at once                         */
  /* 2) compressed files need to be read in full, then decompressed, before    */
  /*    accessing the data, so it is inevitable to read all data               */
  
  nifti_image
    *inputimage = NULL;
  
  if ( input_compressed )
    inputimage  = nifti_image_read ( inputfilename.c_str(), DO_READ_DATA );
  else
    inputimage  = nifti_image_read ( inputfilename.c_str(), NO_READ_DATA );
  
  if ( inputimage == NULL ) {
    
    cout << "There is a problem with reading input image " << inputfilename << endl;
    nifti_image_free ( inputimage );
    return 1;
    
  }
  
  /* if file is not 4D (time series of 3D volumes) then exit                   */
  if ( inputimage->dim[0] < 4 ) {
    
    cout << "Image is not a 4D dataset: " << inputfilename << endl;
    nifti_image_free ( inputimage );
    return 1;
    
  }
  
  nifti_image
    *outputimage = nifti_copy_nim_info ( inputimage );
  
  /* compute the size of the whole 4D data set                                 */
  long unsigned
    totsize = 1;
  
  for ( size_t i = 1; i < ( inputimage->dim[0] ); i++ )
    totsize *= ( inputimage->dim[i] );
  
  /* and the size of 1 volume                                                  */
  unsigned
    volsize = inputimage->dim[1] * inputimage->dim[2] * inputimage->dim[3];
  vector < size_t >
    mask ( volsize );
  
  /* initialise the 4d data structure, with each timepoint data empty, size 0  */
  vector < vector < num > >
    data ( inputimage->dim[4], vector < num > ( 0 ) );
  unsigned
    masksize = 0;
  
  if ( !input_compressed ) {
    
    /* uncompressed file -> voxel data have not yet been read -->               */
    /*   read 1st vol as mask; assume nonzeros here also in others              */
    /*   read all the volumes in the time series separately                     */
    /*   after reading each volume, the in-mask voxels are used to a contiguous */
    /*   block starting at index 0 after which the data is resized to mask size */
    
    nifti_brick_list
      inputbrick;
    
    /* build mask while reading data                                            */
    /* mask is the minimum of each voxel's time series, binarised               */
    for ( int t = 0; t < size_t ( data.size() ); t++ )
      
      if ( nifti_image_load_bricks ( inputimage, 1, &t, &inputbrick ) ) {
	
	cout << " reading brick " << t << " of " << inputimage->dim[4] << " ... ";

	/* resize to whole volume (brick - including zeroes) and read the brick */
	data[t].resize ( volsize );
	getNiftiBricks ( inputimage, inputbrick.bricks[0], data[t].size(), &data[t] );
	
	cout << " done\r" << flush;
	
      } else {
	
	cout << "There is a problem with reading input image " << inputfilename << endl;
	nifti_image_free ( inputimage );
	return 1;
	
      }
    
  } else {
    
    /* compressed file -> voxel data have been read into a 1D array -->         */
    /*   each volume in the time series is a block of <volsize> points          */
    /*   use the 1st volume as mask; assume nonzeros here also in others        */
    /*   after reading each volume, the in-mask voxels are used to a contiguous */
    /*   block starting at index 0 after which the data is resized to mask size */
    
    vector < vector <num> >
      allvolumes ( 1, vector < num > ( data.size() * volsize ) );
    
    getNiftiBricks (inputimage, inputimage->data, allvolumes[0].size(), &allvolumes[0] );
    
    num
      *dataptr=&allvolumes[0][0];
    
    for ( int t = 0; t < ( data.size() ); t++ ) {
      
      data[t].resize ( volsize );
      for ( size_t i = 0; i < data[t].size(); i++ )
	data[t][i] = *dataptr++;
      
    }
    
  }

  /* increase each voxel's count for each volume in which it is nonzero        */
  for ( int t = 0; t < ( data.size() ); t++ ) 
    for ( size_t i = 0; i < data[t].size(); i++ )
      if ( data[t][i] > 0 )
	mask[i]++;
  
  /* now only store indices of voxels that have data in all volumes            */
  for ( size_t i = 0; i < mask.size(); i++ )
    if ( mask[i]==data.size() )
      mask[masksize++]=i; // amazingly this is safe (masksize never overtakes i)
  
  mask.resize(masksize);
  
  /* apply the mask to the data, storing in-mask voxels only                   */
  for ( size_t t = 0; t < size_t ( data.size() ); t++ ) {
    
    for ( size_t i = 0; i < masksize; i++ )
      data[t][i] = data[t][mask[i]];
    
    /* resize to only include the in-mask voxels                             */
    data[t].resize ( masksize );
    
  }
  
  /* now the data is in a 'matrix' of size <timepoints> x <brainvoxels         */
  /* (it's actually a vector of vectors -slightly less efficient, but all STL  */
  
  /* first, compute mean at every in-mask voxel                                */
  /* see http://stackoverflow.com/questions/14924912                           */
  cout << " computing mean image     ... ";
  
  num
    matlab_eps = pow ( 2., -52. );
  vector < num >
    colmean ( data[0].size(), 0. );
  
  /* first compute vector of column sums (add all rows [timepoints] together)  */
  for_each (  data.begin(),
	      data.end(),
	      [&] ( const vector < num >& row ) {
		transform ( row.begin(), row.end(),
			    colmean.begin(),
			    colmean.begin(),
			    [] ( num d1, num d2 ) { return d1 + d2; } ); }
	      );
  
  /* divide row vector of column sums by the number of rows                    */
  transform ( colmean.begin(), colmean.end(),
	      colmean.begin(),
	      bind2nd ( divides < num > (), data.size() )
	      );
  
  cout << " done \n";
    
  /* then compute variance at every in-mask voxel                              */
  cout << " computing variance image ... ";
  
  /* first compute  the sum of squares                                         */
  /* see http://stackoverflow.com/questions/14924912  */
  vector < num >
    colvar ( colmean.size(), 0. );
  
  for_each (  data.begin(), data.end(),
	      [&] ( const vector < num >& row ) { transform ( row.begin(), row.end(),
							      colvar.begin(),
							      colvar.begin(),
							      [] ( num d1, num d2 ) { return ( d1 * d1 + d2 ); } ); }
	      );
  
  /* divide row vector of sums of squares by the number of rows                */
  transform ( colvar.begin(), colvar.end(),
	      colvar.begin(),
	      bind2nd ( divides < num > (), data.size() )
	      );
  
  /* colmean now contains E(x) and colvar contains E(x^2) for every column     */
  /* population variance = E(x^2) - (E(x))^2                                   */
  transform ( colvar.begin(), colvar.end(),
	      colmean.begin(),
	      colvar.begin(),
	      [] (    num d1, num d2 ) { return ( d1 - d2 * d2 ); }
	      );
  
  /* finally, population var divides by n, sample var by                       */
  /* n-1 ->  to get sample var, multiply by n/(n-1)                            */
  transform ( colvar.begin(), colvar.end(),
	      colvar.begin(),
	      bind2nd ( multiplies < num > (), num(data.size())/num(data.size()-1) )
	      );
  
  cout << " done \n";
  
  /* then standard deviation at every in-mask voxel                            */
  cout << " computing std dev image  ... ";
  
  vector < num >
    colstd ( colvar.size() );
  
  /* see http://stackoverflow.com/questions/3280541                            */
  transform ( colvar.begin(), colvar.end(),
	      colstd.begin(),
	      static_cast< num ( * ) ( num ) > ( sqrt )
	      );
  
  /* add a tiny number to prevent divisions by 0                               */
  transform ( colstd.begin(), colstd.end(), colstd.begin(),
	      bind2nd ( plus < num > (), matlab_eps ) );
  
  cout << " done \n";
  
  /* before we start the eigenvector computation we prepare the data           */
  /* for fast ECM to work the data needs to have mean 0 and variance 1         */
  /* so we can just divide by the standard devation and subtract the mean      */
  cout << " normalising input data   ... ";
  
  for_each (  data.begin(), data.end(),
	      [&] ( vector < num >& row ) { transform ( row.begin(), row.end(),
							colmean.begin(),
							row.begin(),
							[] ( num d1, num d2 ) { return d1 - d2; } ); }
	      );
    
  for_each (  data.begin(), data.end(),
	      [&] ( vector < num >& row ) { transform ( row.begin(), row.end(),
							colstd.begin(),
							row.begin(),
							[] ( num d1, num d2 ) { return d1 / d2; } ); }
	      );
    
  // cout << "mean          = " << colmean[masksize-1] << " " << endl;
  // cout << "var           = " << colvar[masksize-1] << " " << endl;
  // cout << "mask(0)       = " << mask[0] <<endl;
  // cout << "data(mask(0)) = "; for ( auto d: data ) cout << d[0] << " ";
  
  /* also divide data by sqrt ( data.size() -1 ) to have unit diagonals in M   */
  for_each ( data.begin(), data.end(),
	     [&] (   vector < num >& row ) { transform ( row.begin(), row.end(),
							 row.begin(),
							 bind2nd ( divides < num > (), sqrt ( data.size() -1 ) ) ); }
	     );
  
  cout << " done \n";
  
  /* now the power iteration algorithm can start                               */
  cout << " running power iteration  ..." << flush;
  
  /* initialise 'map' vector with ones and scale to unit L2-norm               */
  vector < num >
    vcurr ( colstd.size(), 1. / sqrt ( ( num ) colstd.size() ) );
  
  unsigned
    iter = 0;
  
  num
    cnorm   = 0.,
    dnorm   = 1.,
    prevsum = 0.;
  
  while ( ( iter < 100 ) && ( dnorm > cnorm ) ) {
    
    /* temporary vector vprev [1xN] = vcurr, prevsum is its sum                */
    vector < num >
      vprev ( vcurr );
    
    prevsum = accumulate ( vprev.begin(), vprev.end(), 0. );
    
    /* temporary vector vcurr_1 [1xT] = product vprev * data'                  */
    vector < num >
      vcurr_1 ( data.size(), 0. );
    
    transform ( data.begin(), data.end(), vcurr_1.begin(),
		[&] ( vector < num > const & d ) {
		  return inner_product ( d.begin(), d.end(), vprev.begin(), 0.0 ); }
		);
    
    /* temporary vector vcurr_2 [1xN] = product vcurr_1 * data                 */
    /* see http://stackoverflow.com/questions/15212343                         */
    vector < num >
      vcurr_2 ( vcurr.size(), 0. );
    
    auto
      vcurr1_iter = vcurr_1.begin();
    
    for_each( data.begin(), data.end(),
	      [&](const vector < num >& row) {
		num vcurr1_val = *vcurr1_iter++;
		transform( row.begin(), row.end(), vcurr_2.begin(), vcurr_2.begin(),
			   [=] ( num a,  num b ) { return a * vcurr1_val + b; } ); }
	      );
    
    /* temporary vector vcurr_3 [1xN] = vcurr_2 + prevsum                      */
    /* equivalent to vprev * ( ( data' * data ) + 1 )                          */
    vector < num >
      vcurr_3 ( vcurr.size(), 0. );
    
    transform ( vcurr_2.begin(), vcurr_2.end(), vcurr_3.begin(),
		bind2nd ( plus < num > (), prevsum )
		);
    
    /* L2 norm of current estimate                                             */
    cnorm = sqrt ( inner_product ( vcurr_3.begin(), vcurr_3.end(), vcurr_3.begin(), 0. ) );
    
    /* normalise current estimate                                              */
    transform ( vcurr_3.begin(), vcurr_3.end(), vcurr.begin(),
		bind2nd ( divides < num > (), cnorm ) );
    
    /* store the difference between vcurr and vprev in vprev                   */
    transform ( vcurr.begin(), vcurr.end(), vprev.begin(), vprev.begin(),
		minus < num > () );
    
    /* L2 norm of difference                                                   */
    dnorm = sqrt ( inner_product ( vprev.begin(), vprev.end(), vprev.begin(), 0. ) );
    
    /* L2 norm of current estimate                                             */
    cnorm = matlab_eps *
      sqrt ( inner_product ( vcurr.begin(), vcurr.end(), vcurr.begin(), 0. ) );
    
    iter++;
    
    if (verbose)
      cout << "\n iteration " << setw(2) << iter <<
	", || v_i - v_(i-1) || / || v_i * epsilon || = " <<
	fixed << setprecision(16) << setfill('0') << dnorm << " / " <<
	fixed << setprecision(16) << setfill('0') << cnorm;
    
    
  } /* while ( ( iter < 100 ) && ( dnorm > numeric_limits< num >::min() ) )    */
  
  cout << "  done" << endl;
  
  /* sort vcurr to get intensities for percentiles 5 and 95                    */
  cout << " computing percentiles    ... ";
  
  vector < num >
    vcsort (vcurr);
  
  make_heap ( vcsort.begin(), vcsort.end() );
  sort_heap ( vcsort.begin(), vcsort.end() );
  
  num
    int_min = vcsort [ (unsigned) ( floor ( 0.05 * vcurr.size() ) ) ],
    int_max = vcsort [ (unsigned) ( floor ( 0.95 * vcurr.size() ) ) ];
  
  cout << " done \n";
  
  /* put the map voxels at the right place in the mask                         */
  cout << " resizing                 ... ";
  
  vcurr.resize(volsize, 0.);
  
  for (int i=(vcsort.size()-1); i>=0; i--) {
    
    vcurr[mask[i]]=vcurr[i];
    vcurr[i]=0;

  }
  
  cout << " done \n";
  
  /* write vcurr to the output image                                           */
  cout << " writing                  ... ";
  
  outputimage->dim[4] = 1;                      /* volume not a time series    */
  outputimage->cal_min = (float) (int_min);     /* lower bound for viewing     */
  outputimage->cal_max = (float) (int_max);     /* upper bound for viewing     */
  outputimage->datatype = DT_DOUBLE;            /* 64-bit float for storage    */
  nifti_datatype_sizes( outputimage->datatype,
			&outputimage->nbyper,
			&outputimage->swapsize);
  setNiftiBricks( outputimage, &vcurr );        /* cast/export nifti data      */
  nifti_update_dims_from_array ( outputimage ); /* fix the nifti record        */
  nifti_set_filenames ( outputimage, outputfilename.c_str(), 0, 0); /* write   */
  nifti_image_write ( outputimage );
  
  cout << " done \n";
  
  /* free memory                                                               */
  nifti_image_free ( inputimage );
  nifti_image_free ( outputimage );
  
  return 0;
  
} /* main                                                                      */




