#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 19:36:21 2020

@author: amwink
"""


# import general libraries os      (for handling files/directories), 
#                          time    (for reporting processing times),
#                          numpy   (for all the numerical routines),
#                          nibabel (for reading/writing NIfTI files)
import os;
import time;
import argparse;
import numpy as np;
import nibabel as nib;

# import routines from the dipy package (dipy.org)
# for performing 3D image registration / resamplng
from dipy.io.image import ( load_nifti, save_nifti );
from dipy.align.imaffine import ( transform_centers_of_mass, MutualInformationMetric, AffineRegistration );
from dipy.align.transforms import ( TranslationTransform3D, RigidTransform3D, AffineTransform3D );



################################################################################
################################################################################



def main():

    # parse options
    parser = argparse.ArgumentParser();
    parser.add_argument ( "-m", "--mrifiles",   help="text file with paths to NIfTIs of T1 scans", 
                          dest="mri_file", default="mr_images.txt" );
    parser.add_argument ( "-g", "--greymatter", help="text file with paths extracted GM density maps", 
                          dest="gm_file", default="gm_masks.txt" );
    parser.add_argument ( "-t", "--template",   help="path to MNI-space T1 template", 
                          default="/usr/share/fsl/data/standard/MNI152lin_T1_2mm.nii.gz" );
    args = parser.parse_args();

    # MNI-space template
    mni_t1_template = args.template;
    
    # load MR image filenames from a list    
    mri_file = args.mri_file;    
    with open( mri_file, 'r' ) as filehandle:    
        mri_files = [ current_place.rstrip() for current_place in filehandle.readlines() ];
     
    # load GM density map filenames from a list
    gm_file = args.gm_file;    
    with open( gm_file, 'r' ) as filehandle:
        gm_files = [ current_place.rstrip() for current_place in filehandle.readlines() ];

    # process them in pairs to build GM networks
    for mr, gm in zip ( mri_files, gm_files ):
        print ( 'processing MR scan and GM density map \n\t  {0:s}\n\t{1:s}'.format( mr, gm ) );

        start = time.process_time();
        obs_ran_net, net_file_name = gm_network ( mr, gm, mni_t1_template );
        print ( 'finished in {:.2f}s'.format( time.process_time() - start) );

        print ( "filename: {}  ".format ( net_file_name     ) );
        print ( "# points: {}\n".format ( obs_ran_net.shape ) );
        
        # not yet implemented
        # process_atlas_networks ( obs_ran_net, net_file_name, native_atlas );

################################################################################
################################################################################



def gm_network ( mr_filename, gm_filename, template_mr ):
       
    new_gmfilename = os.path.abspath (gm_filename).replace (".nii","_mni.nii");
    networks = 0;
    
    if not ( ( os.path.isfile ( new_gmfilename )  ) or ( os.path.islink ( new_gmfilename ) ) ):

        print ( 'file {} does not exist'.format ( new_gmfilename ) );
        print ( 'performing registration to MNI ... ', end = '' );
        start = time.process_time();
                
        # see https://dipy.org/documentation/1.2.0./examples_built/ ..
        #           .. affine_registration_3d/#example-affine-registration-3d
        static, static_affine = load_nifti( template_mr );      
        moving, moving_affine = load_nifti( mr_filename );       
        grey,   grey_affine   = load_nifti( gm_filename );

        # first initialise by putting centres of mass on top of each other
        c_of_mass  = transform_centers_of_mass ( static, static_affine, moving, moving_affine );

        # initialise transform parameters (e.g. the mutual information criterion)
        # these parameters won' need to be changed between the different stages
        nbins           = 64
        sampling_prop   = None
        metric          = MutualInformationMetric ( nbins, sampling_prop );
        level_iters     = [ 25, 15, 5 ];
        sigmas          = [ 2,  1,  0 ];
        factors         = [ 4,  2,  1 ];
        affreg          = AffineRegistration ( metric = metric, level_iters = level_iters, sigmas = sigmas, factors = factors );

        # give slightly more degrees of freedom, by allowing translation of centre of gravity
        print ( '\nTranslation only:' );
        transform   = TranslationTransform3D();
        params0     = None;
        translation = affreg.optimize ( static, moving, transform, params0, static_affine, moving_affine, starting_affine = c_of_mass.affine );

        # refine further by allowing all rigid transforms (rotations/translations around the centre of gravity)
        print ( 'Rigid transform:' );
        transform = RigidTransform3D();
        params0   = None;
        rigid     = affreg.optimize ( static, moving, transform, params0, static_affine, moving_affine, starting_affine = translation.affine );

        full_affine = False; # the GM networks method is based on keeping the cortical shape intact
        
        if ( full_affine ):

            # refine to a full affine transform by adding scaling and shearing
            print ( 'Affine transform:' );
            transform = AffineTransform3D();
            params0   = None;
            affine    = affreg.optimize(static, moving, transform, params0, static_affine, moving_affine, starting_affine = rigid.affine );
            final     = affine;
        
        else:   
            
            final     = rigid;
        
        # transform the grey matter data instead of the MRI itself       
        resampled = final.transform ( grey );
        save_nifti ( new_gmfilename, resampled, static_affine );
        
        print ( 'finished in {:.2f}s'.format( time.process_time() - start ) );
       
    if ( ( os.path.isfile ( new_gmfilename )  ) or ( os.path.islink ( new_gmfilename ) ) ):
        
        # only cube size implemented so far
        cubesize    = 3;
        
        # load the grey matter map and the template to which it was registered
        gm_img        = nib.load   ( new_gmfilename );        
        template_data = np.asarray ( nib.load ( template_mr ).dataobj );
        gm_data       = np.asarray ( gm_img.dataobj );        
        
        # find the best cube grid position (with the most nonzero cubes)
        cube_nonzeros, cube_offsets, gm_incubes = cube_grid_position ( gm_data, template_data, cubesize);
        gm_shape  = gm_incubes.shape;
   
        # write out the cube map, where each voxel in a cube is labelled with its cube index
        # be aware of the @ operator, this is a true matrix product A[n*m] x B[m*p] = C [n*p]
        cubes_data                  = np.zeros ( template_data.shape ).flatten();
        cubes_data [ cube_offsets ] = np.ones( cubesize ** 3) [ :, np.newaxis ] @ np.arange (cube_nonzeros).reshape(1,cube_nonzeros);
        cubes_data                  = cubes_data.reshape ( template_data.shape );
        cubes_file                  = os.path.abspath (gm_filename).replace (".nii","_cubes.nii");
        cubes_map                   = nib.Nifti1Image ( cubes_data, gm_img.affine );
        cubes_map.to_filename ( cubes_file );       

        # make a randomised version of the grey matter densities in the cubes
        # 1: exchange between and inside cubes (could be too many degrees of freedom!)
        gm_random = gm_incubes.flatten();
        gm_random = gm_random [ np.random.permutation ( len ( gm_random ) ).reshape ( gm_shape ) ];
        # 2: exchange cubes only ( this won't change the values in the correlation matrix, only positions )
        # gm_random = gm_incubes [ :, np.random.permutation ( gm_shape [1] ) ];
        # 3: exchange cubes and shuffle inside cubes
        # gm_random = gm_incubes [ np.random.permutation ( gm_shape [0] ), np.random.permutation ( gm_shape [1] )[ :, np.newaxis ] ]; 
        
        add_diag = True;

        # name of the NIfTI file with networks
        networks_file = os.path.abspath (gm_filename).replace (".nii","_gmnet.nii");

        if not ( ( os.path.isfile ( networks_file )  ) or ( os.path.islink ( networks_file ) ) ):

            # compute the cross correlation for observed and randomised cubes
            networks = cube_cross_correlation ( gm_incubes, gm_random, cubesize, add_diag );
    
            # save the networks to a file
            networks_map  = nib.Nifti1Image( networks, np.eye(4) );
            networks_map.to_filename ( networks_file );
            
        else:
            
            print ( "loading already existing file" );
            networks = np.asarray ( nib.load ( networks_file ).dataobj );        

    return networks, networks_file;



################################################################################
################################################################################



# compute the 1D positions of coefficients inside all cubes
# of a given grid size, for every 3D grid starting position
def cube_offsets ( template_sizes, cube_size ):
    
    # number of x, y and z positions in the template space
    dimx = template_sizes [0];
    dimy = template_sizes [1];   
    dimz = template_sizes [2];   
    
    # distance between offsets of one cube is the same everywhere
    # (independent of cube position) - given by dim[0] and dim[1]
    ind_1 =     np.arange ( cube_size ) ;
    ind_2 = ( ( np.arange ( cube_size ) * dimx)         [ :, np.newaxis ] + ind_1 ).flatten();
    ind_3 = ( ( np.arange ( cube_size ) * dimx * dimy ) [ :, np.newaxis ] + ind_2 ).flatten();    
    # print ( 'indices of cube at position (0,0,0) :\n{}'.format(ind_3) );

    # initialise the array of offset arrays (which may have different sizes)
    offsets = np.empty ( ( cube_size, cube_size, cube_size ), dtype = object );
    
    # now build all the offset lists and put them in
    for startx in range ( cube_size ):
        for starty in range ( cube_size ):
            for startz in range ( cube_size ):
                
                # indices found for an image of multiple slices slice: row * column * sli
                sliz = np.arange ( startz, dimz - cube_size, step = cube_size );
                isli = ( sliz * ( dimx * dimy ) ) [ :, np.newaxis ] + ind_3;
                isli = isli.flatten();
                                
                # indices found for an image of 1 slice: row * column
                coly = np.arange ( starty, dimy - cube_size, step = cube_size );
                icol = (coly * dimx) [ :, np.newaxis ] + isli;
                icol = icol.flatten();
                
                # indices found for an image of 1 row
                rowx = np.arange ( startx, dimx - cube_size, step = cube_size );
                irow = rowx [ :, np.newaxis ] + icol;
                irow = irow.flatten();
              
                offsets [startx] [starty] [startz] = irow;
                
    #for startx in range ( cube_size ):
    #    for starty in range ( cube_size ):
    #        for startz in range ( cube_size ):
    #            print ( 'size of cube grid {},{},{}: {}'.format(startx,starty,startz,len(offsets[startx][starty][startz]) ) );
                
    return offsets;
    
    
    
################################################################################
################################################################################



# given a grey matter density image, find the best cube grid position 
# given a grid size, which is the grid with the most non-empty cubes    
def cube_grid_position ( gm_data, template_data, cube_side ):

    # number of coefficients in a cube
    cubevol    = cube_side ** 3;    

    # get all 1D cube coefficient indices 
    # for every 3D grid starting position
    alloffsets = cube_offsets ( template_data.shape, cube_side );
    
    # we have not fount any non-empty cubes
    numgoodcubes=0;
    
    for posx in range ( alloffsets.shape[0] ):
        for posy in range ( alloffsets.shape[1] ):
            for posz in range ( alloffsets.shape[2] ):

                # use gm data at all cube offsets to start with
                offset = alloffsets[posx][posy][posz];
                cubes  = gm_data.flatten() [ offset ];
                
                # reshape (matlab / Fortran style) as matrix where each cube is a column
                offset = offset.reshape ( cubevol, len ( cubes ) // cubevol, order = 'F' );
                cubes  = cubes.reshape  ( cubevol, len ( cubes ) // cubevol, order = 'F' );
                
                # limit to cubes with nonzero variance
                cvar = np.var ( cubes, 0 ) ;     
                cvar = np.nonzero ( cvar > .001 )[0];                
                lcvar = len ( cvar );                

                # update is the current offset has more valid cubes
                if ( lcvar > numgoodcubes ):
                    numgoodcubes = lcvar;
                    greymatter = cubes  [ :, cvar ];
                    bestoffset = offset [ :, cvar ];                    

    return  numgoodcubes, bestoffset, greymatter;



################################################################################
################################################################################



def cube_rotations ( cube_size, diagonals ):
#
# Takes the dimensions of a cube as argument and returns the array
# a which contains the indices to rotate a cube of that dimension.
#
# An example for a square whose indices are 1..9:
#
# a = [ 0 1 2; ...
#       3 4 5; ...
#       6 7 8 ];
#
# For nearest neighbour interpolation, rotating the square by 45
# degrees corresponds to shuffling the indices:
#
#           2
#         1   5
# a45 = 0   4   8 
#         3   7 
#           6  
#           
# a45 = [ 1 2 5; ...
#         0 4 8; ...
#         3 6 7 ];
#
# Because this is algebraically correct; a45(a45) = a90, and so on,
# this corresponds to rotating around the Z axis. other rotations 
# are harder to visualise. 
#
# This function returns the new positions for all rotions that are
# multiples of 45 degrees around the X, Y and Z axes.
#
    if not ( cube_size == 3 ):
        print ( 'sorry, only implemented for 3x3x3 cubes ...' );
        rot_shifts = [];

    # original cube, in 3D and as a column vector
    rcube = np.arange ( 27 ).reshape ( 3, 3, 3 );
    vcube = rcube.flatten();
        
    # The example is rotation about the Z axis, top slice
    # (the 3 slices of the cube together have indices 0:26)
    zfl  = np.flip ( rcube, 0 ).flatten();
    z45  = np.asarray ( [ 1, 2, 5, 0, 4, 8, 3, 6, 7, 10, 11, 14, 9, 13, 17, 12, 15, 16, 19, 20, 23, 18, 22, 26, 21, 24, 25 ] );
    z90  = np.rot90 ( rcube, 1, axes = ( 1, 2 ) ).flatten();
    z135 = z45 [ z90  ];
    z180 = z90 [ z90  ];
    z225 = z45 [ z180 ];
    z270 = z90 [ z180 ];
    z315 = z45 [ z270 ];

    # rotation over x goes similarly
    yfl  = np.flip ( rcube, 1 ).flatten();
    y45  = np.asarray ( [ 21, 22, 23, 18, 19, 20, 9, 10, 11, 24, 25, 26, 12, 13, 14, 0, 1, 2, 15, 16, 17, 6, 7, 8, 3, 4, 5 ] );
    y90  = np.rot90 ( rcube, 1, axes = ( 0, 1 ) ).flatten();
    y135 = y45 [ y90  ];
    y180 = y90 [ y90  ];
    y225 = y45 [ y180 ];
    y270 = y90 [ y180 ];
    y315 = y45 [ y270 ];
    
    # rotation over x goes similarly
    xfl  = np.flip ( rcube, 2 ).flatten();
    x45  = np.asarray ( [ 1, 2, 11, 4, 5, 14, 7, 8, 17, 0, 10, 20, 3, 13, 23, 6, 16, 26, 9, 18, 19, 12, 21, 22, 15, 24, 25 ] );
    x90  = np.rot90 ( rcube, 1, axes = ( 0, 2 ) ).flatten();
    x135 = x45 [ x90  ];
    x180 = x90 [ x90  ];
    x225 = x45 [ x180 ];
    x270 = x90 [ x180 ];
    x315 = x45 [ x270 ];

    rot_shifts = np.vstack ( ( vcube, zfl, yfl, xfl, z90, z180, z270, y90, y180, y270, z90, z180, z270 ) );
    if ( diagonals ):
        rot_shifts = np.vstack ( ( rot_shifts, z45, z135, z225, z315, y45, y135, y225, y315, x45, x135, x225, x315 ) );

    return rot_shifts;



################################################################################
################################################################################



def cube_cross_correlation ( observed_columns, random_columns, cube_side, use_diagonals ):
    
    # each column is a cube    
    num_columns = observed_columns.shape[1];
    
    # compute different orderings of coefficients inside a cube after rotations
    angles     = cube_rotations ( cube_side, use_diagonals );
    num_angles = angles.shape [ 0 ];
        
    ori_mean = np.mean ( observed_columns, 0, keepdims = False );
    ori_mean = observed_columns - ori_mean;
    ori_std  = np.sqrt( np.sum ( np.square ( ori_mean ), 0, keepdims = False ) ); 
    
    ran_mean = np.mean (   random_columns, 0, keepdims = False );
    ran_mean = random_columns - ran_mean;
    ran_std  = np.sqrt( np.sum ( np.square ( ran_mean ), 0, keepdims = False ) );     

    # initialise the matrices for observed and permuted data
    # 'slice' 0 is from the observed data, slice 1 from the permutation
    networks = np.zeros ( ( num_columns, num_columns, 2 ), dtype='float' );

    for index_i in range ( num_columns ):
        
        c_ori   = ori_mean [ :, index_i ];
        s_ori   = ori_std [ index_i ];
        rot_ori = c_ori [ angles ];
        
        c_ran   = ran_mean [ :, index_i ];
        s_ran   = ran_std [ index_i ];        
        rot_ran = c_ran [ angles ];
        
        if ( not ( index_i % 100 ) ):
            print ( 'cube {}'.format( index_i ), end='\r' );
        
        for index_j in range ( index_i + 1, num_columns ):
            
            c_ori2  = ori_mean [ :, index_j ];
            s_ori2  = ori_std  [    index_j ];
            # be aware of the @ operator, this is a true matrix product A[n*m] x B[m*p] = C [n*p]
            cor_ori = ( rot_ori @ c_ori2 ) / ( ( s_ori * s_ori2 ) * np.ones ( num_angles ) ); 
            networks [ index_i, index_j, 0 ] = cor_ori.max();

            c_ran2  = ran_mean [ :, index_j ];          
            s_ran2  = ran_std  [    index_j ];            
            # be aware of the @ operator, this is a true matrix product A[n*m] x B[m*p] = C [n*p]
            cor_ran = ( rot_ran @ c_ran2 ) / ( ( s_ran * s_ran2 ) * np.ones ( num_angles ) );   
            networks [ index_i, index_j, 1 ] = cor_ran.max();
    
    return networks;
    


################################################################################
################################################################################



# if nothing else has been done yet, call main()    
if __name__ == '__main__': 
    main() 