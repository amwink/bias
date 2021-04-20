# file: ventrmask.py
#
# making a standard hull of a (standard space) brain
# pre:  the file standard.nii is present in the current 
# directory
# post: the file standardhull.nii contans the hull of
# the brain
#
# this image is thresholded from the outside inwards
# until the first above-threshold voxel is found, i.e.,
# it does not remove voxels inside the hull.
#
# (c) Alle Meije Wink 2021
# a.m.wink@gmail.com

# libraries for file/nifity I/O and numerics
import nibabel as nib;
import numpy as np;
import itertools;
import os;

# load the standard space image 
N  = nib.load ( os.path.abspath( os.path.curdir ) + os.path.sep + "standard.nii.gz" );
n  = N.get_fdata();
sn = n.shape;
m  = np.zeros ( sn, dtype=np.int8 );

thresh = 3300;

# fill from inferior to superior
for x, y in itertools.product( range( sn [ 0 ] ), range( sn [ 1 ] ) ):
    z = 0;
    while ( ( z < sn[2]-1 ) & ( m[x][y][z] >= 0 ) & ( n[x][y][z] < thresh ) ):
        m[x][y][z] = -1;
        z += 1;
        
# fill from superior to inferior
for x, y in itertools.product( range( sn [ 0 ] ), range( sn [ 1 ] ) ):
    z = sn[2] - 1;
    while ( ( z >= 0 ) & ( m[x][y][z] >= 0 ) & ( n[x][y][z] < thresh ) ):
        m[x][y][z] = -1;
        z -= 1;

# fill from anterior to posterior
for x, z in itertools.product( range( sn [ 0 ] ), range( sn [ 2 ] ) ):
    y = 0;
    while ( ( y < sn[1] ) & ( m[x][y][z] >= 0 ) & ( n[x][y][z] < thresh ) ):
        m[x][y][z] = -1;
        y += 1;
        
# fill from posterior to anterior
for x, z in itertools.product( range( sn [ 0 ] ), range( sn [ 2 ] ) ):
    y = sn[1] - 1;
    while ( ( y >= 0 ) & ( m[x][y][z] >= 0 ) & ( n[x][y][z] < thresh ) ):
        m[x][y][z] = -1;
        y -= 1;

# fill from left to righ
for y, z in itertools.product( range( sn [ 1 ] ), range( sn [ 2 ] ) ):
    x = 0;
    while ( ( x < sn[0] ) & ( m[x][y][z] >= 0 ) & ( n[x][y][z] < thresh ) ):
        m[x][y][z] = -1;
        x += 1;
        
# fill from right to left
for y, z in itertools.product( range( sn [ 1 ] ), range( sn [ 2 ] ) ):
    x = sn[0] - 1;
    while ( ( x >= 0 ) & ( m[x][y][z] >= 0 ) & ( n[x][y][z] < thresh ) ):
        m[x][y][z] = -1;
        x -= 1;
        
Nout = nib.Nifti1Image( m + 1, N.affine );
Nout.to_filename ( os.path.dirname( N.get_filename() ) + os.path.sep + "standardhull.nii.gz" );
