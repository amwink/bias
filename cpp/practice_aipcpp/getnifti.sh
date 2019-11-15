rm -rf nifticlib*
wget http://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
tar zxvf nifticlib-2.0.0.tar.gz
ln -sf nifticlib-2.0.0 nifticlib
mkdir -p nifti/include
mkdir -p nifti/src
cd nifti/include
ln -sf ../../nifticlib/niftilib/nifti1.h .
ln -sf ../../nifticlib/niftilib/nifti1_io.h
cd ../src
ln -sf ../../nifticlib/niftilib/nifti1_io.c
cd ../..
mkdir -p znz/include
mkdir -p znz/src
cd znz/include 
ln -sf ../../nifticlib/znzlib/znzlib.h
cd ../src
ln -sf ../../nifticlib/znzlib/znzlib.c
cd ../..
