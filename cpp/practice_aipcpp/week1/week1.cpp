#include <iostream>

// downloaded from https://framagit.org/dtschump/CImg/raw/master/CImg.h
#include "CImg/CImg.h"

// our own image class
#include "image.hpp"

// using this for NIfTI file I/O
#include "niftiio.hpp"

// used for NIfTI I/O only do read is relevant now
enum readdata {
  NO_READ_DATA,
  DO_READ_DATA
};

int main() {

    // 4. initialise an image object
    //    make sure the constructor that
    //    initialises the pixels, is called
    aip::floatimage my_floatimage;

    // the following code loads the image
    //     ../data/brainT1.nii.gz
    //     (original http://www.jannin.org/mritemplate)
    //
    // it uses a mix of C functions from the
    // original nifti library and a modern
    // C++ interface from niftiio.hpp
    nifti_image
        *inputimage = nifti_image_read ( "../data/brainT1.nii.gz",
                                         DO_READ_DATA );
    unsigned char*
        pixeldata = (unsigned char *) ( inputimage->data );
    int
        slice_offset = 90;

    // 5. put 1 slice into our image
    //    what happens if you change slice_offset?
    //    what is the meaning of += 182?
    for (int i=0; i < (my_floatimage.getwidth() * my_floatimage.getheight()); i++ )
        my_floatimage.setpixel( i,
                                float ( pixeldata[ slice_offset += 182 ] ) );

    // this code writes the slice to a BMP file
    cimg_library::CImg<float> my_bitmap(my_floatimage.getdata(),
                                        256,
                                        256,
                                        1,
                                        1,
                                        true);
    my_bitmap.rotate(180);
    my_bitmap.save_bmp("../data/brainT1_week1.bmp");

    // 6. check the BMP file with an image viewer.
    //    is it what you expected?
    //    and if you uncomment the 'rotate' line?

    return ( 0 );

}
