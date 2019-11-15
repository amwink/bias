#include <iostream>
#include <utility>
#include <cmath>

// our own image class
#include "image.hpp"

int main() {

    // read slice 85
    aip::image2d<float> my2dimageA ( "../data/brainT1.nii.gz", 85 );

    // read slice 95 and add it to slice 85 using operator+()
    aip::image2d<float> my2dimageB ( "../data/brainT1.nii.gz", 95 );
    my2dimageA += my2dimageB;

    // make the average by dividing by the appropriate number (slices you have used)
    my2dimageA /= 2;

    // make a square of white pixels (wit intensity 255)
    // of width 100 and the upper left corner at (77,77)
    // using the class's 2-dimensional brackets indexing
    size_t
        zz = 77,
        hz = zz + 100;

    // draw 4 lines in 1 for-loop
    for ( size_t i = 0; i <= 100; i++ ) {
        my2dimageA[zz]  [zz+i] = 255.;
        my2dimageA[zz+i][zz]   = 255.;
        my2dimageA[hz]  [zz+i] = 255.;
        my2dimageA[zz+i][hz]   = 255.;
    }

    // draw a largest triangle ** possible inside the
    // square that shares at most one side with it
    //
    // ** there are more than 1 possibilities (right answers)
    for (size_t i = 0; i < 100; i++) {
        my2dimageA[zz+i][zz+50-size_t(round(i/2))] = 255.;
        my2dimageA[zz+i][zz+50+size_t(round(i/2))] = 255.;
    }

    // draw the largest circle possible inside the square
    for (size_t i = 0; i < 71; i++) {
        size_t distance = round ( sqrt ( 2500 - i*i/4 ) );
        my2dimageA[zz+50+distance][zz+50+i/2] = 255.;
        my2dimageA[zz+50+distance][zz+50-i/2] = 255.;
        my2dimageA[zz+50-distance][zz+50+i/2] = 255.;
        my2dimageA[zz+50-distance][zz+50-i/2] = 255.;
        my2dimageA[zz+50+i/2][zz+50+distance] = 255.;
        my2dimageA[zz+50+i/2][zz+50-distance] = 255.;
        my2dimageA[zz+50-i/2][zz+50+distance] = 255.;
        my2dimageA[zz+50-i/2][zz+50-distance] = 255.;
    }

    // flip upside down so that the head is upright in the BMP
    for ( size_t i = 0; i < (my2dimageA.getwidth() / 2); i++ )
        for ( size_t j = 0; j < my2dimageA.getheight() ; j++ )
            std::swap ( my2dimageA[i][j], my2dimageA[my2dimageA.getwidth() -i -1][j] );

    // write to BMP
    my2dimageA.saveBMP ( "../data/brainT1_week2.bmp" );

    return ( 0 );

}
