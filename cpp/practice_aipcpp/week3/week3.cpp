#include <iostream>
#include <utility>

// our own image class
#include "image.hpp"

int main() {

    // read a 3D image
    aip::imageNd<float> myNdimageA ( "../data/brainT1.nii.gz" );

    // make a sphere of radius 50 around voxel 91,128,128
    // by setting the intensities on its surface to 255
    for ( size_t h = 0; h <= 50; h++ ) {

        // h is the distance from slice 91;
        // the sphere is a cricle in this slice
        // its radius is given by sqrt ( 2500 - h*h )
        float radius = sqrt ( 2500. - h*h );

        // als we need to do now is draw a circle with
        // size |radius| in the 2 slices with distance
        // h from slice 91
        size_t slices [2] = { 91 - h, 91 + h };

        for (auto s : slices) {

            float circlesegment = 1.42 * radius; // just above sqrt(2)
            float csquare = radius * radius;

            for (size_t i = 0; i < circlesegment; i++) {

                size_t distance = round ( sqrt ( csquare - i*i/4. ) );

                myNdimageA( {s, 128-distance, 128-i/2} ) = 0.;
                myNdimageA( {s, 128-distance, 128+i/2} ) = 0.;
                myNdimageA( {s, 128+distance, 128-i/2} ) = 0.;
                myNdimageA( {s, 128+distance, 128+i/2} ) = 0.;
                myNdimageA( {s, 128-i/2, 128-distance} ) = 0.;
                myNdimageA( {s, 128-i/2, 128+distance} ) = 0.;
                myNdimageA( {s, 128+i/2, 128-distance} ) = 0.;
                myNdimageA( {s, 128+i/2, 128+distance} ) = 0.;

            } // for i (points on the circle in that slice

        } // for s (two slices equally far from 91)

    } // for h (distance to slice 91)

    // add lots of normally distributed noise
    myNdimageA.addNormalNoise( 25, 5 );

    // save as another image
    myNdimageA.saveNII( "../data/football.nii.gz" );

    // get X slice 91  and write as BMP
    // get Y slice 127 and write as BMP
    // get Z slice 160 and write as BMP
    myNdimageA.getSlice( 0, 91,  "../data/xslice.bmp" );
    myNdimageA.getSlice( 1, 127, "../data/yslice.bmp" );
    myNdimageA.getSlice( 2, 127, "../data/zslice.bmp" );

    return ( 0 );

}
