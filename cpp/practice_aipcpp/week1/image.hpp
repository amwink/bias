#ifndef IMAGE_HPP_INCLUDED
#define IMAGE_HPP_INCLUDED

// 1. define a namespace for your own functions
//    so that you can use function names that
//    make sense without mixing existing code
namespace aip {

    // 2. define a 256x256 image class that has
    //    a) a container that contains <float> pixel data
    //    b) constant values for the width and height
    class floatimage {

        public:

            typedef float value_type;

        private:

            const int
                height = 256,        // fixed width
                width  = 256;        // fixed height

            value_type                    // not the best of names but
                data [256 * 256];    // for NIfTI IO compatibility

        public:
        // 3. create a default constructor that
        //    sets all voxels to zero
            floatimage() {

                for (int i=0; i < (this->width * this->height); i++ )
                    this->data [i] = 0;

            }

        // 4. create functions that return the
        //    width and height, respectively
            const int
                getwidth()  { return this->width; }
            const int
                getheight() { return this->width; }

        // 5. create a function that sets a data pixel
            void
                setpixel ( int position, value_type value ) {
                    this -> data [ position ] = value;
                }

        // 6. create a function that returns the pointer to the data
            const value_type*
                getdata () { return &this->data[0]; }
    };

};

#endif // IMAGE_HPP_INCLUDED
