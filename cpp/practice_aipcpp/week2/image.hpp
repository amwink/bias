#ifndef IMAGE_HPP_INCLUDED
#define IMAGE_HPP_INCLUDED

#include <vector>

// downloaded from https://framagit.org/dtschump/CImg/raw/master/CImg.h
#include "CImg/CImg.h"

// using this for NIfTI file I/O
#include "niftiio.hpp"

// used for NIfTI I/O only do read is relevant now
enum readdata {
  NO_READ_DATA,
  DO_READ_DATA
};

namespace aip {

    // week 1

    class floatimage {

        public:

            typedef float value_type;

        private:

            const int
                height = 256,        // fixed width
                width  = 256;        // fixed height

            value_type               // not the best of names but
                data [256 * 256];    // for NIfTI IO compatibility

        public:

            floatimage() {

                for (int i=0; i < (this->width * this->height); i++ )
                    this->data [i] = 0;

            }

            const int
                getwidth()  { return this->width; }
            const int
                getheight() { return this->width; }

            void
                setpixel ( int position, value_type value ) {
                    this -> data [ position ] = value;
                }

            const value_type*
                getdata () { return &this->data[0]; }
    };

    // week 2

    template <typename T>
    class image2d {

        public:

        typedef
            T value_type;

        // data are private -- only accessible via member functions
        private:

        size_t
            width  = 0,
            height = 0;
        std::vector <value_type>
            data;

        public:

        // constructors
        image2d ();

        // copy-constructor
        image2d ( image2d <value_type> const &rhs ) :
                width  ( rhs.width ),
                height ( rhs.height ),
                data   ( rhs.data ) {}

        // construct from old-fashioned C-style array
        image2d (  value_type *dataref,
                    size_t w, size_t h) :
                    width ( w ), height ( h ) {

                            data.assign ( dataref,
                            dataref + ( w * h ) );

            }

        // constructor that reads a 2D slice from a 3D volume
        image2d ( std::string filename, size_t slice ) {

            nifti_image
                *inputimage = nifti_image_read ( filename.c_str(),
                                                 DO_READ_DATA );
            size_t
                step      = inputimage -> dim[1];

            width  = inputimage -> dim[2];
            height = inputimage -> dim[3];
            data.resize ( width * height );

            vector<value_type> datatmp ( width * height * step );

            aip::getNiftiBricks ( inputimage,
                                  inputimage -> data,
                                  datatmp.size(),
                                  &datatmp );

            value_type *v = &datatmp[ slice ];

            for ( size_t s = 0;
                  s < ( width * height );
                  s ++ )
                data[s] = *(v += step);

        }

        // destructor
        // ~image2d();

        // operators [] for rows and columns
        value_type const* operator[] (size_t const r) const
                   { return &data [ r * width ]; }
        value_type*       operator[] (size_t const r)
                   { return &data [ r * width ]; }

        // assignment operator for image2d
        const image2d<T> operator= ( image2d <T>& rhs ) {
            image2d<T> out(rhs);
            width  = rhs.getwidth();
            height = rhs.getheight();
            data   = rhs.data; // std::vector handles vector assignment
        }

        // operators += for scalar and image2d
        const image2d<T>& operator+= ( const value_type rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] += rhs.data[s];
            return (*this);
        }
        template <typename U>
        const image2d<T>& operator+= ( const image2d<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] += rhs.data[s];
            return (*this);
        }

        // operators + for scalar and image2d
        const image2d<T> operator+ ( const value_type rhs ) {
            image2d out(*this);
            out += rhs;
            return out;
        }
        const image2d<T> operator+ ( const image2d<T>& rhs ) {
            image2d out(*this);
            out += rhs;
            return out;
        }



        // operators *= for scalar and image2d
        const image2d<T>& operator*= ( const value_type rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] *= rhs;
            return (*this);
        }
        template <typename U>
        const image2d<T>& operator*= ( const image2d<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] *= rhs.data[s];
            return (*this);
        }

        // operators * for scalar and image2d
        const image2d<T> operator* ( const value_type rhs ) {
            image2d out(*this);
            out *= rhs;
            return out;
        }
        const image2d<T> operator* ( const image2d<T>& rhs ) {
            image2d out(*this);
            out *= rhs;
            return out;
        }



        // operators -= for scalar and image2d
        const image2d<T>& operator-= ( const value_type rhs ){
            for (size_t s = 0; s < data.size(); s++ )
                data[s] -= rhs;
            return (*this);
        }
        template <typename U>
        const image2d<T>& operator-= ( const image2d<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] -= rhs.data[s];
            return (*this);
        }

        // operators - for scalar and image2d
        const image2d<T> operator- ( const value_type rhs ) {
            image2d out(*this);
            out -= rhs;
            return out;
        }
        const image2d<T> operator- ( const image2d<T>& rhs ) {
            image2d out(*this);
            out -= rhs;
            return out;
        }



        // operators /= for scalar and image2d
        const image2d<T>& operator/= ( const value_type rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] /= rhs;
            return (*this);
        }
        template <typename U>
        const image2d<T>& operator/= ( const image2d<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] /= rhs.data[s];
            return (*this);
        }

        // operators / for scalar and image2d
        const image2d<T> operator/ ( const value_type rhs ) {
            image2d out(*this);
            out /= rhs;
            return out;
        }
        const image2d<T> operator/ ( const image2d<T>& rhs ) {
            image2d out(*this);
            out /= rhs;
            return out;
        }



        size_t getwidth  ()   { return width;       }
        size_t getheight ()   { return height;      }
        value_type* getdata() { return data.data(); }

        // other functions
        void saveBMP ( std::string filename ) {

            cimg_library::CImg<value_type> my_bitmap(getdata(),
                                                width,
                                                height,
                                                1,
                                                1,
                                                true);

            my_bitmap.save_bmp( filename.c_str() );

        }

};

};

#endif // IMAGE_HPP_INCLUDED
