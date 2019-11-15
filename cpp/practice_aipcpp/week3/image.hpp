#ifndef IMAGE_HPP_INCLUDED
#define IMAGE_HPP_INCLUDED

#include <vector>
#include <random>
#include <numeric>
#include <initializer_list>

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
    }; // class



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

        // += operators for scalar and image2d
        image2d& operator+= ( const value_type rhs ) {
                for (size_t s = 0; s < data.size(); s++ )
                    data[s] += rhs.data[s];
                return *this;
        }
        template <typename U>
        image2d& operator+= ( const image2d<U> rhs ) {
                for (size_t s = 0; s < data.size(); s++ )
                    data[s] += rhs.data[s];
                return *this;
        }

        // *= operators for scalar and image2d
        image2d& operator*= ( const value_type rhs ) {
                for (size_t s = 0; s < data.size(); s++ )
                    data[s] *= rhs;
                return *this;
        }        template <typename U>
        image2d& operator*= ( const image2d<U> rhs ) {
                for (size_t s = 0; s < data.size(); s++ )
                    data[s] *= rhs.data[s];
                return *this;
        }

        // -= operators for scalar and image2d
        image2d& operator-= ( const value_type rhs ){
                for (size_t s = 0; s < data.size(); s++ )
                    data[s] -= rhs;
                return *this;
        }
        template <typename U>
        image2d& operator-= ( const image2d<U> rhs ) {
                for (size_t s = 0; s < data.size(); s++ )
                    data[s] -= rhs.data[s];
                return *this;
        }

        // /= operators for scalar and image2d
        image2d& operator/= ( const value_type rhs )
                { for (size_t s = 0; s < data.size(); s++ )
                      data[s] /= rhs;
                  return *this;
        }
        template <typename U>
        image2d& operator/= ( const image2d<U> rhs ) {
                for (size_t s = 0; s < data.size(); s++ )
                    data[s] /= rhs.data[s];
                return *this;
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

}; // class


    // week 3

    template <typename T>
    class imageNd {



        public:

        typedef
            T value_type;



        // data are private -- only accessible via member functions
        private:

        std::vector <size_t>
            sizes;

        std::vector <size_t>
            strides;

        nifti_image
            *header = NULL;

        std::vector <value_type>
            data;




        public:

        // destructor: clear vectors and free pointer too nifti_image
        ~imageNd() {

            data.resize    (0);
            sizes.resize   (0);
            strides.resize (0);
            if ( header != NULL )
                free (header);

        }

        // deep copy constructor
        imageNd ( const imageNd& rhs ) {

                header = nifti_copy_nim_info (rhs.header);

                std::copy ( rhs.data.begin(),    rhs.data.end(),    data.begin()    );
                std::copy ( rhs.sizes.begin(),   rhs.sizes.end(),   sizes.begin()   );
                std::copy ( rhs.strides.begin(), rhs.strides.end(), strides.begin() );

        }

        // assignment operator for imageNd
        const imageNd<T>& operator= ( imageNd <T>& rhs ) {

            if ( this != &rhs ) {

                // just to make sure we don't leave stuff
                if ( this->header != NULL ) {
                    free ( header );
                    header = NULL;
                }

                header = nifti_copy_nim_info( rhs.header );

                // just to make sure we don't leave stuff
                data.resize    (0);
                sizes.resize   (0);
                strides.resize (0);

                std::copy ( rhs.data.begin(),    rhs.data.end(),    data.begin()    );
                std::copy ( rhs.sizes.begin(),   rhs.sizes.end(),   sizes.begin()   );
                std::copy ( rhs.strides.begin(), rhs.strides.end(), strides.begin() );

            } // if this != rhs

            return *this;

        }

        // constructor that reads an n-dimensional image ( for NIfTI 1 <= n <= 7 )
        imageNd ( std::string filename ) {

            header = nifti_image_read ( filename.c_str(),
                                        DO_READ_DATA );

            sizes.resize   ( header -> dim[0]   );
            strides.resize ( header -> dim[0]+1 );

            // make the array 'strides' so that it uses the last dimension as well
            strides [0] = 1;
            for ( size_t i=1; i<=sizes.size(); i++ ) {
                sizes [i-1] = header -> dim [i];
                strides [i] = strides[i-1] * sizes[i-1];
            }

            data.resize ( *(strides.rbegin()) ); // the end of strides holds the image's size

            aip::getNiftiBricks ( header,
                                  header -> data,
                                  data.size(),
                                  &data );

        }



        // operator() for positional addressing
        value_type const operator() ( std::initializer_list < size_t > const& indices ) const {
            size_t const offset =
                std::inner_product ( indices.begin(), indices.end(),
                                     strides.begin(),
                                     0 );
                return data [ offset ];
        }

        value_type& operator() ( std::initializer_list < size_t > const& indices ) {
            size_t const offset =
                std::inner_product ( indices.begin(), indices.end(),
                                     strides.begin(),
                                     0 );
                return data [ offset ] ;
        }



        // operators += for scalar and imageNd
        const imageNd<T>& operator+= ( const value_type& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] += rhs;
            return (*this);
        }
        template <typename U>
        const imageNd<T>& operator+= ( const imageNd<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] += rhs.data[s];
            return (*this);
        }

        // operator + for scalar and imageNd => value_type and imageNd are now templated
        template <typename U>
        const imageNd<T>& operator+ ( const U& rhs ) {
            imageNd out(*this);
            out += rhs;
            return out;
        }



        // operators *= for scalar and imageNd
        const imageNd<T>& operator*= ( const value_type& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] *= rhs;
            return (*this);
        }
        template <typename U>
        const imageNd<T>& operator*= ( const imageNd<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] *= rhs.data[s];
            return (*this);
        }

        // operator * for scalar and imageNd => value_type and imageNd are now templated
        template <typename U>
        const imageNd<T>& operator* ( const U& rhs ) {
            imageNd out(*this);
            out *= rhs;
            return out;
        }



        // operators -= for scalar and imageNd
        const imageNd<T>& operator-= ( const value_type rhs ){
            for (size_t s = 0; s < data.size(); s++ )
                data[s] -= rhs;
            return (*this);
        }
        template <typename U>
        const imageNd<T>& operator-= ( const imageNd<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] -= rhs.data[s];
            return (*this);
        }

        // operator - for scalar and imageNd => value_type and imageNd are now templated
        template <typename U>
        const imageNd<T>& operator- ( const U& rhs ) {
            imageNd out(*this);
            out -= rhs;
            return out;
        }



        // operators /= for scalar and imageNd
        const imageNd<T>& operator/= ( const value_type rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] /= rhs;
            return (*this);
        }
        template <typename U>
        const imageNd<T>& operator/= ( const imageNd<U>& rhs ) {
            for (size_t s = 0; s < data.size(); s++ )
                data[s] /= rhs.data[s];
            return (*this);
        }

        // operator / for scalar and imageNd => value_type and imageNd are now templated
        template <typename U>
        const imageNd<T>& operator/ ( const U& rhs ) {
            imageNd<T> out(*this);
            out /= rhs;
            return out;
        }



        // return an image with 1/c for all coefficients c
        const imageNd<T>& reciprocal ( const imageNd<T>& rhs ) {
            imageNd<T> out(*rhs);
            for (size_t s = 0; s < out.data.size(); s++ )
                out.data[s] = 1 / out.data[s];
            return out;
        }



        std::vector<size_t> getsize (          )
            { return sizes;       }

        size_t              getsize ( size_t s )
            { return sizes[s];    }

        value_type* getdata()
            { return data.data(); }

        void saveNII ( std::string filename ) {
            strcpy ( header -> fname, filename.c_str() );
            setNiftiBricks( header, &data );
            nifti_image_write( header );
        }



        // get a slice from a volume
        std::vector<value_type> getSlice( size_t dim, size_t sli ) {

            std::vector <value_type> ovec;

            if (sizes.size() != 3) {

                std::cout << "slices can only be selected from a 3D image\n";

            } else {

                value_type* dptr = NULL;
                value_type* optr = NULL;

                switch (dim) {

                case (0):

                    ovec.resize( sizes[1]*sizes[2] ); //output vector
                    optr = ovec.data(); // pointer in ovec

                    dptr = &data [ sli ]; // pointer in imageNd

                    for (size_t c=0; c<ovec.size(); c++, dptr+=sizes[0])
                        ovec[c] = *dptr; // copy x-slice
                    break;

                case (1):

                    ovec.resize ( sizes[0] * sizes[2] ); //output vector

                    dptr = &data [ sli * strides[1] ]; // pointer in imageNd
                    optr = ovec.data(); // pointer in ovec

                    for ( size_t y=0; y<sizes[2]; y++) {
                        for ( size_t x=0; x<sizes[0]; x++, optr++, dptr++ )
                            *optr = *dptr; // copy y-slice
                        dptr+=strides[2]; // jump forwards to next slice
                        dptr-=strides[1]; // jump backwards to start row
                    }

                    break;

                case (2):

                    ovec.resize ( sizes[0] * sizes[1] ); //output vector

                    dptr = &data [ sli * strides[2] ]; // pointer in imageNd
                    optr = ovec.data(); // pointer in ovec

                    for (size_t c=0; c<ovec.size(); c++)
                        *optr++ = *dptr++; // copy z-slice

                    break;

                }

            } // if sizes

            return ovec;

        } // getSlice


        void getSlice ( size_t dim, size_t sli, std::string filename ) {

            std::vector<value_type> slice ( getSlice( dim, sli) );

            if ( !slice.empty() ) {

                size_t slicex = (dim>0) ? 0 : 1;
                size_t slicey = (dim>1) ? 1 : 2;

                cimg_library::CImg<value_type>*
                my_bitmap = new cimg_library::CImg<value_type> (slice.data(),
                                                                sizes[slicex],
                                                                sizes[slicey],
                                                                1,
                                                                1,
                                                                true);

                my_bitmap->rotate(180);
                my_bitmap->save_bmp( filename.c_str() );

                delete ( my_bitmap );

            } // if !empty

        } // getSlice



        void addNormalNoise ( double mu, double sigma ) {

            // random device class instance, source of 'true' randomness for initializing random seed
            std::random_device randomdevice{};

            // Mersenne twister PRNG, initialized with seed from previous random device instance
            std::mt19937 engine { randomdevice() };

            // set up the distribution
            std::normal_distribution <double> normaldist ( mu, sigma );

            // add noise to the data
            for ( size_t i = 0; i < data.size(); i++ )
                data[i] += normaldist ( engine );

        }



}; // class

// we also want the +, *, -, / operators to work from right to left
// which is a bit more work for the (non-commutative) - and /
template <typename U, typename T>
inline imageNd <T> operator+ ( U x, imageNd <T> y) { return y             + x; }
template <typename U, typename T>
inline imageNd <T> operator* ( U x, imageNd <T> y) { return y             * x; }
template <typename U, typename T>
inline imageNd <T> operator- ( U x, imageNd <T> y) { return -y            + x; }
template <typename U, typename T>
inline imageNd <T> operator/ ( U x, imageNd <T> y) { return reciprocal(y) * x; }


}; // namespace

#endif // IMAGE_HPP_INCLUDED
