%module af
%{
#include <cstdint>
#include "af.h"
using namespace ArrayFire;
%}

%include "std_complex.i"

%template(float_complex) std::complex<float>;
%template(double_complex) std::complex<double>;

typedef long long   dim_t;
#define AFDLL

%include<af/defines.h>


#undef AFAPI
#define AFAPI 
#undef AF_DEPRECATED
#define AF_DEPRECATED(x)

%include<af/constants.h>
%include<af/complex.h>
%include<af/dim4.hpp>
%include<af/index.h>
%include<af/seq.h>
%include<af/util.h>
%include<arrayfire.h>
%include<af/random.h>
%include<af/algorithm.h>
%include<af/arith.h>
%include<af/random.h>
%include<af/blas.h>
%include<af/features.h>
%include<af/graphics.h>
%include<af/image.h>
%include<af/lapack.h>
%include<af/ml.h>
%include<af/signal.h>
%include<af/sparse.h>
%include<af/statistics.h>
%include<af/vision.h>

%include "af.h"


%template(float_array) ArrayFire::Array<float>;
%template(double_array) ArrayFire::Array<double>;

%template(complex32) std::complex<float>;
%template(complex64) std::complex<double>;

%template(complex32_array) ArrayFire::Array<af::cfloat>;
%template(complex64_array) ArrayFire::Array<af::cdouble>;

// missing af_write
//%template(af_bool) ArrayFire::Array<bool,b8>;
%template(uint8_array) ArrayFire::Array<uint8_t>;
%template(int16_Array) ArrayFire::Array<int16_t>;
%template(uint16_Array) ArrayFire::Array<uint16_t>;
%template(int32_Array) ArrayFire::Array<int32_t>;
%template(uint32_Array) ArrayFire::Array<uint32_t>;
// not sure why int64_t is a problem.
%template(int64_Array)  ArrayFire::Array<long long int>;
%template(uint64_Array) ArrayFire::Array<unsigned long long int>;
//%template(HalfArray)   ArrayFire::Array<float,fp16>;

%template(float_vector) ArrayFire::Vector<float>;
%template(double_vector) ArrayFire::Vector<double>;

%template(float_matrix) ArrayFire::Matrix<float>;
%template(double_matrix) ArrayFire::Matrix<double>;

%inline %{
    void  set_float(float * p, size_t i, float v) { p[i] = v; }
    float get_float(float * p, size_t i) { return p[i]; }
%}
%extend ArrayFire::Array {
    double __getitem(size_t i) { return (*$self)[ArrayFire::Dim4(i); ]}
    void   __setitem(size_t i, double val) { (*$self)[ArrayFire::Dim4(i)] = val; }
}

%template(float_scalar) ArrayFire::Scalar<float>;
%template(double_scalar) ArrayFire::Scalar<double>;
/*
%template(char_scalar) ArrayFire::Scalar<char>;
%template(uchar_scalar) ArrayFire::Scalar<unsigned char>;
%template(short_scalar) ArrayFire::Scalar<short>;
%template(ushort_scalar) ArrayFire::Scalar<unsigned short>;
%template(int_scalar) ArrayFire::Scalar<int>;
%template(uint_scalar) ArrayFire::Scalar<unsigned int>;
%template(long_scalar) ArrayFire::Scalar<long>;
%template(ulong_scalar) ArrayFire::Scalar<unsigned long>;
%template(llong_scalar) ArrayFire::Scalar<long long>;
%template(ullong_scalar) ArrayFire::Scalar<unsigned long long>;
*/
%template(float_image) ArrayFire::ImageProcessing<float>;
%template(double_image) ArrayFire::ImageProcessing<double>;

%template(float_vision) ArrayFire::ComputerVision<float>;
%template(double_vision) ArrayFire::ComputerVision<double>;

%template(float_signal) ArrayFire::SignalProcessing<float>;
%template(double_signal) ArrayFire::SignalProcessing<double>;