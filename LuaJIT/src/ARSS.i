%module ARSS
%{
extern "C" {
#include "dsp.h"
#include "image_io.h"
#include "sound_io.h"
}
%}

%include "stdint.i"
%include "std_vector.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(ldouble_vector) std::vector<long double>;
%template(char_vector) std::vector<signed char>;
%template(uchar_vector) std::vector<unsigned char>;
%template(short_vector) std::vector<signed short>;
%template(ushort_vector) std::vector<unsigned short>;
%template(int_vector) std::vector<signed int>;
%template(uint_vector) std::vector<unsigned int>;
%template(long_vector) std::vector<signed long>;
%template(ulong_vector) std::vector<unsigned long>;
%template(llong_vector) std::vector<signed long long int>;
%template(ullong_vector) std::vector<unsigned long long int>;





%include "dsp.h"
%include "image_io.h"
%include "sound_io.h"

