%module Pyvc_tests.lua
%{
typedef float DspFloatType;
#include "SoundObject.hpp"
#include "Analog/vc_tests.lua.hpp"
%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"


typedef float DspFloatType;
%include "SoundObject.hpp"
%include "Analog/vc_tests.lua.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

