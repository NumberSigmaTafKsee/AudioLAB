%module audiodsp_pink_noise
%{
typedef float DspFloatType;
#include "SoundObject.hpp"
#include "audiodsp_pink_noise.hpp"
%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

typedef float DspFloatType;
%include "SoundObject.hpp"
%include "audiodsp_pink_noise.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

