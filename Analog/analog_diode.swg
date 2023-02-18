%module analog_diode
%{
typedef float DspFloatType;
#include "GenericSoundObject.hpp"
#include "analog_diode.hpp"
%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

typedef float DspFloatType;
%include "GenericSoundObject.hpp"
%include "analog_diode.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%template(Diode) Analog::Distortion::Diode::Dioder<DspFloatType>;
