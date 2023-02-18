%module analog_diode_ladder_filter
%{
typedef float DspFloatType;
#include "GenericSoundObject.hpp"
#include "analog_diode_ladder_filter.hpp"
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
%include "analog_diode_ladder_filter.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%template(DiodeLadderFilter) Analog::Filters::DiodeLadderFilter2::DiodeLadderFilter<DspFloatType>;
