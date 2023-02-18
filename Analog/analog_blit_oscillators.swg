%module analog_blit_oscillators
%{
typedef float DspFloatType;
#include "GenericSoundObject.hpp"
#include "analog_blit_oscillators.hpp"
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
%include "analog_blit_oscillators.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%template(Blit2SawOscillator) Analog::Oscillators::Blit2SawOscillator<DspFloatType>;
%template(Blit2SquareOscillator) Analog::Oscillators::Blit2SquareOscillator<DspFloatType>;
%template(Blit2TriangleOscillator) Analog::Oscillators::Blit2TriangleOscillator<DspFloatType>;
