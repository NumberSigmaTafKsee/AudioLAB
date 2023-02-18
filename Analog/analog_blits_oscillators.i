%module analog_blits_oscillators
%{
#include "GenericSoundObject.hpp"
#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>


#include "analog_blit_saw_oscillator.hpp"
#include "analog_blit_square_oscillator.hpp"
#include "analog_blit_triangle_oscillator.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "GenericSoundObject.hpp"
%include "analog_blit_saw_oscillator.hpp"
%include "analog_blit_square_oscillator.hpp"
%include "analog_blit_triangle_oscillator.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%template(BlitSaw) Analog::Oscillators::BlitSaw<DspFloatType>;
%template(BlitSquare) Analog::Oscillators::BlitSquare<DspFloatType>;
%template(BlitTriangle) Analog::Oscillators::BlitTriangle<DspFloatType>;
