#ifdef PYTHON
%module PyAnalogBlitOscillators
#else
%module analog_blit_oscillators
#endif
%{
#define DSPFLOATDOUBLE    
#include "GenericSoundObject.hpp"
#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "analog_blit_oscillators.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
#ifdef LUA
%include "lua_fnptr.i"
#endif

#define DSPFLOATDOUBLE
%include "GenericSoundObject.hpp"
%include "analog_blit_oscillators.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%template(Blit2SawOscillator_f32) Analog::Oscillators::Blit2SawOscillator<float>;
%template(Blit2SawOscillator_f64) Analog::Oscillators::Blit2SawOscillator<double>;
%template(Blit2SquareOscillator_f32) Analog::Oscillators::Blit2SquareOscillator<float>;
%template(Blit2SquareOscillator_f64) Analog::Oscillators::Blit2SquareOscillator<double>;
%template(Blit2TriangleOscillator_f32) Analog::Oscillators::Blit2TriangleOscillator<float>;
%template(Blit2TriangleOscillator_f64) Analog::Oscillators::Blit2TriangleOscillator<double>;