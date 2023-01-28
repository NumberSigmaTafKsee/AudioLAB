%module analog_diode_clipper
%{
#define DSPFLOATDOUBLE
#include "SoundObject.hpp"
#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "analog_diode_clipper.hpp"
using namespace Analog::Distortion::Diode;
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"


#define DSPFLOATDOUBLE
%include "SoundObject.hpp"
%include "analog_diode_clipper.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%template(DiodeClipper) Analog::Distortion::Diode::DiodeClipper<DspFloatType>;