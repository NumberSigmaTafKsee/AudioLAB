%module VAOscillators
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VAOscillators.hpp"
#include "Analog/VAPolyBLEPOscillator.hpp"
#include "Analog/VAPolyBlepOscillators.hpp"
#include "Analog/VAMinBlepOscillators.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

%rename Oscillators::BlitSaw        VABlitSawOsc;
%rename Oscillators::BlitSquare     VABlitSquareOsc;
%rename Oscillators::BlitTriangle   VABlitTriangleOsc;
%rename Oscillators::BlitDSF        VABlitDSFOsc;
%rename Oscillators::blitSaw        VABlitSawOsc2;
%rename Oscillators::blitSquare     VABlitSquareOsc2;
%rename Oscillators::blitTriangle   VABlitTriangleOsc2;
%rename Oscillators::DPWSaw         VADPWSawOsc;
%rename Oscillators::DPWPulse       VADPWPulseOsc;
%rename Oscillators::DPWTriangle    VADPWTriangleOsc;;
%include "Analog/VAOscillators.hpp"

%ignore Analog::Oscillators::PolyBLEPOsc::blep;
%ignore Analog::Oscillators::PolyBLEPOsc::blamp;
%include "Analog/VAPolyBLEPOscillator.hpp"
//%include "Analog/VAPolyBlepOscillators.hpp"
%include "FX/minBLEP.hpp"

/*
#include "Analog/VABlitOscillators.hpp"
#include "Analog/VADPWOscillators.hpp"
#include "Analog/VAOscillators.hpp"
#include "Analog/VAPolyBlepOscillators.hpp"
*/


%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

