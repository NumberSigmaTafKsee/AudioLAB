%module Envelopes
%{
typedef float DspFloatType;
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>


#include "FX/ADSR.hpp"
#include "FX/ADSR2.hpp"
#include "FX/GammaEnvelope.hpp"
//#include "FX/qmADSR.hpp"
#include "FX/CEnvelopeDetector.hpp"


%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

%include "FX/ADSR.hpp"
%include "FX/ADSR2.hpp"
%include "FX/GammaEnvelope.hpp"
//%include "FX/qmADSR.hpp"
%include "FX/CEnvelopeDetector.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%inline %{
    const int BufferSize = 256;
    Default noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}
