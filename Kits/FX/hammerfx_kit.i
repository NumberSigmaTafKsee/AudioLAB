%module HammerFX
%{
#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

typedef float DspFloatType;
#include "SoundObject.hpp"

//#include "FX/HammerFX.cpp
//#include "FX/HammerFX.hpp"
#include "FX/HammerFXAutoWah.hpp"
#include "FX/HammerFXChorus.hpp"
#include "FX/HammerFXDelay.hpp"
#include "FX/HammerFXDistortion.hpp"
#include "FX/HammerFXEcho.hpp"
#include "FX/HammerFXPhaser.hpp"
#include "FX/HammerFXPitch.hpp"
#include "FX/HammerFXRotary.hpp"
#include "FX/HammerFXSustain.hpp"
#include "FX/HammerFXTremolo.hpp"
#include "FX/HammerFXTubeAmp.hpp"
#include "FX/HammerFXVibrato.hpp"
%}

typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"



%include "FX/HammerFX.hpp"
%include "FX/HammerFXAutoWah.hpp"
%include "FX/HammerFXChorus.hpp"
%include "FX/HammerFXDelay.hpp"
%include "FX/HammerFXDistortion.hpp"
%include "FX/HammerFXEcho.hpp"
%include "FX/HammerFXPhaser.hpp"
%include "FX/HammerFXPitch.hpp"
%include "FX/HammerFXRotary.hpp"
%include "FX/HammerFXSustain.hpp"
%include "FX/HammerFXTremolo.hpp"
%include "FX/HammerFXTubeAmp.hpp"
%include "FX/HammerFXVibrato.hpp"


%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%inline %{
    const int BufferSize = 256;
    Std::RandomMersenne noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}