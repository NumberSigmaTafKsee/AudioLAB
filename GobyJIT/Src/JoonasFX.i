%module JoonasFX
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

#include "FX/KHDelayReverb.hpp"
#include "FX/KHCrossDelay.hpp"
#include "FX/KHDelay.hpp"
#include "FX/KHMoorerReverb.hpp"
#include "FX/KHPingPongDelay.hpp"
#include "FX/KHSchreoderImproved.hpp"
#include "FX/KHSchreoderReverb.hpp"
#include "FX/KHSyncedTapDelayLine.hpp"

%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

%include "FX/AudioDSP_Chorus.hpp"
%include "FX/AudioDSP_DCBlocker.hpp"
%include "FX/AudioDSP_Delay.hpp"
%include "FX/AudioDSP_DelayLine.hpp"
%include "FX/AudioDSP_FirstOrderAllPass.hpp"
%include "FX/AudioDSP_Lfo.hpp"
%include "FX/AudioDSP_ModDelay.hpp"
%include "FX/AudioDSP_Phaser.hpp"
%include "FX/AudioDSP_VibraFlange.hpp"


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

