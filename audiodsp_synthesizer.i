%module Synthesizer
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

#include "AaronStaticMusicLib.hpp"
#include "AMSynth.hpp"
#include "AMSynthEffects.hpp"
#include "Arpeggiator.hpp"
#include "Chords.hpp"
#include "EGProcessors.hpp"
#include "FMProcessor.hpp"
//#include "libMTSClient.cpp
//#include "libMTSClient.h
//#include "LorisSynthesizer.cpp
//#include "LorisSynthesizer.hpp"
#include "ObXd.hpp"
#include "SF2Processors.hpp"
#include "StepSequencer.hpp"

%}

typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"


%include "AaronStaticMusicLib.hpp"
%include "AMSynth.hpp"
%include "AMSynthEffects.hpp"
%include "Arpeggiator.hpp"
%include "Chords.hpp"
%include "EGProcessors.hpp"
%include "FMProcessor.hpp"
//%include "libMTSClient.cpp"
//%include "libMTSClient.h"
//%include "LorisSynthesizer.cpp
//%include "LorisSynthesizer.hpp"
%include "ObXd.hpp"
//%include "SF2Processors.hpp"
%include "StepSequencer.hpp"


%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

/*
%inline %{
    const int BufferSize = 256;
    Default noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}
*/