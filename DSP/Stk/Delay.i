%module Delay
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

//#include "FX/ringbuffer.hpp"
//#include "FX/RingBuffer.hpp"
#include "FX/RingBufferProcessor.hpp"
//#include "FX/RingBuffers.hpp"
#include "FX/BasicDelay.hpp"
#include "FX/BasicDelayLine.hpp"
#include "FX/BasicDelayLineStereo.hpp"
//#include "FX/CircularDelay.hpp"
#include "FX/Comb.hpp"
#include "FX/CombFilter.hpp"
//#include "FX/CombFilters.hpp"
#include "FX/CrossDelayLine.hpp"
#include "FX/DelayFilters.hpp"
#include "FX/DelayLine.hpp"
#include "FX/DelayLines.hpp"
//#include "FX/DelayProcessors.hpp"
//#include "FX/Delays.hpp"
#include "FX/DelaySmooth.hpp"
#include "FX/DelaySyncedTapDelayLine.hpp"
#include "FX/Mu45.hpp"
//#include "FX/Mu45FilterCalc.hpp"
#include "FX/ERTapDelayLine.hpp"
#include "FX/Moorer.hpp"
//#include "FX/MoorerStereo.hpp"
#include "FX/PPDelayLine.hpp"
#include "FX/reverb.hpp"
#include "FX/StereoDelay.hpp"
#include "FX/Schroeder.hpp"
#include "FX/SchroederAllpass.hpp"
#include "FX/SchroederImproved.hpp"
//#include "FX/SchroederStereo.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

//%include "FX/ringbuffer.hpp"
//%include "FX/RingBuffer.hpp"
%include "FX/RingBufferProcessor.hpp"
//%include "FX/RingBuffers.hpp"
%include "FX/BasicDelay.hpp"
%include "FX/BasicDelayLine.hpp"
%include "FX/BasicDelayLineStereo.hpp"
//%include "FX/CircularDelay.hpp"
%include "FX/Comb.hpp"
%include "FX/CombFilter.hpp"
//%rename FX::CombFilters::CombFilter FXCombFilter;
//%include "FX/CombFilters.hpp"
%include "FX/CrossDelayLine.hpp"
%include "FX/DelayFilters.hpp"
%include "FX/DelayLine.hpp"
%include "FX/DelayLines.hpp"
//%include "FX/DelayProcessors.hpp"
//%include "FX/Delays.hpp"
%include "FX/DelaySmooth.hpp"
%include "FX/DelaySyncedTapDelayLine.hpp"
%include "FX/Mu45.hpp"
//%include "FX/Mu45FilterCalc.hpp"
%include "FX/ERTapDelayLine.hpp"
%include "FX/Moorer.hpp"
//%include "FX/MoorerStereo.hpp"
%include "FX/PPDelayLine.hpp"
%include "FX/reverb.hpp"
%include "FX/StereoDelay.hpp"
%include "FX/Schroeder.hpp"
%include "FX/SchroederAllpass.hpp"
%include "FX/SchroederImproved.hpp"
//%include "FX/SchroederStereo.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%inline %{
    int BufferSize = 256;
    Std::RandomMersenne noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;

    void Init(DspFloatType sr, int size=256) {
        sampleRate = sr;
        inverseSampleRate = 1.0/sr;
        invSampleRate = 1.0/sr;
        BufferSize = size;
        //noise.seed()
    }
%}
