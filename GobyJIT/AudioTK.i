%module AudioTK
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "FX/ATK.hpp"
#include "FX/ATKAdaptiveFilters.hpp"
#include "FX/ATKAttackRelease.hpp"
#include "FX/ATKAttackReleaseHysterisis.hpp"
#include "FX/ATKBesselFilters.hpp"
#include "FX/ATKBlockLMSFilter.hpp"
#include "FX/ATKButterworthFilters.hpp"
#include "FX/ATKChamberlinFilter.hpp"
#include "FX/ATKChebyshev1Filter.hpp"
#include "FX/ATKChebyshev2Filter.hpp"
#include "FX/ATKDelays.hpp"
#include "FX/ATKDistortionProcessors.hpp"
#include "FX/ATKDynamicProcessors.hpp"
#include "FX/ATKEqProcessors.hpp"
#include "FX/ATKExpander.hpp"
#include "FX/ATKFeedbackDelayNetwork.hpp"
#include "FX/ATKFIRFilter.hpp"
#include "FX/ATKFixedDelayLine.hpp"
#include "FX/ATKFollowerTransistorClassAProcessor.hpp"
#include "FX/ATKGainColoredCompressor.hpp"
#include "FX/ATKGainColoredExpander.hpp"
#include "FX/ATKGainCompressor.hpp"
//#include "FX/ATKGainFilter.hpp"
#include "FX/ATKGainLimiter.hpp"
#include "FX/ATKGainMaxColoredExpander.hpp"
#include "FX/ATKGainSwell.hpp"
#include "FX/ATKIIRFilter.hpp"
#include "FX/ATKLinkwitzReillyFilters.hpp"
#include "FX/ATKLMSFilter.hpp"
#include "FX/ATKMaxCompressor.hpp"
#include "FX/ATKMaxExpander.hpp"
#include "FX/ATKMultipleUniversalFixedDelayLine.hpp"

#include "FX/ATKPowerFilter.hpp"
//#include "FX/ATKPreampProcessors.hpp"
#include "FX/ATKRBJFilters.hpp"
#include "FX/ATKRelativePowerFilter.hpp"
#include "FX/ATKRemezFilter.hpp"

#include "FX/ATKReverbProcessors.hpp"
#include "FX/ATKRIAAFilters.hpp"
#include "FX/ATKRLSFilter.hpp"
#include "FX/ATKSecondOrderFilters.hpp"
#include "FX/ATKSecondOrderSVFFilters.hpp"
#include "FX/ATKTimeVaryingFilters.hpp"
#include "FX/ATKTimeVaryingSVFFilters.hpp"
#include "FX/ATKToneFilters.hpp"
#include "FX/ATKToneStackFilter.hpp"
#include "FX/ATKToolProcessors.hpp"
#include "FX/ATKTransistorClassAProcessor.hpp"
#include "FX/ATKTriode2Processor.hpp"
#include "FX/ATKTriodeProcessor.hpp"
#include "FX/ATKUniversalFixedDelayLine.hpp"
#include "FX/ATKUniversalVariableDelayLine.hpp"
#include "FX/ATKVariableDelayLine.hpp"


%}
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

%include "FX/ATK.hpp"
%include "FX/ATKAttackRelease.hpp"
%include "FX/ATKAttackReleaseHysterisis.hpp"
%include "FX/ATKBesselFilters.hpp"
%include "FX/ATKBlockLMSFilter.hpp"
%include "FX/ATKButterworthFilters.hpp"
%include "FX/ATKChamberlinFilter.hpp"
%include "FX/ATKChebyshev1Filter.hpp"
%include "FX/ATKChebyshev2Filter.hpp"
%include "FX/ATKExpander.hpp"
%include "FX/ATKFeedbackDelayNetwork.hpp"
%include "FX/ATKFIRFilter.hpp"
%include "FX/ATKFixedDelayLine.hpp"
%include "FX/ATKFollowerTransistorClassAProcessor.hpp"
%include "FX/ATKGainColoredCompressor.hpp"
%include "FX/ATKGainColoredExpander.hpp"
%include "FX/ATKGainCompressor.hpp"
//%include "FX/ATKGainFilter.hpp"
%include "FX/ATKGainLimiter.hpp"
%include "FX/ATKGainMaxColoredExpander.hpp"
%include "FX/ATKGainSwell.hpp"
%include "FX/ATKIIRFilter.hpp"
%include "FX/ATKLinkwitzReillyFilters.hpp"
%include "FX/ATKLMSFilter.hpp"
%include "FX/ATKMaxCompressor.hpp"
%include "FX/ATKMaxExpander.hpp"
%include "FX/ATKMultipleUniversalFixedDelayLine.hpp"
%include "FX/ATKPowerFilter.hpp"
%include "FX/ATKPreampProcessors.hpp"
%include "FX/ATKRBJFilters.hpp"
%include "FX/ATKRelativePowerFilter.hpp"
//%include "FX/ATKRemezFilter.hpp"
%include "FX/ATKReverbProcessors.hpp"
%include "FX/ATKRIAAFilters.hpp"
%include "FX/ATKRLSFilter.hpp"
%include "FX/ATKSecondOrderFilters.hpp"
%include "FX/ATKSecondOrderSVFFilters.hpp"
%include "FX/ATKTimeVaryingFilters.hpp"
%include "FX/ATKTimeVaryingSVFFilters.hpp"
%include "FX/ATKToneFilters.hpp"
%include "FX/ATKToneStackFilter.hpp"
%include "FX/ATKToolProcessors.hpp"
%include "FX/ATKTransistorClassAProcessor.hpp"
%include "FX/ATKTriode2Processor.hpp"
%include "FX/ATKTriodeProcessor.hpp"
%include "FX/ATKUniversalFixedDelayLine.hpp"
%include "FX/ATKUniversalVariableDelayLine.hpp"
%include "FX/ATKVariableDelayLine.hpp"


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
