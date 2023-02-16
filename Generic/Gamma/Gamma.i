%module Gamma
%{

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "SoundObject.hpp"

//#include "Synthesizer/Gamma.hpp"
//#include "Synthesizer/GammaAnalysis.hpp"
#include "Synthesizer/GammaAnalysisInspectorProcessor.hpp"
#include "Synthesizer/GammaAnalysisMaxAbsProcessor.hpp"
#include "Synthesizer/GammaAnalysisPCounterProcessor.hpp"
#include "Synthesizer/GammaAnalysisProcessors.hpp"
#include "Synthesizer/GammaAnalysisThresholdProcessor.hpp"
#include "Synthesizer/GammaAnalysisZeroCrossProcessor.hpp"
#include "Synthesizer/GammaAnalysisZeroCrossRateProcessor.hpp"
//#include "Synthesizer/GammaAudio.hpp"
#include "Synthesizer/GammaCombFilterProcessors.hpp"
#include "Synthesizer/GammaDelayProcessors.hpp"
//#include "Synthesizer/GammaDelays.hpp"
#include "Synthesizer/GammaDFTProcessors.hpp"
//#include "Synthesizer/GammaEnvelopeProcessors.hpp"
//#include "Synthesizer/GammaEnvelopes.hpp"
#include "Synthesizer/GammaEnvelopesADProcessor.hpp"
#include "Synthesizer/GammaEnvelopesADSRProcessor.hpp"
#include "Synthesizer/GammaEnvelopesDecayProcessor.hpp"
#include "Synthesizer/GammaEnvelopesGateProcessor.hpp"
#include "Synthesizer/GammaEnvelopesSegmentProcessors.hpp"
#include "Synthesizer/GammaEnvFollowerProcessor.hpp"
#include "Synthesizer/GammaFFTProcessors.hpp"
//#include "Synthesizer/GammaFilterProcessors.hpp"
//#include "Synthesizer/GammaFilters.hpp"
#include "Synthesizer/GammaFiltersAllpass1.hpp"
#include "Synthesizer/GammaFiltersAllpass2.hpp"
#include "Synthesizer/GammaFiltersBiquad.hpp"
#include "Synthesizer/GammaFiltersBlockNyquist.hpp"
#include "Synthesizer/GammaFiltersDCBlock.hpp"
#include "Synthesizer/GammaFiltersDifferencer.hpp"
#include "Synthesizer/GammaFiltersHilbert.hpp"
#include "Synthesizer/GammaFiltersIntegrator.hpp"
#include "Synthesizer/GammaFiltersMovingAverage.hpp"
#include "Synthesizer/GammaFiltersNotch.hpp"
#include "Synthesizer/GammaFiltersOnePole.hpp"
#include "Synthesizer/GammaFiltersReson.hpp"
//#include "Synthesizer/GammaFX.hpp"
#include "Synthesizer/GammaFXAMProcessor.hpp"
#include "Synthesizer/GammaFXBiquad3.hpp"
#include "Synthesizer/GammaFXBurst.hpp"
#include "Synthesizer/GammaFXChebyshevPolynomials.hpp"
#include "Synthesizer/GammaFXChirp.hpp"
#include "Synthesizer/GammaFXChorus.hpp"
#include "Synthesizer/GammaFXFreqShift.hpp"
#include "Synthesizer/GammaFXMonoSynth.hpp"
#include "Synthesizer/GammaFXPanProcessor.hpp"
#include "Synthesizer/GammaFXPluck.hpp"
//#include "Synthesizer/GammaFXProcessors.hpp"
#include "Synthesizer/GammaFXQuantizer.hpp"
#include "Synthesizer/GammaNoiseProcessors.hpp"
#include "Synthesizer/GammaOscillatorProcessors.hpp"
//#include "Synthesizer/GammaOscillators.hpp"
#include "Synthesizer/GammaRecorderProcessor.hpp"
#include "Synthesizer/GammaSampleProcessors.hpp"
#include "Synthesizer/GammaSDFTProcessors.hpp"
#include "Synthesizer/GammaSilenceDetectorProcessor.hpp"
#include "Synthesizer/GammaSpatialProcessors.hpp"
#include "Synthesizer/GammaSTFTProcessors.hpp"
#include "Synthesizer/GammaTransferFunctions.hpp"
#include "Synthesizer/GammaVoiceProcessors.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

#

%include "SoundObject.hpp"
//%include "Synthesizer/Gamma.hpp"
//%include "Synthesizer/GammaAnalysis.hpp"
%include "Synthesizer/GammaAnalysisInspectorProcessor.hpp"
%include "Synthesizer/GammaAnalysisMaxAbsProcessor.hpp"
%include "Synthesizer/GammaAnalysisPCounterProcessor.hpp"
%include "Synthesizer/GammaAnalysisProcessors.hpp"
%include "Synthesizer/GammaAnalysisThresholdProcessor.hpp"
%include "Synthesizer/GammaAnalysisZeroCrossProcessor.hpp"
%include "Synthesizer/GammaAnalysisZeroCrossRateProcessor.hpp"
//%include "Synthesizer/GammaAudio.hpp"
%include "Synthesizer/GammaCombFilterProcessors.hpp"
%include "Synthesizer/GammaDelayProcessors.hpp"
//%include "Synthesizer/GammaDelays.hpp"
%include "Synthesizer/GammaDFTProcessors.hpp"
//%include "Synthesizer/GammaEnvelopeProcessors.hpp"
//%include "Synthesizer/GammaEnvelopes.hpp"
%include "Synthesizer/GammaEnvelopesADProcessor.hpp"
%include "Synthesizer/GammaEnvelopesADSRProcessor.hpp"
%include "Synthesizer/GammaEnvelopesDecayProcessor.hpp"
%include "Synthesizer/GammaEnvelopesGateProcessor.hpp"
%include "Synthesizer/GammaEnvelopesSegmentProcessors.hpp"
%include "Synthesizer/GammaEnvFollowerProcessor.hpp"
%include "Synthesizer/GammaFFTProcessors.hpp"
//%include "Synthesizer/GammaFilterProcessors.hpp"
//%include "Synthesizer/GammaFilters.hpp"
%include "Synthesizer/GammaFiltersAllpass1.hpp"
%include "Synthesizer/GammaFiltersAllpass2.hpp"
%include "Synthesizer/GammaFiltersBiquad.hpp"
%include "Synthesizer/GammaFiltersBlockNyquist.hpp"
%include "Synthesizer/GammaFiltersDCBlock.hpp"
%include "Synthesizer/GammaFiltersDifferencer.hpp"
%include "Synthesizer/GammaFiltersHilbert.hpp"
%include "Synthesizer/GammaFiltersIntegrator.hpp"
%include "Synthesizer/GammaFiltersMovingAverage.hpp"
%include "Synthesizer/GammaFiltersNotch.hpp"
%include "Synthesizer/GammaFiltersOnePole.hpp"
%include "Synthesizer/GammaFiltersReson.hpp"
//%include "Synthesizer/GammaFX.hpp"
%include "Synthesizer/GammaFXAMProcessor.hpp"
%include "Synthesizer/GammaFXBiquad3.hpp"
%include "Synthesizer/GammaFXBurst.hpp"
%include "Synthesizer/GammaFXChebyshevPolynomials.hpp"
%include "Synthesizer/GammaFXChirp.hpp"
%include "Synthesizer/GammaFXChorus.hpp"
%include "Synthesizer/GammaFXFreqShift.hpp"
%include "Synthesizer/GammaFXMonoSynth.hpp"
%include "Synthesizer/GammaFXPanProcessor.hpp"
%include "Synthesizer/GammaFXPluck.hpp"
//%include "Synthesizer/GammaFXProcessors.hpp"
%include "Synthesizer/GammaFXQuantizer.hpp"
%include "Synthesizer/GammaNoiseProcessors.hpp"
//%include "Synthesizer/GammaOscillatorProcessors.hpp"
%include "Synthesizer/GammaOscillators.hpp"
%include "Synthesizer/GammaRecorderProcessor.hpp"
%include "Synthesizer/GammaSampleProcessors.hpp"
%include "Synthesizer/GammaSDFTProcessors.hpp"
%include "Synthesizer/GammaSilenceDetectorProcessor.hpp"
%include "Synthesizer/GammaSpatialProcessors.hpp"
%include "Synthesizer/GammaSTFTProcessors.hpp"
%include "Synthesizer/GammaTransferFunctions.hpp"
%include "Synthesizer/GammaVoiceProcessors.hpp"
%include "Synthesizer/GranularProcessor.hpp"


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
