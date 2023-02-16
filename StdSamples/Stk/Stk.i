%module Stk
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>
   
#include "Synthesizer/StkHeaders.hpp"
#include "Synthesizer/StkADSRProcessor.hpp"
#include "Synthesizer/StkAsympProcessor.hpp"
#include "Synthesizer/StkBiquadProcessor.hpp"
#include "Synthesizer/StkBlit.hpp"
#include "Synthesizer/StkDrummer.hpp"
#include "Synthesizer/StkEnvelopeProcessor.hpp"
#include "Synthesizer/StkEnvelopeProcessors.hpp"
#include "Synthesizer/StkFilterProcessors.hpp"
#include "Synthesizer/StkFirProcessor.hpp"
#include "Synthesizer/StkFMBeeThree.hpp"
#include "Synthesizer/StkFMHevyMetl.hpp"
#include "Synthesizer/StkFMPercFlute.hpp"
#include "Synthesizer/StkFMProcessors.hpp"
#include "Synthesizer/StkFMRhodey.hpp"
#include "Synthesizer/StkFMTubeBell.hpp"
#include "Synthesizer/StkFMVoices.hpp"
#include "Synthesizer/StkFMWurley.hpp"
#include "Synthesizer/StkFormSweep.hpp"
#include "Synthesizer/StkFXChorus.hpp"
#include "Synthesizer/StkFXCubic.hpp"
#include "Synthesizer/StkFXDelay.hpp"
#include "Synthesizer/StkFXDelayA.hpp"
#include "Synthesizer/StkFXDelayL.hpp"
#include "Synthesizer/StkFXEcho.hpp"
#include "Synthesizer/StkFXFreeVerb.hpp"
#include "Synthesizer/StkFXJCRev.hpp"
#include "Synthesizer/StkFXLentPitchShift.hpp"
#include "Synthesizer/StkFXNRev.hpp"
#include "Synthesizer/StkFXPitchShift.hpp"
#include "Synthesizer/StkFXPRCRev.hpp"
#include "Synthesizer/StkFXProcessors.hpp"
#include "Synthesizer/StkFXTapDelay.hpp"
#include "Synthesizer/StkGenerators.hpp"
#include "Synthesizer/StkGranulate.hpp"
#include "Synthesizer/StkIirProcessor.hpp"
#include "Synthesizer/StkModulate.hpp"
#include "Synthesizer/StkMoog.hpp"
#include "Synthesizer/StkNoise.hpp"
#include "Synthesizer/StkOnePoleProcessor.hpp"
#include "Synthesizer/StkOneZeroProcessor.hpp"
#include "Synthesizer/StkPhysicalModels.hpp"
#include "Synthesizer/StkPoleZeroProcessor.hpp"
#include "Synthesizer/StkResonate.hpp"
#include "Synthesizer/StkSimple.hpp"
#include "Synthesizer/StkSineWave.hpp"
#include "Synthesizer/StkSingWave.hpp"
#include "Synthesizer/StkTwoPoleProcessor.hpp"
#include "Synthesizer/StkTwoZeroProcessor.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"
//%include "Synthesizer/StkHeaders.hpp"
%include "Synthesizer/StkADSRProcessor.hpp"
%include "Synthesizer/StkAsympProcessor.hpp"
%include "Synthesizer/StkBiquadProcessor.hpp"
%include "Synthesizer/StkBlit.hpp"
%include "Synthesizer/StkDrummer.hpp"
%include "Synthesizer/StkEnvelopeProcessor.hpp"
//%include "Synthesizer/StkEnvelopeProcessors.hpp"
//%include "Synthesizer/StkFilterProcessors.hpp"
%include "Synthesizer/StkFirProcessor.hpp"
%include "Synthesizer/StkFMBeeThree.hpp"
%include "Synthesizer/StkFMHevyMetl.hpp"
%include "Synthesizer/StkFMPercFlute.hpp"
%include "Synthesizer/StkFMProcessors.hpp"
%include "Synthesizer/StkFMRhodey.hpp"
%include "Synthesizer/StkFMTubeBell.hpp"
%include "Synthesizer/StkFMVoices.hpp"
%include "Synthesizer/StkFMWurley.hpp"
%include "Synthesizer/StkFormSweep.hpp"
%include "Synthesizer/StkFXChorus.hpp"
%include "Synthesizer/StkFXCubic.hpp"
%include "Synthesizer/StkFXDelay.hpp"
%include "Synthesizer/StkFXDelayA.hpp"
%include "Synthesizer/StkFXDelayL.hpp"
%include "Synthesizer/StkFXEcho.hpp"
%include "Synthesizer/StkFXFreeVerb.hpp"
%include "Synthesizer/StkFXJCRev.hpp"
%include "Synthesizer/StkFXLentPitchShift.hpp"
%include "Synthesizer/StkFXNRev.hpp"
%include "Synthesizer/StkFXPitchShift.hpp"
%include "Synthesizer/StkFXPRCRev.hpp"
%include "Synthesizer/StkFXProcessors.hpp"
%include "Synthesizer/StkFXTapDelay.hpp"
%include "Synthesizer/StkGenerators.hpp"
%include "Synthesizer/StkGranulate.hpp"
///%include "Synthesizer/StkHeaders.hpp"
%include "Synthesizer/StkIirProcessor.hpp"
%include "Synthesizer/StkModulate.hpp"
%include "Synthesizer/StkMoog.hpp"
%include "Synthesizer/StkNoise.hpp"
%include "Synthesizer/StkOnePoleProcessor.hpp"
%include "Synthesizer/StkOneZeroProcessor.hpp"
%include "Synthesizer/StkPhysicalModels.hpp"
%include "Synthesizer/StkPoleZeroProcessor.hpp"
%include "Synthesizer/StkResonate.hpp"
%include "Synthesizer/StkSimple.hpp"
%include "Synthesizer/StkSineWave.hpp"
%include "Synthesizer/StkSingWave.hpp"
%include "Synthesizer/StkTwoPoleProcessor.hpp"
%include "Synthesizer/StkTwoZeroProcessor.hpp"

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
