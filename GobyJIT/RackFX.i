%module RackFX
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


//#include "FX/racka3.cpp
//#include "FX/racka3.hpp"
//#include "FX/RackFX.hpp"
#include "FX/RackFXAlienWah.hpp"
#include "FX/RackFXAnalogPhaser.hpp"
#include "FX/RackFXArpie.hpp"
#include "FX/RackFXBeatTracker.hpp"
#include "FX/RackFXChorus.hpp"
#include "FX/RackFXCoilCrafter.hpp"
#include "FX/RackFXCompBand.hpp"
#include "FX/RackFXCompressor.hpp"
#include "FX/RackFXConvolotron.hpp"
#include "FX/RackFXDFlange.hpp"
#include "FX/RackFXDynamicFilter.hpp"
#include "FX/RackFXEcho.hpp"
#include "FX/RackFXEchotron.hpp"
#include "FX/RackFXEqualizer.hpp"
#include "FX/RackFXExciter.hpp"
#include "FX/RackFXExpander.hpp"
#include "FX/RackFXFormantFilter.hpp"
#include "FX/RackFXGate.hpp"
#include "FX/RackFXHarmEnhancer.hpp"
#include "FX/RackFXHarmonizer.hpp"
#include "FX/RackFXInfinity.hpp"
#include "FX/RackFXLooper.hpp"
#include "FX/RackFXMBDist.hpp"
#include "FX/RackFXMBVol.hpp"
#include "FX/RackFXMetronome.hpp"
#include "FX/RackFXMusicDelay.hpp"
#include "FX/RackFXNewDIst.hpp"
#include "FX/RackFXOpticalTrem.hpp"
#include "FX/RackFXPan.hpp"
#include "FX/RackFXPhaser.hpp"
#include "FX/RackFXPitchShifter.hpp"
#include "FX/RackFXRBEcho.hpp"
#include "FX/RackFXRBFilter.hpp"
//#include "FX/RackFXRecChord.cpp
#include "FX/RackFXRecognize.hpp"
#include "FX/RackFXResample.hpp"
#include "FX/RackFXReverb.hpp"
#include "FX/RackFXSustainer.hpp"
#include "FX/RackFXSVFilter.hpp"
#include "FX/RackFXSynthFilter.hpp"
#include "FX/RackFXValve.hpp"
#include "FX/RackFXVibe.hpp"
#include "FX/RackFXVocoder.hpp"
#include "FX/RackFXWaveshaper.hpp"
%}

typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

//%include "FX/RackFX.hpp"
%include "FX/RackFXAlienWah.hpp"
%include "FX/RackFXAnalogPhaser.hpp"
%include "FX/RackFXArpie.hpp"
%include "FX/RackFXBeatTracker.hpp"
%include "FX/RackFXChorus.hpp"
%include "FX/RackFXCoilCrafter.hpp"
%include "FX/RackFXCompBand.hpp"
%include "FX/RackFXCompressor.hpp"
%include "FX/RackFXConvolotron.hpp"
%include "FX/RackFXDFlange.hpp"
%include "FX/RackFXDynamicFilter.hpp"
%include "FX/RackFXEcho.hpp"
%include "FX/RackFXEchotron.hpp"
%include "FX/RackFXEqualizer.hpp"
%include "FX/RackFXExciter.hpp"
%include "FX/RackFXExpander.hpp"
%include "FX/RackFXFormantFilter.hpp"
%include "FX/RackFXGate.hpp"
%include "FX/RackFXHarmEnhancer.hpp"
%include "FX/RackFXHarmonizer.hpp"
%include "FX/RackFXInfinity.hpp"
%include "FX/RackFXLooper.hpp"
%include "FX/RackFXMBDist.hpp"
%include "FX/RackFXMBVol.hpp"
%include "FX/RackFXMetronome.hpp"
%include "FX/RackFXMusicDelay.hpp"
%include "FX/RackFXNewDIst.hpp"
%include "FX/RackFXOpticalTrem.hpp"
%include "FX/RackFXPan.hpp"
%include "FX/RackFXPhaser.hpp"
%include "FX/RackFXPitchShifter.hpp"
%include "FX/RackFXRBEcho.hpp"
%include "FX/RackFXRBFilter.hpp"
%include "FX/RackFXRecognize.hpp"
%include "FX/RackFXResample.hpp"
%include "FX/RackFXReverb.hpp"
%include "FX/RackFXSustainer.hpp"
%include "FX/RackFXSVFilter.hpp"
%include "FX/RackFXSynthFilter.hpp"
%include "FX/RackFXValve.hpp"
%include "FX/RackFXVibe.hpp"
%include "FX/RackFXVocoder.hpp"
%include "FX/RackFXWaveshaper.hpp"


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
