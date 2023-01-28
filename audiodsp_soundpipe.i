%module Soundpipe
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
#include "SoundPipe.hpp"
#include "SoundPipeADSR.hpp"
#include "SoundPipeAllPass.hpp"
#include "SoundPipeATone.hpp"
#include "SoundPipeAutoWah.hpp"
#include "SoundPipeBal.hpp"
#include "SoundPipeBar.hpp"
#include "SoundPipeBiquad.hpp"
#include "SoundPipeBiScale.hpp"
#include "SoundPipeBitCrush.hpp"
#include "SoundPipeBlSaw.hpp"
#include "SoundPipeBlSquare.hpp"
#include "SoundPipeBlTriangle.hpp"
#include "SoundPipeBrownNoise.hpp"
#include "SoundPipeButterworths.hpp"
#include "SoundPipeClip.hpp"
#include "SoundPipeClock.hpp"
#include "SoundPipeCombFilter.hpp"
#include "SoundPipeCompressor.hpp"
#include "SoundPipeConvolve.hpp"
#include "SoundPipeCount.hpp"
#include "SoundPipeCrossfade.hpp"
#include "SoundPipeDCBlock.hpp"
#include "SoundPipeDelay.hpp"
#include "SoundPipeDiode.hpp"
#include "SoundPipeDiskIn.hpp"
#include "SoundPipeDistortion.hpp"
#include "SoundPipeDrip.hpp"
#include "SoundPipeDTrig.hpp"
#include "SoundPipeDust.hpp"
#include "SoundPipeEQFilter.hpp"
#include "SoundPipeExponential.hpp"
#include "SoundPipeFOF.hpp"
#include "SoundPipeFOFilt.hpp"
#include "SoundPipeFog.hpp"
#include "SoundPipeMetronome.hpp"

%}

typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"


%include "SoundPipe.hpp"
%include "SoundPipeADSR.hpp"
%include "SoundPipeAllPass.hpp"
%include "SoundPipeATone.hpp"
%include "SoundPipeAutoWah.hpp"
%include "SoundPipeBal.hpp"
%include "SoundPipeBar.hpp"
%include "SoundPipeBiquad.hpp"
%include "SoundPipeBiScale.hpp"
%include "SoundPipeBitCrush.hpp"
%include "SoundPipeBlSaw.hpp"
%include "SoundPipeBlSquare.hpp"
%include "SoundPipeBlTriangle.hpp"
%include "SoundPipeBrownNoise.hpp"
%include "SoundPipeButterworths.hpp"
%include "SoundPipeClip.hpp"
%include "SoundPipeClock.hpp"
%include "SoundPipeCombFilter.hpp"
%include "SoundPipeCompressor.hpp"
%include "SoundPipeConvolve.hpp"
%include "SoundPipeCount.hpp"
%include "SoundPipeCrossfade.hpp"
%include "SoundPipeDCBlock.hpp"
%include "SoundPipeDelay.hpp"
%include "SoundPipeDiode.hpp"
%include "SoundPipeDiskIn.hpp"
%include "SoundPipeDistortion.hpp"
%include "SoundPipeDrip.hpp"
%include "SoundPipeDTrig.hpp"
%include "SoundPipeDust.hpp"
%include "SoundPipeEQFilter.hpp"
%include "SoundPipeExponential.hpp"
%include "SoundPipeFOF.hpp"
%include "SoundPipeFOFilt.hpp"
%include "SoundPipeFog.hpp"
%include "SoundPipeMetronome.hpp"


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