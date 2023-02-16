%module ReissFX
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

#include "FX/AudioEffectAmplifier.hpp"
#include "FX/AudioEffectAutoWah.hpp"
#include "FX/AudioEffectChorus.hpp"
#include "FX/AudioEffectCompressor.hpp"
#include "FX/AudioEffectDelay.hpp"
#include "FX/AudioEffectDistortion.hpp"
#include "FX/AudioEffectEarlyDelay.hpp"
#include "FX/AudioEffectEarlyReflection.hpp"
#include "FX/AudioEffectFlanger.hpp"
#include "FX/AudioEffectHallReverb.hpp"
#include "FX/AudioEffectLFO.hpp"
#include "FX/AudioEffectParametricEQ.hpp"
#include "FX/AudioEffectPhaser.hpp"
#include "FX/AudioEffectPingPong.hpp"
#include "FX/AudioEffectPitchShifter.hpp"
#include "FX/AudioEffectPlateReverb.hpp"
#include "FX/AudioEffectRingMod.hpp"
#include "FX/AudioEffectRoomReverb.hpp"
#include "FX/AudioEffects.hpp"
#include "FX/AudioEffectsSuite.hpp"
#include "FX/AudioEffectStereoDelay.hpp"
#include "FX/AudioEffectSVF.hpp"
#include "FX/AudioEffectTremolo.hpp"
#include "FX/AudioEffectVibrato.hpp"
#include "FX/AudioProcessor.hpp"


%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

%include "FX/AudioEffectAmplifier.hpp"
%include "FX/AudioEffectAutoWah.hpp"
%include "FX/AudioEffectChorus.hpp"
%include "FX/AudioEffectCompressor.hpp"
%include "FX/AudioEffectDelay.hpp"
%include "FX/AudioEffectDistortion.hpp"
%include "FX/AudioEffectEarlyDelay.hpp"
%include "FX/AudioEffectEarlyReflection.hpp"
%include "FX/AudioEffectFlanger.hpp"
%include "FX/AudioEffectHallReverb.hpp"
%include "FX/AudioEffectLFO.hpp"
%include "FX/AudioEffectParametricEQ.hpp"
%include "FX/AudioEffectPhaser.hpp"
%include "FX/AudioEffectPingPong.hpp"
%include "FX/AudioEffectPitchShifter.hpp"
%include "FX/AudioEffectPlateReverb.hpp"
%include "FX/AudioEffectRingMod.hpp"
%include "FX/AudioEffectRoomReverb.hpp"
%include "FX/AudioEffectStereoDelay.hpp"
%include "FX/AudioEffectSVF.hpp"
%include "FX/AudioEffectTremolo.hpp"
%include "FX/AudioEffectVibrato.hpp"

/*
%inline %{
    const int BufferSize = 256;
    Default noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}
*/