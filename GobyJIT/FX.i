%module FX
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



#include "FX/AllPass.hpp"
#include "FX/AllpassFilter.hpp"
#include "FX/AllPassFilters.hpp"
#include "FX/CSmoothFilters.hpp"
#include "FX/DCBlock.hpp"
#include "FX/DCFilter.hpp"
#include "FX/effects_biquad.h"
#include "FX/cppfilters.hpp"

#include "FX/ExpSmootherCascade.hpp"

#include "FX/HardLimiter.hpp"
#include "FX/Limiter.hpp"
#include "FX/LimiterDsp.hpp"
#include "FX/PeakHoldCascade.hpp"
#include "FX/OutputLimiter.hpp"
#include "FX/PeakLimiter.hpp"

#include "FX/IIRDCBlock.hpp"
#include "FX/IIRDCFilter.hpp"

#include "FX/LFO.hpp"
#include "FX/LFO9000.hpp"
#include "FX/LowFrequencyOscillator.hpp"
%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"



%include "FX/Amplifier.hpp"
%include "FX/AmplifierFold.hpp"
%include "FX/Amplifiers.hpp"
%include "FX/AmplifiersUdo1.hpp"
%include "FX/Diode.hpp"
%include "FX/DiodeClipper.hpp"
%include "FX/DiodeSim.hpp"
%include "FX/DiodeSimulator.hpp"
%include "FX/DistortionFunctions.hpp"
%include "FX/DistortionProcessors.hpp"
%include "FX/MDFM-1000.hpp"
%include "FX/ClipFunctions.hpp"
%include "FX/ClipperCircuit.hpp"
%include "FX/ClippingFunctions.hpp"
%include "FX/ClipSerpentCurve.hpp"

%include "FX/SstWaveshaper.hpp"
%include "FX/SstFilters.hpp"


/*
%include "FX/ringbuffer.hpp"
%include "FX/RingBuffer.hpp"
%include "FX/RingBufferProcessor.hpp"
%include "FX/RingBuffers.hpp"


%include "FX/BasicDelay.hpp"
%include "FX/BasicDelayLine.hpp"
%include "FX/BasicDelayLineStereo.hpp"
%include "FX/CircularDelay.hpp"
%include "FX/Comb.hpp"
%include "FX/CombFilter.hpp"
%include "FX/CombFilters.hpp"
%include "FX/CrossDelayLine.hpp"
%include "FX/DelayFilters.hpp"
%include "FX/DelayLine.hpp"
%include "FX/DelayLines.hpp"
%include "FX/DelayProcessors.hpp"
%include "FX/Delays.hpp"
%include "FX/DelaySmooth.hpp"
%include "FX/DelaySyncedTapDelayLine.hpp"
%include "FX/Mu45.hpp"
%include "FX/Mu45FilterCalc.hpp"
%include "FX/ERTapDelayLine.hpp"
%include "FX/Moorer.hpp"
%include "FX/MoorerStereo.hpp"
%include "FX/PPDelayLine.hpp"
%include "FX/reverb.hpp"
%include "FX/StereoDelay.hpp"
%include "FX/Schroeder.hpp"
%include "FX/SchroederAllpass.hpp"
%include "FX/SchroederImproved.hpp"
%include "FX/SchroederStereo.hpp"
*/

/*
//%include "FX/DynamicProcessors.hpp"
%include "FX/Compreezor.hpp"
%include "FX/Compressor.hpp"
//%include "FX/CTAGDRCCompressor.hpp"
%include "FX/DistortionCompressor.hpp"
%include "FX/Limiter.hpp"
%include "FX/LimiterDsp.hpp"
%include "FX/OutputLimiter.hpp"
%include "FX/PeakHoldCascade.hpp"
%include "FX/PeakLimiter.hpp"
%include "FX/PKLimiter.hpp"
%include "FX/HardLimiter.hpp"
%include "FX/SimpleDynamics.hpp"


%include "FX/AllPass.hpp"
%include "FX/AllpassFilter.hpp"
%include "FX/AllPassFilters.hpp"
%include "FX/onepole.hpp"
%include "FX/OnePole.hpp"
%include "FX/OnePoleAllPassFilter.hpp"
%include "FX/OnePoleEnvelope.hpp"
%include "FX/OnePoleEnvelopes.hpp"
%include "FX/OrfanidisEQ.hpp"
%include "FX/paraeq.hpp"
%include "FX/cppfilters.hpp"
%include "FX/CSmoothFilters.hpp"
%include "FX/DCBlock.hpp"
%include "FX/DCFilter.hpp"
%include "FX/effects_biquad.h"
%include "FX/EMAFilter.hpp"
%include "FX/FilterBank.h"
%include "FX/Filters.h"
%include "FX/Filters.hpp"
%include "FX/BodeShifter.hpp"
%include "FX/IIRDCBlock.hpp"
%include "FX/IIRDCFilter.hpp"
%include "FX/LiquidMoog.hpp"
%include "FX/LiquidNeuron.hpp"
%include "FX/MoogFilters.hpp"
%include "FX/KocMocPhasor.hpp"
%include "FX/Mu45.hpp"
%include "FX/Mu45FilterCalc.hpp"
%include "FX/PolyPhaseFilter.hpp"
%include "FX/PolyPhaseFilterBank.hpp"
%include "FX/OctaveFilterBank.hpp"
%include "FX/RBJ.hpp"
%include "FX/RBJFilter.hpp"
%include "FX/RCFilter.hpp"
%include "FX/MovingAverageFilter.hpp"
%include "FX/SVF.hpp"
%include "FX/TwoPoleEnvelopes.hpp"
%include "FX/StateVariableFilters.hpp"
%include "FX/ThirdOrderAllPole.hpp"
%include "FX/ToneStack.hpp"
%include "FX/RKLadderFilter.hpp"
%include "FX/RKSimulationModel.h"
*/

//%include "FX/DCAProcessors.hpp"
//%include "FX/DCFProcessors.hpp"
//%include "FX/DCOProcessors.hpp"

/*
%include "FX/FXCE2Chorus.hpp"
%include "FX/FXChorus.hpp"
%include "FX/FXChorus2.hpp"
%include "FX/FXChorusDelay.hpp"
%include "FX/FXFlanger.hpp"
%include "FX/FXFreeVerb3.hpp"
%include "FX/FXLeslie.hpp"
%include "FX/FXLimiterDsp.hpp"
%include "FX/FXOutputLimiter.hpp"
%include "FX/FXProcessor.hpp"
%include "FX/FXProcessors.hpp"
%include "FX/FXPsola.hpp"
%include "FX/FXRoboVerb.hpp"
%include "FX/FXYKChorus.hpp"
*/

/*
%include "FX/RandomSmoother.hpp"
%include "FX/ResamplerProcessors.hpp"
%include "FX/RMCircularBuffer.hpp"
%include "FX/rt-wdf.hpp"
%include "FX/Rubberband.hpp"
%include "FX/ScaledMapParam.hpp"
%include "FX/SlewLimiter.hpp"
%include "FX/Stereofiyer.hpp"
%include "FX/StereoVector.hpp"
%include "FX/TinyTricksOscillators.hpp"
%include "FX/TriggerFishNoise.hpp"
%include "FX/TriggerFishVanDerPol.hpp"
%include "FX/TSClipper.hpp"
%include "FX/TSTone.hpp"
*/

//%include "FX/AnalogCircuits.hpp"
//%include "FX/AnalogSVF.hpp"
//%include "FX/AudioTK.hpp"

//%include "FX/HammerFX.cpp
//%include "FX/RackFXRecChord.cpp
//%include "FX/spline.h
//%include "FX/uri_table.h
//%include "FX/ladspa-util.h
//%include "FX/LagrangeInterpolation.hpp"
//%include "FX/LV2Host.cpp
//%include "FX/Matsuko5Pendulum.hpp"
//%include "FX/Modulation.hpp"
//%include "FX/PartialSynth.hpp"
//%include "FX/racka3.cpp
//%include "FX/racka3.hpp"
//%include "FX/AudioProcessor.hpp"
//%include "FX/ReissFX.hpp"
//%include "FX/AudioEffects.hpp"
/*
%include "FX/BandLimitedOscillators.hpp"
%include "FX/Bezier.hpp"
%include "FX/BezierDistortion.hpp"
%include "FX/Blit.hpp"
%include "FX/BlitSaw.hpp"
%include "FX/BlitSquare.hpp"
%include "FX/CrystalAliasFreeOscillator.hpp"
%include "FX/DPW.hpp"
%include "FX/DPWOscillators.hpp"
%include "FX/DSF.hpp"
%include "FX/minBLEP.hpp"
%include "FX/MinBLEP.hpp"
%include "FX/MinBLEPOsc.hpp"
%include "FX/NOBSNonlinearOscillator.hpp"
%include "FX/NOHMADStrangeAttractors.hpp"
%include "FX/LFO9000.hpp"
%include "FX/Phasor.hpp"
%include "FX/PolyBLEP.hpp"
%include "FX/PolygonalOscillator.hpp"
%include "FX/VanDerPolOscillator.hpp"
%include "FX/Walsh.hpp"
%include "FX/WaveFile.hpp"
%include "FX/WaveFourierWave.hpp"
%include "FX/WaveGenerators.hpp"
%include "FX/WaveguideLibrary.hpp"
%include "FX/WaveShaperATanSoftClip.hpp"
%include "FX/Waveshapers.hpp"
%include "FX/Waveshaping.hpp"
//%include "FX/pinknoise.hpp"
//%include "FX/PinkNoise.hpp"
//%include "FX/AdaptiveFilterProcessors.hpp"
//%include "FX/DigitalFilters.hpp"
//%include "FX/ExpSmootherCascade.hpp"
//%include "FX/DAFXProcessors.hpp"
//%include "FX/DataStructures.hpp"
//%include "FX/DSFWALSH.cpp
//%include "FX/effects_biquad.cpp
//%include "FX/FilterBank.cpp
//%include "FX/Fir1.cpp
//%include "FX/Fir1.h"
//%include "FX/fir_filter.cpp
//%include "FX/fir_filter.h"
*/

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
