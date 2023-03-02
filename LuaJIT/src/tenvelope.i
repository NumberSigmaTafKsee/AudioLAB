%module tenvelope
%{
#include "Envelope/TEnvelope.hpp"
%}
%include "stdint.i"
%include "std_vector.i"

%include "Envelope/TEnvelope.hpp"

%template(FloatADSR) SoundWave::ADSR<float>;
%template(FloatGammaEnv) SoundWave::GammaEnv<float>;
%template(FloatParameters) SoundWave::parameter_struct<float>;
%template(FloatParametricEnv) SoundWave::ParametricEnvelope<float>;
