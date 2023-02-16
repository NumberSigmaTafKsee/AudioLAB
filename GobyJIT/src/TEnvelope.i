%module TEnvelope
%{
#include "TEnvelope.hpp"
%}

%include "TEnvelope.hpp"

%template (FloatADSR) SoundWave::ADSR<float>;
%template (FloatGammaEnv) SoundWave::GammaEnv<float>;
%template (FloatParametricEnv) SoundWave::ParametricEnvelope<float>;


%template (DoubleADSR) SoundWave::ADSR<double>;
%template (DoubleGammaEnv) SoundWave::GammaEnv<double>;
%template (DoubleParametricEnv) SoundWave::ParametricEnvelope<double>;
