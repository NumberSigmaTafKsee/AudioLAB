%module TADSR
%{
#include "TGammaEnv.hpp"

using namespace SoundAlchemy;
%}

%include "TGammaEnv.hpp"

%template(FloatGamma)  SoundAlchemy::TGammaEnv<float>;
%template(DoubleGamma) SoundAlchemy::TGammaEnv<double>;
