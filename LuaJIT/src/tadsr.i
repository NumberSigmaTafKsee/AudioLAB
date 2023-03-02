%module TADSR
%{
#include "TADSR.hpp"

using namespace SoundAlchemy;
%}

%include "TADSR.hpp"

%template(FloatADSR)  SoundAlchemy::TADSR<float>;
%template(DoubleADSR) SoundAlchemy::TADSR<double>;
