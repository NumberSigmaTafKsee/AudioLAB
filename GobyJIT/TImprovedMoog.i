%module TImprovedMoog
%{
#include "TImprovedMoog.hpp"
%}

%include "TImprovedMoog.hpp"

%template(FloatImprovedMoog) SoundWave::TImprovedMoog<float>;
%template(DoubleImprovedMoog) SoundWave::TImprovedMoog<float>;
