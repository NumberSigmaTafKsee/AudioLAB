%module TOberheim
%{
#include "TFormantFilter.hpp"
%}

%include "TFormantFilter.hpp"

%template( FloatTFormantFilter ) SoundWave::TFormantFilter<float>;
%template( DoubleTFormantFiler)  SoundWave::TFormantFilter<double>;