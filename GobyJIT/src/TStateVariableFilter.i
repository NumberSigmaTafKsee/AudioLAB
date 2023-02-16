%module TStateVariableFilter
%{
#include "TStateVariableFilter.hpp"
%}

%include "TStateVariableFilter.hpp"

%template( FloatStateVariableFilter ) SoundWave::TStateVariableFilter<float>;
%template( DoubleStateVariableFilter ) SoundWave::TStateVariableFilter<double>;