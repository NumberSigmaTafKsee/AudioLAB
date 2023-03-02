%module TDiodeLadder
%{
#include "TDiodeLadder.hpp"
%}

%include "TDiodeLadder.hpp"

%template(FloatDiodeLadder) SoundWave::TDiodeLadder<float>;
%template(DoubleDiodeLadder) SoundWave::TDiodeLadder<double>;