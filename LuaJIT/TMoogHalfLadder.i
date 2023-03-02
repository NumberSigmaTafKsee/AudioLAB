%module TMoogHalfLadder
%{
#include "TMoogHalfLadder.hpp"
%}

%include "TMoogHalfLadder.hpp"

%template (FloatTMoogHalfLadder) SoundWave::TMoogHalfLadder<float>;
%template (DoubleMoogHalfLadder) SoundWave::TMoogHalfLadder<double>;
