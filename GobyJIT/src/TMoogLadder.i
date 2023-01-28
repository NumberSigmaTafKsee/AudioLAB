%module TMoogLadder
%{
#include "TMoogLadder.hpp"
%}

%include "TMoogLadder.hpp"

%template( FloatTMoogLadder ) SoundWave::TMoogLadder<float>;
%template( DoubleTMoogLadder ) SoundWave::TMoogLadder<double>;