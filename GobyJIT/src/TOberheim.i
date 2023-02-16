%module TOberheim
%{
#include "TOberheim.hpp"
%}

%include "TOberheim.hpp"

%template( FloatTOberheim ) SoundWave::TOberheim<float>;
%template( DoubleTOberheim ) SoundWave::TOberheim<double>;