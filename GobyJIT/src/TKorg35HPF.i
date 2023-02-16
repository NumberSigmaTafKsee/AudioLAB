%module TKorg35HPF
%{
#include "TKorg35HPF.hpp"
%}

%include "TKorg35HPF.hpp"

%template (FloatTKorg35HPF) SoundWave::TKorg35HPF<float>;
%template (DoubleTKorg35HPF) SoundWave::TKorg35HPF<double>;