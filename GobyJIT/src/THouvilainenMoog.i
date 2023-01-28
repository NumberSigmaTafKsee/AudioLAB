%module THouvilainenMoog
%{
#include "THouvilainenMoog.hpp"
%}

%include "THouvilainenMoog.hpp"

%template(FloatHouvilainenMoog) SoundWave::THouvilainenMoog<float>;
%template(DoubleHouvilainenMoog) SoundWave::THouvilainenMoog<float>;
