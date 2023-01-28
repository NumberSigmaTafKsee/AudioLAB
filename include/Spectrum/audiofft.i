%module audiofft
%{
#include "AudioFFT.h"
using namespace audiofft;
%}

%include "stdint.i"
%include "std_vector.i"

%template(float_vector) std::vector<float>;

%include "AudioFFT.h"