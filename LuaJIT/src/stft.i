%module stft
%{
#include "STFT.h"
%}
%include "std_vector.i"
%include "STFT.h"

%template(float_vector) std::vector<float>;