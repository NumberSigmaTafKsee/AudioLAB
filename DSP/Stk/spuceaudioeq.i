%module spucebutterlp
%{
typedef float DspFloatType;
#include "SpuceAudioEqualizer.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceAudioEqualizer.hpp"
