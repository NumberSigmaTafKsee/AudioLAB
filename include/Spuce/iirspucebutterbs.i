%module spucebutterbs
%{
typedef float DspFloatType;
#include "SpuceButterBandStop.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceButterBandStop.hpp"