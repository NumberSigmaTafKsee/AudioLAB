%module spucebutterbp
%{
typedef float DspFloatType;
#include "SpuceButterBandPass.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceButterBandPass.hpp"