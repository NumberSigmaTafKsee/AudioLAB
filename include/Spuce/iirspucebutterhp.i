%module spucebutterhp
%{
typedef float DspFloatType;
#include "SpuceButterHighPass.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceButterHighPass.hpp"