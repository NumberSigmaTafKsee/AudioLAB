%module spuceremezfir
%{
typedef float DspFloatType;
#include "SpuceRemezFIR.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceRemezFIR.hpp"