%module simpleresampler
%{
typedef float DspFloatType;
#include "DSP/SimpleResampler.hpp"
%}

%include "std_vector.i"

typedef float DspFloatType;

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
'
%include "DSP/SimpleResampler.hpp"

