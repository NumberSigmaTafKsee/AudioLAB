%module spucecheby2bp
%{
typedef float DspFloatType;
#include "SpuceChebyshev2Bandpass.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceChebyshev2Bandpass.hpp"
