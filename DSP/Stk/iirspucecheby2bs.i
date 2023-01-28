%module spucecheby2bs
%{
typedef float DspFloatType;
#include "SpuceChebyshev2Bandstop.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceChebyshev2Bandstop.hpp"
