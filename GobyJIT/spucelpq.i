%module spucelpq
%{
typedef float DspFloatType;
#include "SpuceLPQ.hpp"
%}
typedef float DspFloatType;
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "SpuceLPQ.hpp"
