%module BiQuad
%{
#include "BiQuad.hpp"
%}

%include "std_complex.i"
%include "std_vector.i"

%include "BiQuad.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(fcomplex) std::complex<float>;
%template(zcomplex) std::complex<double>;
