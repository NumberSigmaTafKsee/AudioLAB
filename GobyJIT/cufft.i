%module cufft 
%{

#include "cufft.h"

%}

%constant int cufft_forward = -1;
%constant int cufft_inverse = 1;

%include "std_math.i"
%include "std_vector.i"
%include "cufft.h"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;