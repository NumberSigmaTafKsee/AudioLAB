%module convolution_fftw
%{
#include "convolution_fftw.h"
%}
%include "std_vector.i"
%template(double_vector) std::vector<double>;
%include "convolution_fftw.h"
