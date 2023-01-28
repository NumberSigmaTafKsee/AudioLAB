%module xcorr
%{
extern "C" {
#include "xcorr.h"
}
%}

%include "stdint.i"
%include "std_vector.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "xcorr.h"
