%module dsplpq
%{
#include "DspButterworthLowPassQ.hpp"
%}
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%include "IIRDspFilters.hpp"
%include "DspButterworthLowPassQ.hpp"
