%module adsr
%{
#include "ADSR.h"
#include <vector>
%}
%include "std_vector.i"
%include "ADSR.h"

%template (float_vector) std::vector<float>;
