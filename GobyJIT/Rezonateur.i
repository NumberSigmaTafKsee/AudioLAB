%module Rezonateur
%{
#include "RezonateurFilter.hpp"
#include <vector>
using namespace Filters::Formant::Rezonateur;
%}
%include "std_math.i"
%include "std_string.i"
%include "std_vector.i"

%template (float_vector) std::vector<float>;
%template (double_vector) std::vector<double>;

%include "RezonateurFilter.hpp"
