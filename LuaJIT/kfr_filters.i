%module kfr_filters
%{
#include "IIRFilters.hpp"
using namespace Filters;
%}
%include "std_vector.i"
%template(float_vector) std::vector<float>;
%template(double_vecotr) std::vector<double>;

%include "IIRFilters.hpp"
