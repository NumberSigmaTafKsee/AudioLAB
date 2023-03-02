%module vector
%{
#include <vector>
%}
%include "std_vector.i"

%template(float_vector) std::vector<float>;