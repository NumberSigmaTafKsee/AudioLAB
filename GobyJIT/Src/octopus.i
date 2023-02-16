%module octopus
%{
#include "Octopus.hpp"
%}
%include "std_math.i"
%include "std_vector.i"
%include "std_complex.i"
%include "std_string.i"

%template(double_vector) std::vector<double>;
%template(cdouble_vector) std::vector<std::complex<double>>;
%template(complex) std::complex<double>;

%ignore oct;
%ignore interpreter;

%include "Octopus.hpp"

