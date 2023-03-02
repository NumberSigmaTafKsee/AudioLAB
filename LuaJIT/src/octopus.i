%module octopus
%{
#include "Octopus.hpp"
%}
%include "std_vector.i"
%include "std_complex.i"

%template(double_vector) std::vector<double>;
%template(cdouble_vector) std::vector<std::complex<double>>;
%template(complex) std::complex<double>;

%ignore oct;
%ignore interpreter;

%include "Octopus.hpp"

%template (octopus_matrix) OctopusMatrix<double>;
%template (octopus_cmatrix) OctopusMatrix<std::complex<double>>;