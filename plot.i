%module plot
%{
    #include "Plot.hpp"
    using namespace std;
%}
%include "stdint.i"
%include "std_math.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_exception.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%include "Plot.hpp"

%template(Plot_Float) Plot<float>;
%template(Plot_Double) Plot<double>;