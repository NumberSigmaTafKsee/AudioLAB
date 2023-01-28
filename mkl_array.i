%module mkl_array
%{
#include "mkl_array.hpp"
%}
%include "mkl_array.hpp"

%extend Sample::Array
{
    T       __getitem__(size_t i) { return (*$self)[i]; }
    void    __setitem__(size_t i, const T &val) { (*$self)[i] = val;}
}

%template(float_array) Sample::Array<float>;
%template(double_array) Sample::Array<double>;

%template(sumf) Sample::sum<float>;
%template(minf) Sample::min<float>;
%template(maxf) Sample::max<float>;
%template(min_elementf) Sample::min_element<float>;
%template(max_elementf) Sample::max_element<float>;

%template(sumd) Sample::sum<double>;
%template(mind) Sample::min<double>;
%template(maxd) Sample::max<double>;
%template(min_elementd) Sample::min_element<double>;
%template(max_elementd) Sample::max_element<double>;
