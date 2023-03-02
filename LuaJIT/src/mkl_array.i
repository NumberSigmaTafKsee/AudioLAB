%module mkl_array
%{
#include "mkl_array.hpp"
%}
%include "mkl_array.hpp"

%extend Sample::Array
{
    T __getitem(size_t i) { return (*$self)[i]; }
    void __setitem(size_t i, const T &val) { (*$self)[i] = val;}
}

%template(float_array) Sample::Array<float>;

%template(sumf) Sample::sum<float>;
%template(minf) Sample::min<float>;
%template(maxf) Sample::max<float>;
%template(min_elementf) Sample::min_element<float>;
%template(max_elementf) Sample::max_element<float>;