%module KfrDSP1
%{
#include "Kfr/KfrDSP1.hpp"
using namespace DSP1;
%}

%include "std_math.i"
%include "kfr.i"

%template(sample_vector_f32) sample_vector<float>;
%template(sample_vector_f64) sample_vector<double>;
%template(sample_vector_f80) sample_vector<long double>;

%template(sample_matrix_f32) sample_matrix<float>;
%template(sample_matrix_f64) sample_matrix<double>;
%template(sample_matrix_f80) sample_matricx<long double>;

%template(complex_f32) complex<float>;
%template(complex_f64) complex<double>;
%template(complex_f80) complex<long double>;

%template(complex_vector_f32) complex_vector<float>;
%template(complex_vector_f64) complex_vector<double>;
%template(complex_vector_f80) complex_vector<long double>;