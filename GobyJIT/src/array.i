%module array
%{
#include "Eigen"
#include <iostream>
using namespace Eigen;
%}

%include "std_complex.i"
%include "std_vector.i"
%include "eigen-array.i"

%template(fcomplex) std::complex<float>;
%template(zcomplex) std::complex<double>;

%template(ArrayXf) Eigen::Array<float,Eigen::Dynamic,1>;
%template(ArrayXd) Eigen::Array<double,Eigen::Dynamic,1>;
%template(ArrayXi) Eigen::Array<int,Eigen::Dynamic,1>;
%template(ArrayXl) Eigen::Array<long,Eigen::Dynamic,1>;
%template(ArrayXc) Eigen::Array<std::complex<float>,Eigen::Dynamic,1>;
%template(ArrayXz) Eigen::Array<std::complex<double>,Eigen::Dynamic,1>;

%template(absf) Ops::abs<float>;
%template(abs2f) Ops::abs2<float>;
%template(inversef) Ops::inverse<float>;
%template(expf) Ops::exp<float>;
%template(logf) Ops::log<float>;
%template(log1pf) Ops::log1p<float>;
%template(log10f) Ops::log10<float>;
%template(powf) Ops::pow<float>;
%template(sqrtf) Ops::sqrt<float>;
%template(rsqrtf) Ops::rsqrt<float>;
%template(square) Ops::square<float>;
%template(cube) Ops::cube<float>;
%template(sinf) Ops::sin<float>;
%template(cosf) Ops::cos<float>;
%template(tanf) Ops::tan<float>;
%template(asinf) Ops::asin<float>;
%template(acosf) Ops::acos<float>;
%template(atanf) Ops::atan<float>;
%template(sinhf) Ops::sinh<float>;
%template(coshf) Ops::cosh<float>;
%template(tanhf) Ops::tanh<float>;
%template(asinhf) Ops::asinh<float>;
%template(acoshf) Ops::acosh<float>;
%template(atanhf) Ops::atanh<float>;
%template(floorf) Ops::floor<float>;
%template(ceilf) Ops::ceil<float>;
%template(roundf) Ops::round<float>;
%template(rintf) Ops::rint<float>;
%template(sizef) Ops::size<float>;
%template(randomf) Ops::random<float>;
%template(fillf) Ops::fill<float>;
%template(colsf) Ops::cols<float>;
%template(resizef) Ops::resize<float>;


