%module valarray
%{
#include <valarray>
#include <algorithm>
%}

%include "stdint.i"
%include "std_valarray.i"

%template(size_array) std::valarray<size_t>;
%template(int_array) std::valarray<int>;
%template(uint_array) std::valarray<unsigned int>;
%template(long_array) std::valarray<long>;
%template(ulong_array) std::valarray<unsigned long>;
%template(float_array) std::valarray<float>;
%template(double_array) std::valarray<double>;

%template(int_gslice_array) std::gslice_array<int>;
%template(uint_gslice_array) std::gslice_array<unsigned int>;
%template(long_gslice_array) std::gslice_array<long>;
%template(ulong_gslice_array) std::gslice_array<unsigned long>;
%template(float_gslice_array) std::gslice_array<float>;
%template(double_gslice_array) std::gslice_array<double>;

/* not sure
%template(int_indirect_array) std::indirect_array<int>;
%template(uint_indirect_array) std::indirect_array<unsigned int>;
%template(long_indirect_array) std::indirect_array<long>;
%template(ulong_indirect_array) std::indirect_array<unsigned long>;
%template(float_indirect_array) std::indirect_array<float>;
%template(double_indirect_array) std::indirect_array<double>;
*/

%template(int_mask_array) std::mask_array<int>;
%template(uint_mask_array) std::mask_array<unsigned int>;
%template(long_mask_array) std::mask_array<long>;
%template(ulong_mask_array) std::mask_array<unsigned long>;
%template(float_mask_array) std::mask_array<float>;
%template(double_mask_array) std::mask_array<double>;

%template(int_slice_array) std::slice_array<int>;
%template(uint_slice_array) std::slice_array<unsigned int>;
%template(long_slice_array) std::slice_array<long>;
%template(ulong_slice_array) std::slice_array<unsigned long>;
%template(float_slice_array) std::slice_array<float>;
%template(double_slice_array) std::slice_array<double>;

// dont know why
//%template(sorti) std::sort<int>;
//%template(swapi) std::swap<int>;
//%template(copyi) std::copy<int>;
//%template(random_shufflei) std::random_shuffle<int>;

%template(absi) std::abs<int>;
%template(powi) std::pow<int>;

%template(absl) std::abs<long>;
%template(powl) std::pow<long>;

%template(fabs) std::abs<float>;
%template(expf) std::exp<float>;
%template(logf) std::log<float>;
%template(log10f) std::log10<float>;
%template(sqrtf) std::sqrt<float>;
%template(sinf) std::sin<float>;
%template(cosf) std::cos<float>;
%template(tanf) std::tan<float>;
%template(asinf) std::asin<float>;
%template(acosf) std::acos<float>;
%template(atanf) std::atan<float>;
%template(sinhf) std::sinh<float>;
%template(coshf) std::cosh<float>;
%template(tanhf) std::tanh<float>;
//%template(asinhf) std::asinh<float>;
//%template(acoshf) std::acosh<float>;
//%template(atanhf) std::atanh<float>;
%template(powf) std::pow<float>;
%template(atan2f) std::atan<float>;

%template(abs) std::abs<double>;
%template(exp) std::exp<double>;
%template(log) std::log<double>;
%template(log10) std::log10<double>;
%template(sqrt) std::sqrt<double>;
%template(sin) std::sin<double>;
%template(cos) std::cos<double>;
%template(tan) std::tan<double>;
%template(asin) std::asin<double>;
%template(acos) std::acos<double>;
%template(atan) std::atan<double>;
%template(sinh) std::sinh<double>;
%template(cosh) std::cosh<double>;
%template(tanh) std::tanh<double>;
//%template(asinhf) std::asinh<double>;
//%template(acoshf) std::acosh<double>;
//%template(atanhf) std::atanh<double>;
%template(pow) std::pow<double>;
%template(atan2) std::atan<double>;