%module scalar
%{
#include "Std/StdScalar.h"
%}

%include "stdint.i"
%include "Std/StdScalar.h"

%inline %{
    typedef float f32;
    typedef double f64;
%}
%template(int_scalar) Std::StdScalar<int>;
%template(uint_scalar) Std::StdScalar<unsigned int>;

%template(long_scalar) Std::StdScalar<long>;
%template(ulong_scalar) Std::StdScalar<unsigned long>;

%template(float_scalar) Std::StdScalar<float>;
%template(double_scalar) Std::StdScalar<double>;

using namespace Std;
%template(absf) abs<f32>;
%template(cubef) cube<f32>;
%template(sqrtf) sqrt<f32>;
%template(expf) exp<f32>;
%template(exp2f) exp2<f32>;
%template(logf) log<f32>;
%template(log10f) log10<f32>;
%template(log2f) log2<f32>;
%template(logbf) logb<f32>;
%template(powf) pow<f32>;
%template(floorf) floor<f32>;
%template(acosf) acos<f32>;
%template(asinf) asin<f32>;
%template(atanf) atan<f32>;
%template(atan2f) atan2<f32>;
%template(cosf) cos<f32>;
%template(sinf) sin<f32>;
%template(tanf) tan<f32>;
%template(coshf) cosh<f32>;
%template(sinhf) sinh<f32>;
%template(tanhf) tanh<f32>;
%template(lgammaf) lgamma<f32>;
%template(acoshf) acosh<f32>;
%template(asinhf) asinh<f32>;
%template(atanhf) atanh<f32>;
%template(cbrtf) cbrt<f32>;
%template(ceilf) cbrt<f32>;
%template(copysignf) copysign<f32>;
%template(erff) erf<f32>;
%template(erfcf) erfc<f32>;
%template(expm1f) expm1<f32>;
%template(fdimf) fdim<f32>;
%template(fmaf) fma<f32>;
%template(fmaxf) fmax<f32>;
%template(fminf) fmin<f32>;
%template(fmodf) fmod<f32>;
%template(fpclassifyf) fpclassify<f32>;
%template(hypotf) hypot<f32>;
%template(ilogbf) ilogb<f32>;
%template(isfinitef) isfinite<f32>;
%template(isgreaterf) isgreater<f32>;
%template(isgreaterequalf) isgreaterequal<f32>;
%template(isinff) isinf<f32>;
%template(islessf) isless<f32>;
%template(islessequalf) islessequal<f32>;
%template(isnanf) isnan<f32>;
%template(isnormalf) isnormal<f32>;
%template(isunorderedf) isunordered<f32>;
%template(ldexpf) ldexp<f32>;
%template(lgammaf) lgamma<f32>;
%template(llrintf) llrint<f32>;
%template(llroundf) llround<f32>;
%template(log1pf) log1p<f32>;
%template(lrintf) lrint<f32>;
%template(lroundf) lround<f32>;
%template(nanf) nan<f32>;
%template(nanff) nanf<f32>;
%template(nanlf) nanl<f32>;
%template(nearbyintf) nearbyint<f32>;
%template(nextafterf) nextafter<f32>;
%template(nexttowardf) nexttoward<f32>;
%template(remainderf) remainder<f32>;
%template(rintf) rint<f32>;
%template(roundf) round<f32>;
%template(scalblnf) scalbln<f32>;
%template(scalbnf) scalbn<f32>;
%template(squaref) square<f32>;
%template(tgammaf) tgamma<f32>;
%template(truncf) trunc<f32>;

%template(absd) abs<f64>;
%template(sqrtd) sqrt<f64>;
%template(expd) exp<f64>;
%template(exp2d) exp2<f64>;
%template(logd) log<f64>;
%template(log10d) log10<f64>;
%template(log2d) log2<f64>;
%template(logbd) logb<f64>;
%template(powd) pow<f64>;
%template(floord) floor<f64>;
%template(acosd) acos<f64>;
%template(asind) asin<f64>;
%template(atand) atan<f64>;
%template(atan2d) atan2<f64>;
%template(cosd) cos<f64>;
%template(sind) sin<f64>;
%template(tand) tan<f64>;
%template(coshd) cosh<f64>;
%template(sinhd) sinh<f64>;
%template(tanhd) tanh<f64>;
%template(lgammad) lgamma<f64>;
%template(acoshd) acosh<f64>;
%template(asinhd) asinh<f64>;
%template(atanhd) atanh<f64>;
%template(cbrtd) cbrt<f64>;
%template(ceild) cbrt<f64>;
%template(copysignd) copysign<f64>;
%template(erfd) erf<f64>;
%template(erfcd) erfc<f64>;
%template(expm1d) expm1<f64>;
%template(fdimd) fdim<f64>;
%template(fmad) fma<f64>;
%template(fmaxd) fmax<f64>;
%template(fmind) fmin<f64>;
%template(fmodd) fmod<f64>;
%template(fpclassifyd) fpclassify<f64>;
%template(hypotd) hypot<f64>;
%template(ilogbd) ilogb<f64>;
%template(isfinited) isfinite<f64>;
%template(isgreaterd) isgreater<f64>;
%template(isgreaterequald) isgreaterequal<f64>;
%template(isinfd) isinf<f64>;
%template(islessd) isless<f64>;
%template(islessequald) islessequal<f64>;
%template(isnand) isnan<f64>;
%template(isnormald) isnormal<f64>;
%template(isunorderedd) isunordered<f64>;
%template(ldexpd) ldexp<f64>;
%template(lgammad) lgamma<f64>;
%template(llrintd) llrint<f64>;
%template(llroundd) llround<f64>;
%template(log1pd) log1p<f64>;
%template(lrintd) lrint<f64>;
%template(lroundd) lround<f64>;
%template(nand) nan<f64>;
%template(nanfd) nanf<f64>;
%template(nanld) nanl<f64>;
%template(nearbyintd) nearbyint<f64>;
%template(nextafterd) nextafter<f64>;
%template(nexttowardd) nexttoward<f64>;
%template(remainderd) remainder<f64>;
%template(rintd) rint<f64>;
%template(roundd) round<f64>;
%template(scalblnd) scalbln<f64>;
%template(scalbnd) scalbn<f64>;
%template(squared) square<f64>;
%template(tgammad) tgamma<f64>;
%template(truncd) trunc<f64>;