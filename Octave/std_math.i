%module stdmath
%{
#include <complex>
#include <valarray>
#include <fenv.h>
%}

%include "std_complex.i"

//%template(float_complex) std::complex<float>;
//%template(double_complex) std::complex<double>;

//%include "std_valarray.i"
//%template(float_valarray) std::valarray<float>;
//%template(double_valarray) std::valarray<double>;

//%include "std_eigen.i"

%include "std_limits.i"

%template(float_numeric_limits) std::numeric_limits<float>;
%template(double_numeric_limits) std::numeric_limits<double>;
%template(int_numeric_limits) std::numeric_limits<int>;
%template(uint_numeric_limits) std::numeric_limits<unsigned int>;
%template(short_numeric_limits) std::numeric_limits<short>;
%template(ushort_numeric_limits) std::numeric_limits<unsigned short>;
%template(long_numeric_limits) std::numeric_limits<long>;
%template(ulong_numeric_limits) std::numeric_limits<unsigned long>;
%template(llong_numeric_limits) std::numeric_limits<long long>;
%template(ullong_numeric_limits) std::numeric_limits<unsigned long long>;


%inline %{
  /*
    typedef float f32;
    typedef double f64;
    typedef signed char i8;
    typedef unsigned char u8;
    typedef signed short i16;
    typedef unsigned short u16;
    typedef signed int i32;
    typedef unsigned int u32;
    typedef signed long ilong;
    typedef unsigned long ulong;
    typedef signed long long i64;
    typedef unsigned long long u64;
  */
    
    namespace Ops
    {
    
      template<typename T> T abs(const T& x) { return std::abs(x); }
      template<typename T> T fabs(const T& x) { return std::fabs(x); }
      template<typename T> T acos(const T& x) { return std::acos(x); }
      template<typename T> T acosh(const T& x) { return std::acosh(x); }
      template<typename T> T asinh(const T& x) { return std::asinh(x); }
      template<typename T> T asin(const T& x) { return std::asinh(x); }
      template<typename T> T atan(const T& x) { return std::atan(x); }
      template<typename T> T atan2(const T& x,const T& y) { return std::atan2(x,y); }
      template<typename T> T atanh(const T& x) { return std::atanh(x); }
      template<typename T> T cbrt(const T& x) { return std::cbrt(x); }
      template<typename T> T ceil(const T& x) { return std::ceil(x); }    
      template<typename T> T copysign(const T& x, const T& y) { return std::copysign(x,y); }
      template<typename T> T cos(const T& x) { return std::cos(x); }
      template<typename T> T cosh(const T& x) { return std::cosh(x); }
      template<typename T> T erf(const T& x) { return std::erf(x); }
      template<typename T> T erfc(const T& x) { return std::erfc(x); }
      template<typename T> T exp(const T& x) { return std::exp(x); }
      template<typename T> T exp2(const T& x) { return std::exp2(x); }
      template<typename T> T expm1(const T& x) { return std::expm1(x); }
      template<typename T> T fdim(const T & x, const T & y) { return std::fdim(x,y); }
      template<typename T> T floor(const T & x) { return std::floor(x); }
      template<typename T> T fma(const T & x, const T & y, const T& z) { return std::fma(x,y,z); }
      template<typename T> T fmax(const T & x, const T & y) { return std::fmax(x,y); }
      template<typename T> T fmin(const T & x, const T & y) { return std::fmax(x,y); }
      template<typename T> T fmod(const T & x, const T & y) { return std::fmod(x,y); }
      template<typename T> int fpclassify(const T & x) { return std::fpclassify(x); }
      template<typename T> T hypot(const T & x, const T & y) { return std::hypot(x,y); }
      template<typename T> int ilogb(const T & x) { return std::ilogb(x); }
      template<typename T> bool isfinite(const T & x) { return std::isfinite(x); }
      template<typename T> bool isgreater(const T & x, const T & y) { return std::isgreater(x,y); }
      template<typename T> bool isgreaterequal(const T & x, const T & y) { return std::isgreaterequal(x,y); }
      template<typename T> bool isinf(const T & x) { return std::isinf(x); }
      template<typename T> bool isless(const T & x, const T & y) { return std::isless(x,y); }
      template<typename T> bool islessequal(const T & x, const T & y) { return std::islessequal(x,y); }
      template<typename T> bool islessgreater(const T & x, const T & y) { return std::islessgreater(x,y); }
      template<typename T> bool isnan(const T & x) { return std::isnan(x); }
      template<typename T> bool isnormal(const T & x) { return std::isnormal(x); }
      template<typename T> bool isunordered(const T & x, const T& y) { return std::isunordered(x,y); }
      template<typename T> T ldexp(const T & x, int exp) { return std::ldexp(x,exp); }
      template<typename T> T lgamma(const T & x) { return std::lgamma(x); }
      template<typename T> T llrint(const T & x) { return std::llrint(x); }
      template<typename T> T llround(const T & x) { return std::llround(x); }
      template<typename T> T log(const T & x) { return std::log(x); }
      template<typename T> T log10(const T & x) { return std::log10(x); }
      template<typename T> T log1p(const T & x) { return std::log1p(x); }
      template<typename T> T log2(const T & x) { return std::log2(x); }
      template<typename T> T logb(const T & x) { return std::logb(x); }
      template<typename T> T lrint(const T & x) { return std::lrint(x); }
      template<typename T> T lround(const T & x) { return std::lround(x); }
      template<typename T> T nan(const char *tagp) { return std::nan(tagp);}
      template<typename T> T nanf(const char *tagp) { return std::nanf(tagp);}
      template<typename T> T nanl(const char *tagp) { return std::nanl(tagp);}
      template<typename T> T nearbyint(const T &x) { return std::nearbyint(x); }
      template<typename T> T nextafter(const T & x, const T & y) { return std::nextafter(x,y); }
      template<typename T> T nexttoward(const T & x, const T & y) { return std::nexttoward(x,y); }
      template<typename T> T pow(const T & b, const T & e) { return std::pow(b,e); }
      template<typename T> T remainder(const T & n, const T & d) { return std::remainder(n,d); }
      template<typename T> T rint(const T& x) { return std::rint(x); }
      template<typename T> T round(const T& x) { return std::round(x); }
      template<typename T> T scalbln(const T& x, long int n) { return std::scalbln(x,n);}
      template<typename T> T scalbn(const T& x, int n) { return std::scalbln(x,n);}
      template<typename T> bool signbit(const T & x) { return signbit(x); }
      template<typename T> T sin(const T& x) { return std::sin(x); }
      template<typename T> T sinh(const T& x) { return std::sinh(x); }    
      template<typename T> T sqrt(const T& x) { return std::sqrt(x); }
      template<typename T> T square(const T& x) { return x*x; }
      template<typename T> T cube(const T& x) { return x*x*x; }
      template<typename T> T tan(const T& x) { return std::tan(x); }
      template<typename T> T tanh(const T& x) { return std::tanh(x); }        
      template<typename T> T tgamma(const T& x) { return std::tgamma(x); }    
      template<typename T> T trunc(const T& x) { return std::trunc(x); }
      double Huge() { return HUGE_VAL; }
      float Hugef() { return HUGE_VALF; }
      double Infinity() { return INFINITY; }
      double NaN() { return NAN; }
    }    
%}

%template(absf)  Ops::abs<float>;
%template(cubef) Ops::cube<float>;
%template(sqrtf) Ops::sqrt<float>;
%template(expf)  Ops::exp<float>;
%template(exp2f) Ops::exp2<float>;
%template(logf)  Ops::log<float>;
%template(log10f) Ops::log10<float>;
%template(log2f) Ops::log2<float>;
%template(logbf) Ops::logb<float>;
%template(powf) Ops::pow<float>;
%template(floorf) Ops::floor<float>;
%template(acosf) Ops::acos<float>;
%template(asinf) Ops::asin<float>;
%template(atanf) Ops::atan<float>;
%template(atan2f) Ops::atan2<float>;
%template(cosf) Ops::cos<float>;
%template(sinf) Ops::sin<float>;
%template(tanf) Ops::tan<float>;
%template(coshf) Ops::cosh<float>;
%template(sinhf) Ops::sinh<float>;
%template(tanhf) Ops::tanh<float>;
%template(lgammaf) Ops::lgamma<float>;
%template(acoshf) Ops::acosh<float>;
%template(asinhf) Ops::asinh<float>;
%template(atanhf) Ops::atanh<float>;
%template(cbrtf) Ops::cbrt<float>;
%template(ceilf) Ops::cbrt<float>;
%template(copysignf) Ops::copysign<float>;
%template(erff) Ops::erf<float>;
%template(erfcf) Ops::erfc<float>;
%template(expm1f) Ops::expm1<float>;
%template(fdimf) Ops::fdim<float>;
%template(fmaf) Ops::fma<float>;
%template(fmaxf) Ops::fmax<float>;
%template(fminf) Ops::fmin<float>;
%template(fmodf) Ops::fmod<float>;
%template(fpclassifyf) Ops::fpclassify<float>;
%template(hypotf) Ops::hypot<float>;
%template(ilogbf) Ops::ilogb<float>;
%template(isfinitef) Ops::isfinite<float>;
%template(isgreaterf) Ops::isgreater<float>;
%template(isgreaterequalf) Ops::isgreaterequal<float>;
%template(isinff) Ops::isinf<float>;
%template(islessf) Ops::isless<float>;
%template(islessequalf) Ops::islessequal<float>;
%template(isnanf) Ops::isnan<float>;
%template(isnormalf) Ops::isnormal<float>;
%template(isunorderedf) Ops::isunordered<float>;
%template(ldexpf) Ops::ldexp<float>;
%template(lgammaf) Ops::lgamma<float>;
%template(llrintf) Ops::llrint<float>;
%template(llroundf) Ops::llround<float>;
%template(log1pf) Ops::log1p<float>;
%template(lrintf) Ops::lrint<float>;
%template(lroundf) Ops::lround<float>;
%template(nanf) Ops::nan<float>;
%template(nanff) Ops::nanf<float>;
%template(nanlf) Ops::nanl<float>;
%template(nearbyintf) Ops::nearbyint<float>;
%template(nextafterf) Ops::nextafter<float>;
%template(nexttowardf) Ops::nexttoward<float>;
%template(remainderf) Ops::remainder<float>;
%template(rintf) Ops::rint<float>;
%template(roundf) Ops::round<float>;
%template(scalblnf) Ops::scalbln<float>;
%template(scalbnf) Ops::scalbn<float>;
%template(squaref) Ops::square<float>;
%template(tgammaf) Ops::tgamma<float>;
%template(truncf) Ops::trunc<float>;

%template(absd) Ops::abs<double>;
%template(sqrtd) Ops::sqrt<double>;
%template(expd) Ops::exp<double>;
%template(exp2d) Ops::exp2<double>;
%template(logd) Ops::log<double>;
%template(log10d) Ops::log10<double>;
%template(log2d) Ops::log2<double>;
%template(logbd) Ops::logb<double>;
%template(powd) Ops::pow<double>;
%template(floord) Ops::floor<double>;
%template(acosd) Ops::acos<double>;
%template(asind) Ops::asin<double>;
%template(atand) Ops::atan<double>;
%template(atan2d) Ops::atan2<double>;
%template(cosd) Ops::cos<double>;
%template(sind) Ops::sin<double>;
%template(tand) Ops::tan<double>;
%template(coshd) Ops::cosh<double>;
%template(sinhd) Ops::sinh<double>;
%template(tanhd) Ops::tanh<double>;
%template(lgammad) Ops::lgamma<double>;
%template(acoshd) Ops::acosh<double>;
%template(asinhd) Ops::asinh<double>;
%template(atanhd) Ops::atanh<double>;
%template(cbrtd) Ops::cbrt<double>;
%template(ceild) Ops::cbrt<double>;
%template(copysignd) Ops::copysign<double>;
%template(erfd) Ops::erf<double>;
%template(erfcd) Ops::erfc<double>;
%template(expm1d) Ops::expm1<double>;
%template(fdimd) Ops::fdim<double>;
%template(fmad)  Ops::fma<double>;
%template(fmaxd) Ops::fmax<double>;
%template(fmind) Ops::fmin<double>;
%template(fmodd) Ops::fmod<double>;
%template(fpclassifyd) Ops::fpclassify<double>;
%template(hypotd) Ops::hypot<double>;
%template(ilogbd) Ops::ilogb<double>;
%template(isfinited) Ops::isfinite<double>;
%template(isgreaterd) Ops::isgreater<double>;
%template(isgreaterequald) Ops::isgreaterequal<double>;
%template(isinfd) Ops::isinf<double>;
%template(islessd) Ops::isless<double>;
%template(islessequald) Ops::islessequal<double>;
%template(isnand) Ops::isnan<double>;
%template(isnormald) Ops::isnormal<double>;
%template(isunorderedd) Ops::isunordered<double>;
%template(ldexpd) Ops::ldexp<double>;
%template(lgammad) Ops::lgamma<double>;
%template(llrintd) Ops::llrint<double>;
%template(llroundd) Ops::llround<double>;
%template(log1pd) Ops::log1p<double>;
%template(lrintd) Ops::lrint<double>;
%template(lroundd) Ops::lround<double>;
%template(nand) Ops::nan<double>;
%template(nanfd) Ops::nanf<double>;
%template(nanld) Ops::nanl<double>;
%template(nearbyintd) Ops::nearbyint<double>;
%template(nextafterd) Ops::nextafter<double>;
%template(nexttowardd) Ops::nexttoward<double>;
%template(remainderd) Ops::remainder<double>;
%template(rintd) Ops::rint<double>;
%template(roundd) Ops::round<double>;
%template(scalblnd) Ops::scalbln<double>;
%template(scalbnd) Ops::scalbn<double>;
%template(squared) Ops::square<double>;
%template(tgammad) Ops::tgamma<double>;
%template(truncd) Ops::trunc<double>;


// not working in Octave
/*
%constant int fe_divbyzero = FE_DIVBYZERO;
%constant int fe_inexact = FE_INEXACT;
%constant int fe_invalid = FE_INVALID;
%constant int fe_overflow = FE_OVERFLOW;
%constant int fe_underflow = FE_UNDERFLOW;
%constant int fe_all_except = FE_ALL_EXCEPT;
%constant int fe_downward = FE_DOWNWARD;
%constant int fe_tonearest = FE_TONEAREST;
%constant int fe_towardzero = FE_TOWARDZERO;
%constant int fe_upward = FE_UPWARD;
%constant int fe_dfl_env = FE_DFL_ENV;
*/

/*
%inline %{
typedef struct
  {
    unsigned short int __control_word;
    unsigned short int __glibc_reserved1;
    unsigned short int __status_word;
    unsigned short int __glibc_reserved2;
    unsigned short int __tags;
    unsigned short int __glibc_reserved3;
    unsigned int __eip;
    unsigned short int __cs_selector;
    unsigned int __opcode:11;
    unsigned int __glibc_reserved4:5;
    unsigned int __data_offset;
    unsigned short int __data_selector;
    unsigned short int __glibc_reserved5;
    // comment out below if 32bit
    unsigned int __mxcsr;
  }
  fenv_t;

  typedef unsigned short int fexcept_t;
%}
*/
/*
// fenv.h
int  feclearexcept(int);
int  fegetexceptflag(fexcept_t *, int);
int  feraiseexcept(int);
int  fesetexceptflag(const fexcept_t *, int);
int  fetestexcept(int);
int  fegetround(void);
int  fesetround(int);
int  fegetenv(fenv_t *);
int  feholdexcept(fenv_t *);
int  fesetenv(const fenv_t *);
int  feupdateenv(const fenv_t *);
*/
%constant int char_bit = CHAR_BIT;
%constant int schar_min = SCHAR_MIN;
%constant int schar_max = SCHAR_MAX;
%constant int uchar_max = UCHAR_MAX;
%constant int char_min = CHAR_MIN;
%constant int char_max = CHAR_MAX;
%constant int mb_len_max = MB_LEN_MAX;
%constant int shrt_min = SHRT_MIN;
%constant int shrt_max = SHRT_MAX;
%constant int ushrt_max = USHRT_MAX;
%constant int int_min = INT_MIN;
%constant int int_max = INT_MAX;
%constant int uint_max = UINT_MAX;
%constant int long_min = LONG_MIN;
%constant int long_max = LONG_MAX;
%constant int ulong_max = ULONG_MAX;
%constant int llong_min = LLONG_MIN;
%constant int llong_max = LLONG_MAX;
%constant int ullong_max = ULLONG_MAX;



