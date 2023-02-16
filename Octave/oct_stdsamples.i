%module stdsamples
%{
#define DSPFLOATDOUBLE
#include "SoundObject.hpp"    
#include "audiodsp_std_samples.hpp"

using namespace AudioDSP;
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

#define DSPFLOATDOUBLE
%include "SoundObject.hpp"

%template(complex) std::complex<DspFloatType>;


%template(float_vector)          std::vector<float,Allocator::aligned_allocator<float,64>>;
%template(double_vector)         std::vector<double,Allocator::aligned_allocator<double,64>>;
%template(complex_float_vector)  std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>>;
%template(complex_double_vector) std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>>;

%template(int8_vector) std::vector<signed char>;
%template(uint8_vector) std::vector<unsigned char>;
%template(int16_vector) std::vector<signed short>;
%template(uint16_vector) std::vector<unsigned short>;
%template(int32_vector) std::vector<signed int>;
%template(uint32_vector) std::vector<unsigned int>;
%template(int64_vector) std::vector<signed long>;
%template(uint64_vector) std::vector<unsigned long>;


%include "audiodsp_std_sample_vector.hpp"
%include "audiodsp_std_sample_matrix.hpp"
%include "audiodsp_std_complex_vector.hpp"
%include "audiodsp_std_complex_matrix.hpp"
%include "audiodsp_std_samples.hpp"

%template(get_left_channel)  AudioDSP::get_left_channel<DspFloatType>;
%template(get_right_channel) AudioDSP::get_right_channel<DspFloatType>;
%template(get_channel)       AudioDSP::get_channel<DspFloatType>;

%template(interleave)     AudioDSP::interleave<DspFloatType>;
%template(deinterleave)   AudioDSP::interleave<DspFloatType>;
%template(copy_vector)    AudioDSP::copy_vector<DspFloatType>;
%template(slice_vector)   AudioDSP::slice_vector<DspFloatType>;
%template(copy_buffer)    AudioDSP::copy_buffer<DspFloatType>;
%template(slice_buffer)   AudioDSP::slice_buffer<DspFloatType>;
%template(stereo_split)   AudioDSP::split_stereo<DspFloatType>;
//%template(insert_front)   AudioDSP::insert_front<DspFloatType>;

%template(fill)   AudioDSP::fill<DspFloatType>;
%template(zeros)  AudioDSP::zeros<DspFloatType>;
%template(ones)   AudioDSP::ones<DspFloatType>;


%template(sample_vector) AudioDSP::sample_vector<DspFloatType>;
%template(complex_vector) AudioDSP::complex_vector<DspFloatType>;

%template(sample_matrix) AudioDSP::sample_matrix<DspFloatType>;
%template(complex_matrix) AudioDSP::complex_matrix<DspFloatType>;

//%template(window) AudioDSP::Window<DspFloatType>;

//%template(abs)  AudioDSP::abs<DspFloatType>;
%template(cube) AudioDSP::cube<DspFloatType>;
%template(sqr) AudioDSP::sqr<DspFloatType>;
%template(sqrt) AudioDSP::sqrt<DspFloatType>;
%template(exp)  AudioDSP::exp<DspFloatType>;
%template(exp2) AudioDSP::exp2<DspFloatType>;
%template(log)  AudioDSP::log<DspFloatType>;
%template(log10) AudioDSP::log10<DspFloatType>;
%template(log2) AudioDSP::log2<DspFloatType>;
%template(logb) AudioDSP::logb<DspFloatType>;
%template(pow) AudioDSP::pow<DspFloatType>;
%template(floor) AudioDSP::floor<DspFloatType>;
%template(acos) AudioDSP::acos<DspFloatType>;
%template(asin) AudioDSP::asin<DspFloatType>;
%template(atan) AudioDSP::atan<DspFloatType>;
%template(atan2) AudioDSP::atan2<DspFloatType>;
%template(cos) AudioDSP::cos<DspFloatType>;
%template(sin) AudioDSP::sin<DspFloatType>;
%template(tan) AudioDSP::tan<DspFloatType>;
%template(cosh) AudioDSP::cosh<DspFloatType>;
%template(sinh) AudioDSP::sinh<DspFloatType>;
%template(tanh) AudioDSP::tanh<DspFloatType>;
%template(lgamma) AudioDSP::lgamma<DspFloatType>;
%template(acosh) AudioDSP::acosh<DspFloatType>;
%template(asinh) AudioDSP::asinh<DspFloatType>;
%template(atanh) AudioDSP::atanh<DspFloatType>;
%template(cbrt) AudioDSP::cbrt<DspFloatType>;
%template(ceil) AudioDSP::cbrt<DspFloatType>;
%template(copysign) AudioDSP::copysign<DspFloatType>;
%template(erf) AudioDSP::erf<DspFloatType>;
%template(erfc) AudioDSP::erfc<DspFloatType>;
%template(expm1) AudioDSP::expm1<DspFloatType>;
%template(fdim) AudioDSP::fdim<DspFloatType>;
%template(fma) AudioDSP::fma<DspFloatType>;
%template(fmax) AudioDSP::fmax<DspFloatType>;
%template(fmin) AudioDSP::fmin<DspFloatType>;
%template(fmod) AudioDSP::fmod<DspFloatType>;
%template(hypot) AudioDSP::hypot<DspFloatType>;

%template(lgamma) AudioDSP::lgamma<DspFloatType>;
%template(remainder) AudioDSP::remainder<DspFloatType>;
%template(round) AudioDSP::round<DspFloatType>;
%template(scalbln) AudioDSP::scalbln<DspFloatType>;
%template(scalbn) AudioDSP::scalbn<DspFloatType>;
%template(tgamma) AudioDSP::tgamma<DspFloatType>;
%template(trunc) AudioDSP::trunc<DspFloatType>;
%template(ilogb) AudioDSP::ilogb<DspFloatType>;

%template(absf)  Ops::abs<DspFloatType>;
%template(cubef) Ops::cube<DspFloatType>;
%template(sqrtf) Ops::sqrt<DspFloatType>;
%template(expf)  Ops::exp<DspFloatType>;
%template(exp2f) Ops::exp2<DspFloatType>;
%template(logf)  Ops::log<DspFloatType>;
%template(log10f) Ops::log10<DspFloatType>;
%template(log2f) Ops::log2<DspFloatType>;
%template(logbf) Ops::logb<DspFloatType>;
%template(powf) Ops::pow<DspFloatType>;
%template(floorf) Ops::floor<DspFloatType>;
%template(acosf) Ops::acos<DspFloatType>;
%template(asinf) Ops::asin<DspFloatType>;
%template(atanf) Ops::atan<DspFloatType>;
%template(atan2f) Ops::atan2<DspFloatType>;
%template(cosf) Ops::cos<DspFloatType>;
%template(sinf) Ops::sin<DspFloatType>;
%template(tanf) Ops::tan<DspFloatType>;
%template(coshf) Ops::cosh<DspFloatType>;
%template(sinhf) Ops::sinh<DspFloatType>;
%template(tanhf) Ops::tanh<DspFloatType>;
%template(lgammaf) Ops::lgamma<DspFloatType>;
%template(acoshf) Ops::acosh<DspFloatType>;
%template(asinhf) Ops::asinh<DspFloatType>;
%template(atanhf) Ops::atanh<DspFloatType>;
%template(cbrtf) Ops::cbrt<DspFloatType>;
%template(ceilf) Ops::cbrt<DspFloatType>;
%template(copysignf) Ops::copysign<DspFloatType>;
%template(erff) Ops::erf<DspFloatType>;
%template(erfcf) Ops::erfc<DspFloatType>;
%template(expm1f) Ops::expm1<DspFloatType>;
%template(fdimf) Ops::fdim<DspFloatType>;
%template(fmaf) Ops::fma<DspFloatType>;
%template(fmaxf) Ops::fmax<DspFloatType>;
%template(fminf) Ops::fmin<DspFloatType>;
%template(fmodf) Ops::fmod<DspFloatType>;
%template(fpclassifyf) Ops::fpclassify<DspFloatType>;
%template(hypotf) Ops::hypot<DspFloatType>;
%template(ilogbf) Ops::ilogb<DspFloatType>;
%template(isfinitef) Ops::isfinite<DspFloatType>;
%template(isgreaterf) Ops::isgreater<DspFloatType>;
%template(isgreaterequalf) Ops::isgreaterequal<DspFloatType>;
%template(isinff) Ops::isinf<DspFloatType>;
%template(islessf) Ops::isless<DspFloatType>;
%template(islessequalf) Ops::islessequal<DspFloatType>;
%template(isnanf) Ops::isnan<DspFloatType>;
%template(isnormalf) Ops::isnormal<DspFloatType>;
%template(isunorderedf) Ops::isunordered<DspFloatType>;
%template(ldexpf) Ops::ldexp<DspFloatType>;
%template(lgammaf) Ops::lgamma<DspFloatType>;
%template(llrintf) Ops::llrint<DspFloatType>;
%template(llroundf) Ops::llround<DspFloatType>;
%template(log1pf) Ops::log1p<DspFloatType>;
%template(lrintf) Ops::lrint<DspFloatType>;
%template(lroundf) Ops::lround<DspFloatType>;
%template(nanf) Ops::nan<DspFloatType>;
%template(nanff) Ops::nanf<DspFloatType>;
%template(nanlf) Ops::nanl<DspFloatType>;
%template(nearbyintf) Ops::nearbyint<DspFloatType>;
%template(nextafterf) Ops::nextafter<DspFloatType>;
%template(nexttowardf) Ops::nexttoward<DspFloatType>;
%template(remainderf) Ops::remainder<DspFloatType>;
%template(rintf) Ops::rint<DspFloatType>;
%template(roundf) Ops::round<DspFloatType>;
%template(scalblnf) Ops::scalbln<DspFloatType>;
%template(scalbnf) Ops::scalbn<DspFloatType>;
%template(squaref) Ops::square<DspFloatType>;
%template(tgammaf) Ops::tgamma<DspFloatType>;
%template(truncf) Ops::trunc<DspFloatType>;


%template(crealf) std::real<DspFloatType>;
%template(cimagf) std::imag<DspFloatType>;
%template(cabsf) std::abs<DspFloatType>;
%template(cargf) std::arg<DspFloatType>;
%template(cexpf) std::exp<DspFloatType>;
%template(clogf) std::log<DspFloatType>;
%template(clog10f) std::log10<DspFloatType>;
%template(cpowf) std::pow<DspFloatType>;
%template(csqrtf) std::sqrt<DspFloatType>;
%template(cnormf) std::norm<DspFloatType>;
%template(cprojf) std::proj<DspFloatType>;
%template(cpolarf) std::polar<DspFloatType>;
%template(csinf) std::sin<DspFloatType>;
%template(ccosf) std::cos<DspFloatType>;
%template(ctanf) std::tan<DspFloatType>;
%template(casinf) std::asin<DspFloatType>;
%template(cacosf) std::acos<DspFloatType>;
%template(catanf) std::atan<DspFloatType>;
%template(csinhf) std::sinh<DspFloatType>;
%template(ccoshf) std::cosh<DspFloatType>;
%template(ctanhf) std::tanh<DspFloatType>;
%template(casinhf) std::asinh<DspFloatType>;
%template(cacoshf) std::acosh<DspFloatType>;
%template(catanhf) std::atanh<DspFloatType>;


%typemap(in) FloatRowVector(float *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_row_vector_value()(i) = in[i];
%}
%typemap(in) RowVector(double *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->row_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexRowVector(std::complex<float> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexRowVector(std::complex<double> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}

%typemap(in) FloatColumnVector(float *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_column_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexColumnVector(std::complex<float> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->float_complex_column_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexColumnVector(std::complex<double> *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ColumnVector(double *in, size_t n) %{
    for (size_t i = 0; i < n; i++)
        $1->column_vector_value()(i) = in[i];
    %}

%typemap(in) FloatRowVector(std::vector<float,Allocator::aligned_allocator<float,64>> &in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_row_vector_value()(i) = in[i];
    %}
%typemap(in) RowVector(std::vector<double,Allocator::aligned_allocator<double,64>> &in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->row_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexRowVector(std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> & in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexRowVector(std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> &in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}

%typemap(in) FloatColumnVector(std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> & in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_complex_column_vector_value()(i) = in[i];
    %}
%typemap(in) ColumnVector(std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> &in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->column_complex_vector_value()(i) = in[i];
    %}
%typemap(in) FloatComplexColumnVector(std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> & in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->float_complex_row_vector_value()(i) = in[i];
    %}
%typemap(in) ComplexColumnVector(std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> &in) %{
    for (size_t i = 0; i < in.size(); i++)
        $1->complex_row_vector_value()(i) = in[i];
    %}

%typemap(out) std::vector<float,Allocator::aligned_allocator<float,64>> in %{
    $result = std::vector<float,Allocator::aligned_allocator<float,64>>(in.size());
    for (size_t i = 0; i < in.size(1); i++)
        $result->float_row_vector()(i) = in(i);
    %}
%typemap(out) std::vector<double,Allocator::aligned_allocator<double,64>> in %{
    $result = std::vector<double,Allocator::aligned_allocator<double,64>>(in.size(1));
    for (size_t i = 0; i < in.size(1); i++)
        $result->row_vector()(i) = in(i);
    %}

%typemap(in) FloatMatrix(float *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->float_matrix_value()(i, j) = in[i * n + j];
    %}
%typemap(in) Matrix(double *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->matrix_value()(i, j) = in[i * n + j];
    %}
%typemap(in) FloatComplexMatrix(std::complex<float> *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->float_complex_matrix_value()(i, j) = in[i * n + j];
    %}
%typemap(in) ComplexMatrix(std::complex<double> *in, size_t m, size_t n) %{
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
            $1->complex_matrix_value()(i, j) = in[i * n + j];
    %}

%inline %{

    octave_value vectorize_row(const std::vector<float,Allocator::aligned_allocator<float,64>> &src)
    {
        FloatRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_row(const float *src, size_t n)
    {
        FloatRowVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<float,Allocator::aligned_allocator<float,64>> 
    vectorize_row_float_vector(const octave_value &src)
    {
        FloatRowVector tmp = src.float_row_vector_value();
        std::vector<float,Allocator::aligned_allocator<float,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_row(const std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> &src)
    {
        FloatComplexRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_row(const std::complex<float> *src, size_t n)
    {
        FloatComplexRowVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> 
    vectorize_row_complex_float_vector(const octave_value &src)
    {
        FloatComplexRowVector tmp = src.float_complex_row_vector_value();
        std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_row(const std::vector<double,Allocator::aligned_allocator<double,64>> &src)
    {
        RowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_row(const double *src, size_t n)
    {
        RowVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<double,Allocator::aligned_allocator<double,64>> 
    vectorize_row_double_vector(const octave_value &src)
    {
        RowVector tmp = src.row_vector_value();
        std::vector<double,Allocator::aligned_allocator<double,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_row(const std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> &src)
    {
        ComplexRowVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_row(const std::complex<double> *src, size_t n)
    {
        ComplexRowVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> 
    vectorize_complex_row_double_vector(const octave_value &src)
    {
        ComplexRowVector tmp = src.complex_row_vector_value();
        std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_col(const std::vector<float,Allocator::aligned_allocator<float,64>> &src)
    {
        FloatColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_col(const float *src, size_t n)
    {
        FloatColumnVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<float,Allocator::aligned_allocator<float,64>> 
    vectorize_col_float_vector(const octave_value &src)
    {
        FloatColumnVector tmp = src.float_column_vector_value();
        std::vector<float,Allocator::aligned_allocator<float,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_col(const std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> &src)
    {
        FloatComplexColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_col(const std::complex<float> *src, size_t n)
    {
        FloatComplexColumnVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> 
    vectorize_complex_col_float_vector(const octave_value &src)
    {
        FloatComplexColumnVector tmp = src.float_complex_column_vector_value();
        std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_col(const std::vector<double,Allocator::aligned_allocator<double,64>> &src)
    {
        ColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_col(const double *src, size_t n)
    {
        ColumnVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<double,Allocator::aligned_allocator<double,64>> 
    vectorize_col_double_vector(const octave_value &src)
    {
        ColumnVector tmp = src.column_vector_value();
        std::vector<double,Allocator::aligned_allocator<double,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value vectorize_complex_col(const std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> &src)
    {
        ComplexColumnVector dst(src.size());
        for (size_t i = 0; i < src.size(); i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    octave_value vectorize_complex_col(const std::complex<double> *src, size_t n)
    {
        ComplexColumnVector dst(n);
        for (size_t i = 0; i < n; i++)
            dst(i) = src[i];
        return octave_value(dst);
    }
    std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> 
    vectorize_complex_col_double_vector(const octave_value &src)
    {
        ComplexColumnVector tmp = src.complex_column_vector_value();
        std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> dst(tmp.size(1));
        for (size_t i = 0; i < tmp.size(1); i++)
            dst[i] = tmp(i);
        return dst;
    }

    octave_value matricize(const std::vector<float,Allocator::aligned_allocator<float,64>> &src, size_t m, size_t n)
    {
        FloatMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value matricize(const float *src, size_t m, size_t n)
    {
        FloatMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<float,Allocator::aligned_allocator<float,64>> 
    matricize_float_vector(const octave_value &src)
    {
        FloatMatrix tmp = src.float_matrix_value();
        std::vector<float,Allocator::aligned_allocator<float,64>> dst(tmp.rows() * tmp.cols());
        for (size_t i = 0; i < tmp.rows(); i++)
            for (size_t j = 0; j < tmp.cols(); j++)
                dst[i * tmp.cols() + j] = tmp(i, j);
        return dst;
    }

    octave_value matricize(const std::vector<double,Allocator::aligned_allocator<double,64>> &src, size_t m, size_t n)
    {
        Matrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value matricize(const double *src, size_t m, size_t n)
    {
        Matrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<double,Allocator::aligned_allocator<double,64>> 
    matricize_double_vector(const octave_value &src)
    {
        Matrix tmp = src.matrix_value();
        std::vector<double,Allocator::aligned_allocator<double,64>> dst(tmp.rows() * tmp.cols());
        for (size_t i = 0; i < tmp.rows(); i++)
            for (size_t j = 0; j < tmp.cols(); j++)
                dst[i * tmp.cols() + j] = tmp(i, j);
        return dst;
    }


    octave_value complex_matricize(const std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> &src, size_t m, size_t n)
    {
        FloatComplexMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value complex_matricize(const std::complex<float> *src, size_t m, size_t n)
    {
        FloatComplexMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> 
    complex_matricize_float_vector(const octave_value &src)
    {
        FloatComplexMatrix tmp = src.float_matrix_value();
        std::vector<std::complex<float>,Allocator::aligned_allocator<std::complex<float>,64>> dst(tmp.rows() * tmp.cols());
        for (size_t i = 0; i < tmp.rows(); i++)
            for (size_t j = 0; j < tmp.cols(); j++)
                dst[i * tmp.cols() + j] = tmp(i, j);
        return dst;
    }
    octave_value complex_matricize(const std::vector < std::complex<double>> & src, size_t m, size_t n)
    {
        ComplexMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    octave_value complex_matricize(const std::complex<double> *src, size_t m, size_t n)
    {
        ComplexMatrix dst(m, n);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                dst(i, j) = src[i * n + j];
        return octave_value(dst);
    }
    std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> 
    complex_matricize_double_vector(const octave_value &src)
    {
        ComplexMatrix tmp = src.complex_matrix_value();
        std::vector<std::complex<double>,Allocator::aligned_allocator<std::complex<double>,64>> dst(tmp.rows() * tmp.cols());
        for (size_t i = 0; i < tmp.rows(); i++)
            for (size_t j = 0; j < tmp.cols(); j++)
                dst[i * tmp.cols() + j] = tmp(i, j);
        return dst;
    }
%}