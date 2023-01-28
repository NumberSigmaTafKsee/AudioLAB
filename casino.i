%module casino
%{
    #include <complex>
    #include <algorithm>
    #include <vector>
    #include <random>
    #include <functional>
    #include <memory>
    #include <chrono>
    #include <cmath>
    #include <cassert>

    #include <ippcore.h>
    #include <ipps.h>

    #include "carlo_samples.hpp"

    #include "carlo_ipp.hpp"
    #include "carlo_ipparray.hpp"
    #include "carlo_random.hpp"
    #include "carlo_autocorr.hpp"
    #include "carlo_xcorr.hpp"
    #include "carlo_convolution.hpp"
    #include "carlo_dct.hpp"
    #include "carlo_dft.hpp"
    #include "carlo_fft.hpp"
    #include "carlo_firlms.hpp"
    #include "carlo_firmr.hpp"
    #include "carlo_firsr.hpp"
    #include "carlo_hilbert.hpp"
    #include "carlo_iir.hpp"
    #include "carlo_resample.hpp"

    #define PRINT(foo) std::cout << foo << std::endl;

    typedef float  DspFloatType;
    typedef std::complex<DspFloatType> ComplexType;
    typedef Ipp32fc IppComplexType;

    using namespace Casino;
    using namespace Casino::IPP;
    using namespace Casino::MKL;
%}
typedef float  DspFloatType;
typedef std::complex<DspFloatType> ComplexType;
typedef Ipp32fc IppComplexType;


%include "std_math.i"
%include "std_vector.i"

%template(complex) std::complex<DspFloatType>;
typedef std::complex<DspFloatType> ComplexType;
%template(float_vector)  std::vector<float, cppmkl::cppmkl_allocator<float>>;
%template(double_vector) std::vector<double, cppmkl::cppmkl_allocator<double>>;
%template(complex_float_vector) std::vector<std::complex<float>, cppmkl::cppmkl_allocator<std::complex<float>>>;
%template(complex_double_vector) std::vector<std::complex<double>, cppmkl::cppmkl_allocator<std::complex<double>>>;

%include "carlo_vector.hpp"
%include "carlo_complex_vector.hpp"
%include "carlo_matrix.hpp"
%include "carlo_complex_matrix.hpp"
%include "carlo_mklfft.hpp"

%include "carlo_samples.hpp"
%include "carlo_ipparray.hpp"
%include "carlo_random.hpp"
%include "carlo_autocorr.hpp"
%include "carlo_xcorr.hpp"
%include "carlo_convolution.hpp"
%include "carlo_dct.hpp"
%include "carlo_dft.hpp"
%include "carlo_fft.hpp"
%include "carlo_firlms.hpp"
%include "carlo_firmr.hpp"
%include "carlo_firsr.hpp"
%include "carlo_hilbert.hpp"
%include "carlo_iir.hpp"
%include "carlo_resample.hpp"

%typemap(in) float* (std::complex<float> *ptr)
{
    $1 = (float*)ptr;
}
%typemap(in) double* (std::complex<double> *ptr)
{
    $1 = (double*)ptr;
}
%typemap(in) Ipp32f* (float *ptr)
{
    $1 = (float*)ptr;
}
%typemap(in) Ipp64f* (double *ptr)
{
    $1 = (double*)ptr;
}
%typemap(in) Ipp32fc* (std::complex<float> *ptr)
{
    $1 = (Ipp32fc*)ptr;
}
%typemap(in) Ipp64fc* (std::complex<double> *ptr)
{
    $1 = (Ipp64fc*)ptr;
}

%template(array) Casino::IPP::IPPArray<DspFloatType>;
%template(vector) Casino::MKL::Vector<DspFloatType>;
%template(complex_vector) Casino::MKL::ComplexVector<DspFloatType>;
%template(matrix_view) Casino::MKL::MatrixView<DspFloatType>;
%template(matrix) Casino::MKL::Matrix<DspFloatType>;
%template(complex_matrix_view) Casino::MKL::ComplexMatrixView<DspFloatType>;
%template(complex_matrix) Casino::MKL::ComplexMatrix<DspFloatType>;

%template(autocorr) Casino::IPP::AutoCorrelation<DspFloatType>;
%template(crosscorr) Casino::IPP::CrossCorrelation<DspFloatType>;
%template(convolver) Casino::IPP::Convolver<DspFloatType>;
%template(convolution_filter) Casino::IPP::ConvolutionFilter<DspFloatType>;

%template(DCT) Casino::IPP::DCT<DspFloatType>;
%template(FIRMR) Casino::IPP::FIRMR<DspFloatType>;
%template(FIRSR) Casino::IPP::FIRSR<DspFloatType>;
%template(IIR) Casino::IPP::IIR<DspFloatType>;
%template(BIQUAD) Casino::IPP::IIRBiquad<DspFloatType>;
%template(random_uniform) Casino::IPP::RandomUniform<DspFloatType>;


%inline
%{
    void DisplayIPP()
    {
        const IppLibraryVersion *vers;
        int cache_size;	
        int mhz;
        int bCache;
        int numThreads;
        
        Ipp64u cpu_clocks;
        IppCache * pCache;
        vers = ippsGetLibVersion();
        ippGetL2CacheSize(&cache_size);
        ippGetCacheParams(&pCache);	
        cpu_clocks = ippGetCpuClocks();
        ippGetCpuFreqMhz(&mhz);
        ippGetMaxCacheSizeB(&bCache);
        // doesn't seem to do anything
        ippSetNumThreads(8);        
        ippGetNumThreads(&numThreads);
        

        PRINT("--------------------------------------")
        PRINT("IPP Library Version");
        PRINT("Minor Version: " << vers->minor)
        PRINT("Major Version: " << vers->major)
        PRINT("Major Build  : " << vers->majorBuild)
        PRINT("Build        : " << vers->build);
        PRINT("TargetCPU[0] : " << vers->targetCpu[0]);
        PRINT("TargetCPU[1] : " << vers->targetCpu[1]);
        PRINT("TargetCPU[2] : " << vers->targetCpu[2]);
        PRINT("TargetCPU[3] : " << vers->targetCpu[3]);
        PRINT("Name         : " << vers->Name);
        PRINT("Version      : " << vers->Version);
        PRINT("Build Date   : " << vers->BuildDate);
        PRINT("--------------------------------------")

        
        PRINT("L2 Cache Size\t: " << cache_size);
        PRINT("Cache Type   \t: " << pCache->type);
        PRINT("Cache Level  \t: " << pCache->level);
        PRINT("Cache Size   \t: " << pCache->size);
        PRINT("CPU Clocks   \t: " << cpu_clocks);
        PRINT("CPU Mhz      \t: " << mhz);
        PRINT("B-Cache Size \t: " << bCache);
        PRINT("Num-Threads  \t: " << numThreads);
    }
%}


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