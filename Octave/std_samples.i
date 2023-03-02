//%module stdsamples
%{
#include "SoundObject.hpp"    
#include "audiodsp_std_samples.hpp"
#include "audiodsp_vtk.hpp"
#include "audiodsp_sndfile.hpp"
#include "audiodsp_fft.hpp"
#include "audiodsp_fftw.hpp"
#include "audiodsp_fftw_convolution.hpp"
#include "audiodsp_samples_functions.hpp"
#include "audiodsp_window.hpp"

using namespace AudioDSP;
%}


%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"




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

%include "SoundObject.hpp"

%include "audiodsp_std_sample_vector.hpp"
%include "audiodsp_std_sample_matrix.hpp"
%include "audiodsp_std_complex_vector.hpp"
%include "audiodsp_std_complex_matrix.hpp"
%include "audiodsp_std_samples.hpp"
%include "audiodsp_sndfile.hpp"
%include "audiodsp_fft.hpp"
%include "audiodsp_fftw.hpp"
%include "audiodsp_fftw_convolution.hpp"
%include "audiodsp_window.hpp"

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

%template(complex) std::complex<DspFloatType>;

%template(sample_vector) AudioDSP::sample_vector<DspFloatType>;
%template(complex_vector) AudioDSP::complex_vector<DspFloatType>;

%template(sample_matrix) AudioDSP::sample_matrix<DspFloatType>;
%template(complex_matrix) AudioDSP::complex_matrix<DspFloatType>;

%template(window) AudioDSP::Window<DspFloatType>;

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

%include "audiodsp_samples_functions.hpp"
%template(generate_noise) AudioDSP::generate_noise<DspFloatType>;
%template(generate_sin) AudioDSP::generate_sin<DspFloatType>;
%template(generate_cos) AudioDSP::generate_cos<DspFloatType>;
%template(generate_tan) AudioDSP::generate_tan<DspFloatType>;
%template(generate_phasor) AudioDSP::generate_phasor<DspFloatType>;
%template(generate_square) AudioDSP::generate_square<DspFloatType>;
%template(generate_saw) AudioDSP::generate_saw<DspFloatType>;    
%template(generate_triangle) AudioDSP::generate_triangle<DspFloatType>;
%template(oscillator) AudioDSP::oscillator<DspFloatType>;
%template(generator) AudioDSP::generator<DspFloatType>;
%template(filter) AudioDSP::filter<DspFloatType>;
%template(function) AudioDSP::function<DspFloatType>;    
    
/*
%template(fpclassify) AudioDSP::fpclassify<DspFloatType>;
%template(ldexp) AudioDSP::ldexp<DspFloatType>;

%template(isfinite) AudioDSP::isfinite<DspFloatType>;
%template(isgreater) AudioDSP::isgreater<DspFloatType>;
%template(isgreaterequal) AudioDSP::isgreaterequal<DspFloatType>;
%template(isinf) AudioDSP::isinf<DspFloatType>;
%template(isless) AudioDSP::isless<DspFloatType>;
%template(islessequal) AudioDSP::islessequal<DspFloatType>;
%template(isnan) AudioDSP::isnan<DspFloatType>;
%template(isnormal) AudioDSP::isnormal<DspFloatType>;
%template(isunordered) AudioDSP::isunordered<DspFloatType>;
%template(llrint) AudioDSP::llrint<DspFloatType>;
%template(llround) AudioDSP::llround<DspFloatType>;
%template(log1p) AudioDSP::log1p<DspFloatType>;
%template(lrint) AudioDSP::lrint<DspFloatType>;
%template(lround) AudioDSP::lround<DspFloatType>;
%template(nan) AudioDSP::nan<DspFloatType>;
%template(nanf) AudioDSP::nanf<DspFloatType>;
%template(nanl) AudioDSP::nanl<DspFloatType>;
%template(nearbyint) AudioDSP::nearbyint<DspFloatType>;
%template(nextafter) AudioDSP::nextafter<DspFloatType>;
%template(nexttoward) AudioDSP::nexttoward<DspFloatType>;
%template(rint) AudioDSP::rint<DspFloatType>;
*/

/*
%include "audiodsp_vtk.hpp"

%template(containsOnlyZeros)   AudioDSP::containsOnlyZeros<DspFloatType>;
%template(isAllPositiveOrZero) isAllPositiveOrZero<DspFloatType>;
%template(isAllNegativeOrZero) isAllNegativeOrZero<DspFloatType>;
%template(contains) contains<DspFloatType>;
%template(max) max<DspFloatType>;
%template(min) min<DspFloatType>;
%template(maxIndex) maxIndex<DspFloatType>;
%template(minIndex) minIndex<DspFloatType>;
%template(printVector) printVector<DspFloatType>;
%template(getFirstElement) getFirstElement<DspFloatType>;
%template(getLastElement) getLastElement<DspFloatType>;
%template(getEvenElements) getEvenElements<DspFloatType>;
%template(getOddElements) getOddElements<DspFloatType>;
%template(getEveryNthElementStartingFromK) getEveryNthElementStartingFromK<DspFloatType>;
%template(fillVectorWith) fillVectorWith<DspFloatType>;
%template(countOccurrencesOf) countOccurrencesOf<DspFloatType>;
%template(sum) sum<DspFloatType>;
%template(product) product<DspFloatType>;
%template(mean) mean<DspFloatType>;
%template(median) median<DspFloatType>;
%template(variance) variance<DspFloatType>;
%template(standardDeviation) standardDeviation<DspFloatType>;
%template(norm1) norm1<DspFloatType>;
%template(norm2) norm2<DspFloatType>;
%template(normP) normP<DspFloatType>;
%template(magnitude) magnitude<DspFloatType>;
%template(multiplyInPlace) multiplyInPlace<DspFloatType>;
%template(divideInPlace) divideInPlace<DspFloatType>;
%template(addInPlace) addInPlace<DspFloatType>;
%template(subtractInPlace) subtractInPlace<DspFloatType>;
%template(absInPlace) absInPlace<DspFloatType>;
%template(squareInPlace) squareInPlace<DspFloatType>;
%template(squareRootInPlace) squareRootInPlace<DspFloatType>;
%template(sort) sort<DspFloatType>;
%template(reverse) reverse<DspFloatType>;
%template(multiply) multiply<DspFloatType>;
%template(divide) divide<DspFloatType>;
%template(add) add<DspFloatType>;
%template(subtract) subtract<DspFloatType>;
%template(abs) abs<DspFloatType>;
%template(square) square<DspFloatType>;
%template(squareRoot) squareRoot<DspFloatType>;
%template(scale) scale<DspFloatType>;
%template(difference) difference<DspFloatType>;
%template(zeros) zeros<DspFloatType>;
%template(ones) ones<DspFloatType>;
%template(range) range<DspFloatType>;
%template(dotProduct) dotProduct<DspFloatType>;
%template(euclideanDistance) euclideanDistance<DspFloatType>;
%template(cosineSimilarity) cosineSimilarity<DspFloatType>;
%template(cosineDistance) cosineDistance<DspFloatType>;
*/