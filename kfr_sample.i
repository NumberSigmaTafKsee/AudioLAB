%module kfr_sample
%{
//#include "SampleVector.h"
#include "samples/kfr_sample.hpp"
#include "samples/kfr_sample_dsp.hpp"
%}

%include "stdint.i"
%include "std_vector.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(int8_vector) std::vector<signed char>;
%template(uint8_vector) std::vector<unsigned char>;
%template(int16_vector) std::vector<signed short>;
%template(uint16_vector) std::vector<unsigned short>;
%template(int32_vector) std::vector<signed int>;
%template(uint32_vector) std::vector<unsigned int>;
%template(int64_vector) std::vector<signed long>;
%template(uint64_vector) std::vector<unsigned long>;

//%include "SampleVector.h"
%include "samples/kfr_sample.hpp"
%include "samples/kfr_sample_dsp.hpp"

//%template(FloatSampleVector) DSP1::SampleVector<float>;

%template(get_left_channel_float) DSP1::get_left_channel<float>;
%template(get_right_channel_float) DSP1::get_right_channel<float>;
%template(get_channel_float) DSP1::get_channel<float>;

%template(interleave_float) DSP1::interleave<float>;
%template(deinterleave_float) DSP1::interleave<float>;
%template(copy_vector_float) DSP1::copy_vector<float>;
%template(slice_vector_float) DSP1::slice_vector<float>;
%template(copy_buffer_float) DSP1::copy_buffer<float>;
%template(slice_buffer_float) DSP1::slice_buffer<float>;
%template(stereo_split_float) DSP1::split_stereo<float>;
%template(insert_front_float) DSP1::insert_front<float>;

%template(containsOnlyZeros_float) DSP1::containsOnlyZeros<float>;
%template(isAllPositiveOrZero_float) DSP1::isAllPositiveOrZero<float>;
%template(isAllNegativeOrZero_float) DSP1::isAllNegativeOrZero<float>;
%template(contains_float) DSP1::contains<float>;
%template(max_float) DSP1::max<float>;
%template(min_float) DSP1::min<float>;
%template(maxIndex_float) DSP1::maxIndex<float>;
%template(minIndex_float) DSP1::minIndex<float>;
%template(printVector_float) DSP1::printVector<float>;
%template(getFirstElement_float) DSP1::getFirstElement<float>;
%template(getLastElement_float) DSP1::getLastElement<float>;
%template(getEvenElements_float) DSP1::getEvenElements<float>;
%template(getOddElements_float) DSP1::getOddElements<float>;
%template(getEveryNthElementStartingFromK_float) DSP1::getEveryNthElementStartingFromK<float>;
%template(fillVectorWith_float) DSP1::fillVectorWith<float>;
%template(countOccurrencesOf_float) DSP1::countOccurrencesOf<float>;
%template(sum_float) DSP1::sum<float>;
%template(product_float) DSP1::product<float>;
%template(mean_float) DSP1::mean<float>;
%template(median_float) DSP1::median<float>;
%template(variance_float) DSP1::variance<float>;
%template(standardDeviation_float) DSP1::standardDeviation<float>;
%template(norm1_float) DSP1::norm1<float>;
%template(norm2_float) DSP1::norm2<float>;
%template(normP_float) DSP1::normP<float>;
%template(magnitude_float) DSP1::magnitude<float>;
%template(multiplyInPlace_float) DSP1::multiplyInPlace<float>;
%template(divideInPlace_float) DSP1::divideInPlace<float>;
%template(addInPlace_float) DSP1::addInPlace<float>;
%template(subtractInPlace_float) DSP1::subtractInPlace<float>;
%template(absInPlace_float) DSP1::absInPlace<float>;
%template(squareInPlace_float) DSP1::squareInPlace<float>;
%template(squareRootInPlace_float) DSP1::squareRootInPlace<float>;
%template(sort_float) DSP1::sort<float>;
%template(reverse_float) DSP1::reverse<float>;
%template(multiply_float) DSP1::multiply<float>;
%template(divide_float) DSP1::divide<float>;
%template(add_float) DSP1::add<float>;
%template(subtract_float) DSP1::subtract<float>;
%template(abs_float) DSP1::abs<float>;
%template(square_float) DSP1::square<float>;
%template(squareRoot_float) DSP1::squareRoot<float>;
%template(scale_float) DSP1::scale<float>;
%template(difference_float) DSP1::difference<float>;
%template(zeros_float) DSP1::zeros<float>;
%template(ones_float) DSP1::ones<float>;
%template(range_float) DSP1::range<float>;
%template(dotProduct_float) DSP1::dotProduct<float>;
%template(euclideanDistance_float) DSP1::euclideanDistance<float>;
%template(cosineSimilarity_float) DSP1::cosineSimilarity<float>;
%template(cosineDistance_float) DSP1::cosineDistance<float>;
