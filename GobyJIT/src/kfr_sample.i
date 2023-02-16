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

//%template(FloatSampleVector) KfrDSP1::SampleVector<float>;

%template(get_left_channel_float) KfrDSP1::get_left_channel<float>;
%template(get_right_channel_float) KfrDSP1::get_right_channel<float>;
%template(get_channel_float) KfrDSP1::get_channel<float>;

%template(interleave_float) KfrDSP1::interleave<float>;
%template(deinterleave_float) KfrDSP1::interleave<float>;
%template(copy_vector_float) KfrDSP1::copy_vector<float>;
%template(slice_vector_float) KfrDSP1::slice_vector<float>;
%template(copy_buffer_float) KfrDSP1::copy_buffer<float>;
%template(slice_buffer_float) KfrDSP1::slice_buffer<float>;
%template(stereo_split_float) KfrDSP1::split_stereo<float>;
%template(insert_front_float) KfrDSP1::insert_front<float>;

%template(containsOnlyZeros_float) KfrDSP1::containsOnlyZeros<float>;
%template(isAllPositiveOrZero_float) KfrDSP1::isAllPositiveOrZero<float>;
%template(isAllNegativeOrZero_float) KfrDSP1::isAllNegativeOrZero<float>;
%template(contains_float) KfrDSP1::contains<float>;
%template(max_float) KfrDSP1::max<float>;
%template(min_float) KfrDSP1::min<float>;
%template(maxIndex_float) KfrDSP1::maxIndex<float>;
%template(minIndex_float) KfrDSP1::minIndex<float>;
%template(printVector_float) KfrDSP1::printVector<float>;
%template(getFirstElement_float) KfrDSP1::getFirstElement<float>;
%template(getLastElement_float) KfrDSP1::getLastElement<float>;
%template(getEvenElements_float) KfrDSP1::getEvenElements<float>;
%template(getOddElements_float) KfrDSP1::getOddElements<float>;
%template(getEveryNthElementStartingFromK_float) KfrDSP1::getEveryNthElementStartingFromK<float>;
%template(fillVectorWith_float) KfrDSP1::fillVectorWith<float>;
%template(countOccurrencesOf_float) KfrDSP1::countOccurrencesOf<float>;
%template(sum_float) KfrDSP1::sum<float>;
%template(product_float) KfrDSP1::product<float>;
%template(mean_float) KfrDSP1::mean<float>;
%template(median_float) KfrDSP1::median<float>;
%template(variance_float) KfrDSP1::variance<float>;
%template(standardDeviation_float) KfrDSP1::standardDeviation<float>;
%template(norm1_float) KfrDSP1::norm1<float>;
%template(norm2_float) KfrDSP1::norm2<float>;
%template(normP_float) KfrDSP1::normP<float>;
%template(magnitude_float) KfrDSP1::magnitude<float>;
%template(multiplyInPlace_float) KfrDSP1::multiplyInPlace<float>;
%template(divideInPlace_float) KfrDSP1::divideInPlace<float>;
%template(addInPlace_float) KfrDSP1::addInPlace<float>;
%template(subtractInPlace_float) KfrDSP1::subtractInPlace<float>;
%template(absInPlace_float) KfrDSP1::absInPlace<float>;
%template(squareInPlace_float) KfrDSP1::squareInPlace<float>;
%template(squareRootInPlace_float) KfrDSP1::squareRootInPlace<float>;
%template(sort_float) KfrDSP1::sort<float>;
%template(reverse_float) KfrDSP1::reverse<float>;
%template(multiply_float) KfrDSP1::multiply<float>;
%template(divide_float) KfrDSP1::divide<float>;
%template(add_float) KfrDSP1::add<float>;
%template(subtract_float) KfrDSP1::subtract<float>;
%template(abs_float) KfrDSP1::abs<float>;
%template(square_float) KfrDSP1::square<float>;
%template(squareRoot_float) KfrDSP1::squareRoot<float>;
%template(scale_float) KfrDSP1::scale<float>;
%template(difference_float) KfrDSP1::difference<float>;
%template(zeros_float) KfrDSP1::zeros<float>;
%template(ones_float) KfrDSP1::ones<float>;
%template(range_float) KfrDSP1::range<float>;
%template(dotProduct_float) KfrDSP1::dotProduct<float>;
%template(euclideanDistance_float) KfrDSP1::euclideanDistance<float>;
%template(cosineSimilarity_float) KfrDSP1::cosineSimilarity<float>;
%template(cosineDistance_float) KfrDSP1::cosineDistance<float>;
