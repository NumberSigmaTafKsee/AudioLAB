%module samples
%{
#include "audiodsp_samples.hpp"
#include "audiodsp_samplesdsp.hpp"
#include "audiodsp_sndfile.hpp"
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

%include "audiodsp_samplevector.hpp"
%include "audiodsp_samples.hpp"
%include "audiodsp_sampledsp.hpp"
%include "audiodsp_sndfile.hpp"

%template(get_left_channel_float) get_left_channel<float>;
%template(get_right_channel_float) get_right_channel<float>;
%template(get_channel_float) get_channel<float>;

%template(interleave_float) interleave<float>;
%template(deinterleave_float) interleave<float>;
%template(copy_vector_float) copy_vector<float>;
%template(slice_vector_float) slice_vector<float>;
%template(copy_buffer_float) copy_buffer<float>;
%template(slice_buffer_float) slice_buffer<float>;
%template(stereo_split_float) split_stereo<float>;
//%template(insert_front_float) insert_front<float>;

%template(containsOnlyZeros_float) AudioDSP::containsOnlyZeros<float>;
%template(isAllPositiveOrZero_float) isAllPositiveOrZero<float>;
%template(isAllNegativeOrZero_float) isAllNegativeOrZero<float>;
%template(contains_float) contains<float>;
%template(max_float) max<float>;
%template(min_float) min<float>;
%template(maxIndex_float) maxIndex<float>;
%template(minIndex_float) minIndex<float>;
%template(printVector_float) printVector<float>;
%template(getFirstElement_float) getFirstElement<float>;
%template(getLastElement_float) getLastElement<float>;
%template(getEvenElements_float) getEvenElements<float>;
%template(getOddElements_float) getOddElements<float>;
%template(getEveryNthElementStartingFromK_float) getEveryNthElementStartingFromK<float>;
%template(fillVectorWith_float) fillVectorWith<float>;
%template(countOccurrencesOf_float) countOccurrencesOf<float>;
%template(sum_float) sum<float>;
%template(product_float) product<float>;
%template(mean_float) mean<float>;
%template(median_float) median<float>;
%template(variance_float) variance<float>;
%template(standardDeviation_float) standardDeviation<float>;
%template(norm1_float) norm1<float>;
%template(norm2_float) norm2<float>;
%template(normP_float) normP<float>;
%template(magnitude_float) magnitude<float>;
%template(multiplyInPlace_float) multiplyInPlace<float>;
%template(divideInPlace_float) divideInPlace<float>;
%template(addInPlace_float) addInPlace<float>;
%template(subtractInPlace_float) subtractInPlace<float>;
%template(absInPlace_float) absInPlace<float>;
%template(squareInPlace_float) squareInPlace<float>;
%template(squareRootInPlace_float) squareRootInPlace<float>;
%template(sort_float) sort<float>;
%template(reverse_float) reverse<float>;
%template(multiply_float) multiply<float>;
%template(divide_float) divide<float>;
%template(add_float) add<float>;
%template(subtract_float) subtract<float>;
%template(abs_float) abs<float>;
%template(square_float) square<float>;
%template(squareRoot_float) squareRoot<float>;
%template(scale_float) scale<float>;
%template(difference_float) difference<float>;
%template(zeros_float) zeros<float>;
%template(ones_float) ones<float>;
%template(range_float) range<float>;
%template(dotProduct_float) dotProduct<float>;
%template(euclideanDistance_float) euclideanDistance<float>;
%template(cosineSimilarity_float) cosineSimilarity<float>;
%template(cosineDistance_float) cosineDistance<float>;
