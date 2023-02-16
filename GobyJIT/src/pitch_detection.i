%module pitch_detection
%{
#include "pitch_detection.h"
%}
%include "std_vector.i"
%include "std_string.i"

typedef unsigned char uint8_t;
typedef char int8_t;
typedef unsigned short uint16_t;
typedef short int16_t;
typedef unsigned int uint32_t;
typedef int int32_t;
typedef unsigned long long int uint64_t;
typedef long long int int64_t;

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
//%template(int8_vector) std::vector<int8_t>;
//%template(uint8_vector) std::vector<uint8_t>;
%template(int16_vector) std::vector<int16_t>;
%template(uint16_vector) std::vector<uint16_t>;
%template(uint32_vector) std::vector<uint32_t>;
%template(int32_vector) std::vector<int32_t>;
///%template(uint64_vector) std::vector<uint64_t>;
///%template(int64_vector) std::vector<int64_t>;
%template(string_vector) std::vector<std::string>;

%include "pitch_detection.h"

//%template(bin_pitchesf) detail::bin_pitches<float>;


%template(yinf) pitch::yin<float>;
%template(mpmf) pitch::mpm<float>;
%template(swipef) pitch::swipe<float>;
%template(pyinf) pitch::pyin<float>;
%template(pmpmf) pitch::pmpm<float>;

%template(FloatBase) pitch_alloc::BaseAlloc<float>;
%template(FloatYin) pitch_alloc::Yin<float>;
%template(FloatMpm) pitch_alloc::Mpm<float>;

%template(parabolic_interpolationf) util::parabolic_interpolation<float>;
%template(acorr_rf) util::acorr_r<float>;
%template(pitch_from_hmmf) util::pitch_from_hmm<float>;
