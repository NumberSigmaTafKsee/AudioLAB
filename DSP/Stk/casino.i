%module casino
%{
#include "carlo_casino.hpp"
#include "carlo_casinodsp.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

%template(short_vector) std::vector<short>;
%template(int_vector) std::vector<int>;
%template(long_vector) std::vector<long>;
%template(llong_vector) std::vector<long long>;

%template(ushort_vector) std::vector<unsigned short>;
%template(uint_vector) std::vector<unsigned int>;
%template(ulong_vector) std::vector<unsigned long>;
%template(ullong_vector) std::vector<unsigned long long>;

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;


%include "carlo_casino.hpp"
%include "carlo_casinodsp.hpp"

%template(RandomUniformF32) Casino::RandomUniform<float>;
%template(IPPArrayF32) Casino::IPPArray<float>;
%template(AutoCorrF32) Casino::AutoCorr<float>;
%template(CrossCorrF32) Casino::CrossCorr<float>;
%template(ConvolutionF32) Casino::Convolver<float>;
%template(ConvolutionFilterF32) Casino::ConvolutionFilter<float>;
%template(CFFTF32) Casino::CFFT<float>;
%template(RFFTF32) Casino::RFFT<float>;
%template(CDFTF32) Casino::CDFT<float>;
%template(RDFTF32) Casino::RDFT<float>;
%template(DCTF32) Casino::DCT<float>;
%template(FIRMRF32) Casino::FIRMR<float>;
%template(FIRSRF32) Casino::FIRSR<float>;
%template(IIRF32) Casino::IIR<float>;
%template(IIRBiquadF32) Casino::IIRBiquad<float>;
