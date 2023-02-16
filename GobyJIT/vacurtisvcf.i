%module vacurtisvcf
%{
typedef float DspFloatType;
#include "Analog/VAHybridCurtisVCF.hpp"
using namespace Analog::Filters::Curtis;
%}

typedef float DspFloatType;

%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"


%include "include/SoundObject.hpp"
%include "Analog/VAHybridCurtisVCF.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

%inline %{
    const int BufferSize = 256;
    Std::RandomMersenne noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}

