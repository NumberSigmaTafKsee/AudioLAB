%module vaanalogsvf
%{
#include "Analog/VAAnalogSVF.hpp"
%}

%include "std_vector.i"

typedef float DspFloatType;

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "include/SoundObject.hpp"
%include "Analog/VAAnalogSVF.hpp"

%inline %{
    const int BufferSize = 256;
    Std::RandomMersenne noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}