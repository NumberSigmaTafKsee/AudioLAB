%module vablitoscillator
%{
#include "Analog/VABlitOscillators.hpp"
using namespace Analog::Oscillators::Blit;

%}

%include "std_vector.i"

typedef float DspFloatType;

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%include "include/SoundObject.hpp"
%include "Analog/VABlitOscillators.hpp"

%inline %{
    const int BufferSize = 256;
    Std::RandomMersenne noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}