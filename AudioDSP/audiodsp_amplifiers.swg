%module audiodsp_amplifiers
%{
typedef float DspFloatType;
#include "SoundObject.hpp"
#include "audiodsp_amplifiers.hpp"
%}
typedef float DspFloatType;
%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

typedef float DspFloatType;
%include "SoundObject.hpp"
%include "FX/ClipFunctions.hpp"
%include "FX/Amplifier.hpp"
%include "FX/AmplifierFold.hpp"
%include "FX/Amplifiers.hpp"
%include "FX/AmplifiersUdo1.hpp"
%include "FX/ClipperCircuit.hpp"
%include "FX/ClippingFunctions.hpp"
%include "FX/ClipSerpentCurve.hpp"
%include "FX/ChebyDistortion.hpp"
%include "FX/DistortionCompressor.hpp"
%include "FX/DistortionFunctions.hpp"
%include "FX/WaveShaperATanSoftClip.hpp"
%include "FX/Waveshapers.hpp"
%include "FX/Waveshaping.hpp"
%include "FX/SstWaveshaper.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

