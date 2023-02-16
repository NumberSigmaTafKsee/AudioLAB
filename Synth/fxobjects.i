%module fxobjects
%{
#define DSPFLOATDOUBLE    
#include "SoundObject.hpp"
#include "FXObjects/FXObjects.hpp"
#include "FXObjects/FXOscillators.hpp"
#include "FXObjects/FXDelays.hpp"
#include "FXObjects/FXFilters.hpp"
#include "FXObjects/FXFirFilters.hpp"
#include "FXObjects/FXDistortion.hpp"
#include "FXObjects/FXDynamics.hpp"
#include "FXObjects/FXFFTSpectrum.hpp"
#include "FXObjects/FXPhaseVocoders.hpp"
#include "FXObjects/FXReverbs.hpp"
#include "FXObjects/FXSampleRate.hpp"
#include "FXObjects/FXWaveDigitalFilters.hpp"
%}

%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

#define DSPFLOATDOUBLE
%include "SoundObject.hpp"

%ignore makeWindow;
%include "FXObjects/FXObjects.hpp"
%include "FXObjects/FXOscillators.hpp"
%include "FXObjects/FXDelays.hpp"
%include "FXObjects/FXFilters.hpp"
%include "FXObjects/FXFirFilters.hpp"
%include "FXObjects/FXDistortion.hpp"
%include "FXObjects/FXDynamics.hpp"
%include "FXObjects/FXFFTSpectrum.hpp"
%include "FXObjects/FXObjects.hpp"
%include "FXObjects/FXPhaseVocoders.hpp"
%include "FXObjects/FXReverbs.hpp"
%include "FXObjects/FXSampleRate.hpp"
%include "FXObjects/FXWaveDigitalFilters.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

/*
%inline %{
    const int BufferSize = 256;
    Default noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}
*/