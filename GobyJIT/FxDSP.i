%module FxDSP
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "FX/FxDSP.hpp"
#include "FX/FxDSPBiquads.hpp"
#include "FX/FxDSPDecimator.hpp"
#include "FX/FxDSPDiodes.hpp"
#include "FX/FxDSPFFT.hpp"
#include "FX/FxDSPFIRFilter.hpp"
#include "FX/FxDSPLadderFilter.hpp"
#include "FX/FxDSPLinkwitzReillyFilter.hpp"
#include "FX/FxDSPMetering.hpp"
#include "FX/FxDSPMIDI.hpp"
#include "FX/FxDSPMultibandFilter.hpp"
#include "FX/FxDSPOnePoleFilter.hpp"
#include "FX/FxDSPOptoCoupler.hpp"
#include "FX/FxDSPPanLaw.hpp"
#include "FX/FxDSPPolySaturator.hpp"
#include "FX/FxDSPProcessors.hpp"
#include "FX/FxDSPRMSEstimator.hpp"
#include "FX/FxDSPSpectrumAnalyzer.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"


//%include "FX/FxDSP.hpp"
%include "FX/FxDSPBiquads.hpp"
%include "FX/FxDSPDecimator.hpp"
%include "FX/FxDSPDiodes.hpp"
%include "FX/FxDSPFFT.hpp"
%include "FX/FxDSPFIRFilter.hpp"
%include "FX/FxDSPLadderFilter.hpp"
%include "FX/FxDSPLinkwitzReillyFilter.hpp"
%include "FX/FxDSPMetering.hpp"
%include "FX/FxDSPMIDI.hpp"
%include "FX/FxDSPMultibandFilter.hpp"
%include "FX/FxDSPOnePoleFilter.hpp"
%include "FX/FxDSPOptoCoupler.hpp"
%include "FX/FxDSPPanLaw.hpp"
%include "FX/FxDSPPolySaturator.hpp"
%include "FX/FxDSPProcessors.hpp"
//%include "FX/FxDSPRMSEstimator.hpp"
%include "FX/FxDSPSpectrumAnalyzer.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;

/*
%inline %{
    const int BufferSize = 256;
    Std::RandomMersenne noise;
    DspFloatType sampleRate = 44100.0f;
    DspFloatType inverseSampleRate = 1 / 44100.0f;
    DspFloatType invSampleRate = 1 / 44100.0f;
%}
*/