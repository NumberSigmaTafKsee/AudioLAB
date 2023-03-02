%module Analog
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>


#include "Analog/VABandLimitedOscillators.hpp"
#include "Analog/VAPolyBLEPOscillator.hpp"
#include "Analog/VAMinBlepOscillators.hpp"

#include "Analog/VADiodeLadderFilter2.hpp"
#include "Analog/VAAnalogSVF.hpp"
#include "Analog/VAGenSVF.hpp"
#include "Analog/VAStateVariableFilter.hpp"
#include "Analog/VAMoogLadders.hpp"
#include "Analog/VASVFFilter.hpp"

#include "Analog/VAOBXDFilter.hpp"
#include "Analog/VAMorphableFilter.hpp"
#include "Analog/VAMS20Filter.hpp"
#include "Analog/VARCFilter.hpp"
#include "Analog/VASstFilters.hpp"
#include "Analog/VASstWaveshaper.hpp"
#include "Analog/VASlewLimiter.hpp"
#include "Analog/VATwoPoleEnvelopes.hpp"
#include "Analog/VAVCS3DiodeFilter.hpp"
#include "Analog/VAVCS3Filter.hpp"
#include "Analog/VAXodFilters.hpp"

/*
#include "Analog/VADinkyFilter.hpp"
#include "Analog/VAHybridCurtisVCF.hpp"
#include "Analog/VAKorg35HPFFilter.hpp"
#include "Analog/VAKorg35HPFilter.cpp
#include "Analog/VAKorg35LPFFilter.cpp
#include "Analog/VAKorg35LPFFilter.hpp"
#include "Analog/VAOberheimFilter.cpp
#include "Analog/VAOberheimFilter.hpp"
#include "Analog/VAVecSVF.hpp"
#include "Analog/VAVectorSVF.hpp"
*/
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

//%include "Analog/VATwoPoleEnvelopes.hpp"

%include "Analog/VABlitSawOscillator.hpp"
%include "Analog/VABlitSquareOscillator.hpp"
%include "Analog/VABlitTriangleOscillator.hpp"
%include "Analog/VablitOscillators.hpp"
%include "Analog/VADPWSawOscillator.hpp"    
%include "Analog/VADPWPulseOscillator.hpp"
%include "Analog/VADPWTriangleOscillator.hpp"

%ignore Analog::Oscillators::PolyBLEPOsc::blep;
%ignore Analog::Oscillators::PolyBLEPOsc::blamp;
%include "Analog/VAPolyBLEPOscillator.hpp"
%include "FX/minBLEP.hpp"


%rename Analog::Filters::StateVariableFilter::SVFLowpass            VASVFLowPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFBandpass           VASVFBandPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFHighpass           VASVFHighPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFUnitGainBandpass   VASVFUnitGainBandpassFilter;
%rename Analog::Filters::StateVariableFilter::SVFBandShelving       VASVFBandShelvingFilter;
%rename Analog::Filters::StateVariableFilter::SVFNotch              VASVFNotchFilter;
%rename Analog::Filters::StateVariableFilter::SVFPeak               VASVFPeakFilter;
%rename Analog::Filters::StateVariableFilter::SVFAllpass            VASVFAllpassFilter;
%ignore Analog::Filters::StateVariableFilter::resonanceToQ;
%ignore Analog::Filters::StateVariableFilter::analogSaturate;
%include "Analog/VAStateVariableFilter.hpp"

%rename Analog::Filters::AnalogSVF              VAAnalogSVF;
%include "Analog/VAAnalogSVF.hpp"

%include "Analog/VAMoogLadders.hpp"


%include "Analog/VASlewLimiter.hpp"
%include "Analog/VARCFilter.hpp"


%include "Analog/VASstFilters.hpp"
%include "Analog/VASstWaveshaper.hpp"

%include "Analog/VAMS20Filter.hpp"
%include "Analog/VAOBXDFilter.hpp"

%include "Analog/VAVCS3DiodeFilter.hpp"
%include "Analog/VAVCS3Filter.hpp"

%include "Analog/VAXodFilters.hpp"

//%ignore FirstOrderFilter;
//%include "Analog/VAMorphableFilter.hpp"
//%include "Analog/VirtualAnalogDiodeLadderFilter2.hpp"

/*
// it is just the musicdsp moog
%include "Analog/VAHybridCurtisVCF.hpp"
%include "Analog/VAAnalogFilters.hpp"
%include "Analog/VADiodeLadderFilter.hpp"
%include "Analog/VADiode.hpp"
%include "Analog/VADiodeClipper.hpp"
%include "Analog/VADinkyFilter.hpp"
%include "Analog/VAKorg35HPFFilter.hpp"
%include "Analog/VAKorg35LPFFilter.hpp"
%include "Analog/VAOberheimFilter.hpp"
*/


%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;
