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

#include "Analog/VAAnalogSVF.hpp"
#include "Analog/VABandLimitedOscillators.hpp"
#include "Analog/VABlitOscillators.hpp"
#include "Analog/VADinkyFilter.hpp"
#include "Analog/VADiode.hpp"
#include "Analog/VADiodeClipper.hpp"
//#include "Analog/VADiodeLadderFilter.cpp
//#include "Analog/VADiodeLadderFilter.hpp"
#include "Analog/VADiodeLadderFilter2.hpp"
//#include "Analog/VADiodeSimulator.hpp"
#include "Analog/VADPWOscillators.hpp"
#include "Analog/VAGenSVF.hpp"
//#include "Analog/VAHybridCurtisVCF.hpp"
#include "Analog/VAImprovedMoogFilter.hpp"
//#include "Analog/VAKorg35HPFFilter.hpp"
//#include "Analog/VAKorg35HPFilter.cpp
//#include "Analog/VAKorg35LPFFilter.cpp
//#include "Analog/VAKorg35LPFFilter.hpp"
#include "Analog/VAKrajeskiMoogFilter.hpp"
#include "Analog/VALadderFilter.hpp"
#include "Analog/VALadderFilter2.hpp"
#include "Analog/VAMicroTrackerMoogFilter.hpp"
#include "Analog/VAMinBlepOscillators.hpp"
#include "Analog/VAMoogCatFilter.hpp"
#include "Analog/VAMoogFilter.hpp"
#include "Analog/VAMoogFilter1.hpp"
#include "Analog/VAMoogFilter2.hpp"
#include "Analog/VAMoogFilter3.hpp"
#include "Analog/VAMoogFilter4.hpp"
#include "Analog/VAMoogFilterI.hpp"
#include "Analog/VAMoogFilterII.hpp"
#include "Analog/VAMoogFilters.hpp"
//#include "Analog/VAMoogHalfLadderFilter.cpp
//#include "Analog/VAMoogHalfLadderFilter.hpp"
//#include "Analog/VAMoogLadder.hpp"
//#include "Analog/VAMoogLadderFilter.cpp
#include "Analog/VAMoogLadderFilter.hpp"
#include "Analog/VAMoogLadderFilters.hpp"
#include "Analog/VAMoogLadders.hpp"
#include "Analog/VAMoogLikeFilter.hpp"
#include "Analog/VAMoogNonLinearFilter.hpp"
#include "Analog/VAMoogNonLinearFilter2.hpp"
#include "Analog/VAMoogRKLadderFilter.hpp"
#include "Analog/VAMoogVCFFilter.hpp"
#include "Analog/VAMorphableFilter.hpp"
#include "Analog/VAMS20Filter.hpp"
//#include "Analog/VAOberheimFilter.cpp
//#include "Analog/VAOberheimFilter.hpp"
#include "Analog/VAOBXDFilter.hpp"
#include "Analog/VAOscillators.hpp"
#include "Analog/VAPolyBLEPOscillator.hpp"
#include "Analog/VAPolyBlepOscillators.hpp"
//#include "Analog/VAPolygonalOscillator.hpp"
#include "Analog/VARCFilter.hpp"
#include "Analog/VARKLadderFilter.hpp"
#include "Analog/VASlewLimiter.hpp"
#include "Analog/VASstFilters.hpp"
//#include "Analog/VAStateVariableCombFilter.hpp"
#include "Analog/VAStateVariableFilter.hpp"
#include "Analog/VAStateVariableFilter1.hpp"
#include "Analog/VAStateVariableFilter2.hpp"
#include "Analog/VAStateVariableFilters.hpp"
#include "Analog/VAStilsonMoogFilter.hpp"
#include "Analog/VAStilsonMoogFilter2.hpp"
#include "Analog/VASVF.hpp"
#include "Analog/VASVFChamberlinFilter.hpp"
#include "Analog/VASVFFilter.hpp"
#include "Analog/VASVFSmoother.hpp"
#include "Analog/VASVSmoothFilter.hpp"
#include "Analog/VASVStateVariableFilter.hpp"
#include "Analog/VATwoPoleEnvelopes.hpp"
#include "Analog/VAVCS3DiodeFilter.hpp"
#include "Analog/VAVCS3Filter.hpp"
#include "Analog/VAVecSVF.hpp"
#include "Analog/VAVectorSVF.hpp"
#include "Analog/VAVoltageControlledFilter.hpp"
//#include "Analog/VAVoltageControlledOscillator.hpp"
/*
#include "Analog/VAWDFCompressor.hpp"
#include "Analog/VAWDFDiodeClipper.hpp"
#include "Analog/VAWDFPassiveLPF.hpp"
#include "Analog/VAWDFSallenKey.hpp"
*/
#include "Analog/VAXodFilters.hpp"
#include "Analog/VCA.hpp"
#include "Analog/VCAProcessor.hpp"
#include "Analog/VCF.hpp"
#include "Analog/VCFProcessor.hpp"
#include "Analog/VCO.hpp"
#include "Analog/VCOProcessor.hpp"
#include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
#include "Analog/VirtualAnalogStateVariableFilter.hpp"
#include "Analog/VoltageControlledFilter.hpp"
//#include "Analog/VoltageControlledOscillator.hpp"
%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

//%include "Analog/VAAnalogFilters.hpp"
%rename Analog::Filters::AnalogSVF              VAAnalogSVF;
%include "Analog/VAAnalogSVF.hpp"
%include "Analog/VABandLimitedOscillators.hpp"
%rename Analog::Oscillators::Blit::BlitSaw      VABlitSawOscillator;
%rename Analog::Oscillators::Blit::BlitSquare   VABlitSquareOscillator;
%include "Analog/VABlitOscillators.hpp"
%include "Analog/VADinkyFilter.hpp"
%include "Analog/VADiode.hpp"
%include "Analog/VADiodeClipper.hpp"
//%include "Analog/VADiodeLadderFilter.hpp"
%include "Analog/VADiodeLadderFilter2.hpp"
%rename Analog::Oscillators::DPW::DPWSaw        VADPWSawOscillator;
%rename Analog::Oscillators::DPW::DPWSquare     VADPWSquareOscillator;
%rename Analog::Oscillators::DPW::DPWTriangle   VADPWTriangleOscillator;
%rename Analog::Oscillators::DPW::DPWPulse      VADPWPulseOscillator;
%include "Analog/VADPWOscillators.hpp"
%include "Analog/VAGenSVF.hpp"
// it is just the musicdsp moog
//%include "Analog/VAHybridCurtisVCF.hpp"
%include "Analog/VAImprovedMoogFilter.hpp"
//%include "Analog/VAKorg35HPFFilter.hpp"
//%include "Analog/VAKorg35LPFFilter.hpp"
%include "Analog/VAKrajeskiMoogFilter.hpp"
%ignore Analog::Filters::LadderFilter2::TEMP;
%ignore Analog::Filters::LadderFilter2::THERMAL_VOLT;
%ignore Analog::Filters::LadderFilter2::OVER_TWO_THERMAL_VOLT;
%ignore Analog::Filters::LadderFilter2::NUMBER_OF_FILTERS;
%include "Analog/VALadderFilter.hpp"
%include "Analog/VALadderFilter2.hpp"
%include "Analog/VAMicroTrackerMoogFilter.hpp"
%include "Analog/VAMinBlepOscillators.hpp"
%include "Analog/VAMoogCatFilter.hpp"
%include "Analog/VAMoogFilter.hpp"
%include "Analog/VAMoogFilter1.hpp"
%include "Analog/VAMoogFilter2.hpp"
%include "Analog/VAMoogFilter3.hpp"
%rename Analog::MoogFilters::MoogFilter4::MoogFilter VAMoogFilter4;
%include "Analog/VAMoogFilter4.hpp"
%include "Analog/VAMoogFilterI.hpp"
%include "Analog/VAMoogFilterII.hpp"
%include "Analog/VAMoogFilters.hpp"
//%include "Analog/VAMoogHalfLadderFilter.hpp"
//%include "Analog/VAMoogLadder.hpp"
//%rename Analog::Filters::MoogLadder::MoogLadder MoogLadder1;
//%include "Analog/VAMoogLadderFilter.hpp"
%ignore Analog::Moog::clip;
%rename Analog::Moog::MoogLike      VAMoogLikeFilter;
%rename Analog::Moog::ImprovedMoog  VAImprovedMoogFilter;
%rename Analog::Moog::KrajeskiMoog  VAKrajeskiMoogFilter;
%rename Analog::Moog::MicrotrackerMoog VAMicrotrackerMoogFilter;
%rename Analog::Moog::MoogFilter2   VAMoogLadderFilter2;
%rename Analog::Moog::MoogVCF       VAMoogLadderVCF;
%rename Analog::Moog::MoogFilter1   VAMoogLadderFilter1;

%include "Analog/VAMoogLadderFilters.hpp"
%include "Analog/VAMoogLikeFilter.hpp"
%include "Analog/VAMoogNonLinearFilter.hpp"

%rename Analog::Filters::Moog::NonLinear2::MoogFilter VANonLinearMoogFilter;
%include "Analog/VAMoogNonLinearFilter2.hpp"
%include "Analog/VAMoogRKLadderFilter.hpp"

%rename Analog::Filters::Moog::MoogVCF::MoogVCF VAMoogVoltageControlledFilter;
%include "Analog/VAMoogVCFFilter.hpp"

%ignore FirstOrderFilter;
%include "Analog/VAMorphableFilter.hpp"
%include "Analog/VAMS20Filter.hpp"
//%include "Analog/VAOberheimFilter.hpp"
%include "Analog/VAOBXDFilter.hpp"

%rename Oscillators::BlitSaw        VABlitSawOsc;
%rename Oscillators::BlitSquare     VABlitSquareOsc;
%rename Oscillators::BlitTriangle   VABlitTriangleOsc;
%rename Oscillators::BlitDSF        VABlitDSFOsc;
%rename Oscillators::blitSaw        VABlitSawOsc2;
%rename Oscillators::blitSquare     VABlitSquareOsc2;
%rename Oscillators::blitTriangle   VABlitTriangleOsc2;
%rename Oscillators::DPWSaw         VADPWSawOsc;
%rename Oscillators::DPWPulse       VADPWPulseOsc;
%rename Oscillators::DPWTriangle    VADPWTriangleOsc;;
%include "Analog/VAOscillators.hpp"

%ignore Analog::Oscillators::PolyBLEPOsc::blep;
%ignore Analog::Oscillators::PolyBLEPOsc::blamp;
%include "Analog/VAPolyBLEPOscillator.hpp"
%include "Analog/VARCFilter.hpp"
%ignore Analog::Filters::RKLadderFilter::clip;
%ignore Analog::Filters::RKLadderFilter::crossfade;
%ignore Analog::Filters::RKLadderFilter::stepRK4;
%rename Analog::Filters::RKLadderFilter::LadderFilter RKLadderFilter;
%include "Analog/VARKLadderFilter.hpp"
%include "Analog/VASlewLimiter.hpp"
%include "Analog/VASstFilters.hpp"

%rename Analog::Filters::StateVariableFilter::SVFLowpass            VASVFLowPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFBandpass           VASVFBandPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFHighpass           VASVFHighPassFilter;
%rename Analog::Filters::StateVariableFilter::SVFUnitGainBandpass   VASVFUnitGainBandpassFilter;
%rename Analog::Filters::StateVariableFilter::SVFBandShelving       VASVFBandShelvingFilter;
%rename Analog::Filters::StateVariableFilter::SVFNotch              VASVFNotchFilter;
%rename Analog::Filters::StateVariableFilter::SVFPeak               VASVFPeakFilter;
%rename Analog::Filters::StateVariableFilter::SVFAllpass            VASVFAllpassFilter;
%ignore Analog::Filters::StateVariableFilter::resonanceToQ;

%include "Analog/VAStateVariableFilter.hpp"
%include "Analog/VAStateVariableFilter1.hpp"

%rename Analog::Filters::StateVariableFilter2::StateVariableFilter VAStateVariableFilter2;
%include "Analog/VAStateVariableFilter2.hpp"
%ignore Analog::Filters::Moog::StilsonMoog::gaintable;
%rename Analog::Filters::Moog::StilsonMoog::StilsonMoog VAStilsonMoogFilter;
%include "Analog/VAStilsonMoogFilter.hpp"
%rename Analog::Filters::Moog::StilsonMoogFilter2::StilsonMoog VAStilsonMoogFilter2;
%include "Analog/VAStilsonMoogFilter2.hpp"

%rename Analog::Filters::SVF::AnalogSVF VASVF;
%include "Analog/VASVF.hpp"
%include "Analog/VASVFChamberlinFilter.hpp"
%rename Analog::Filters::SVF::StateVariableFilter VASVFFilter;
%include "Analog/VASVFFilter.hpp"
%include "Analog/VASVFSmoother.hpp"
%include "Analog/VASVSmoothFilter.hpp"
%include "Analog/VATwoPoleEnvelopes.hpp"
%include "Analog/VAVCS3DiodeFilter.hpp"
%include "Analog/VAVCS3Filter.hpp"
%include "Analog/VAVecSVF.hpp"
%include "Analog/VAVectorSVF.hpp"
%ignore Analog::Filters::VoltageControlledFilter::clip;
%rename Analog::Filters::VoltageControlledFilter::LadderFilter VAVoltageControlledFilter;
%include "Analog/VAVoltageControlledFilter.hpp"
%include "Analog/VAVoltageControlledOscillator.hpp"
%include "Analog/VAXodFilters.hpp"
%include "Analog/VCA.hpp"
%include "Analog/VCAProcessor.hpp"
%include "Analog/VCF.hpp"
%include "Analog/VCFProcessor.hpp"
%include "Analog/VCO.hpp"
%include "Analog/VCOProcessor.hpp"

//%include "Analog/VAPolyBlepOscillators.hpp"
//%include "Analog/VASVStateVariableFilter.hpp"
//%include "Analog/VAStateVariableFilters.hpp"
//%include "Analog/VAPolygonalOscillator.hpp"
//%include "Analog/VASVStateVariableFilter.hpp"
//%include "Analog/VAWDFCompressor.hpp"
//%include "Analog/VAWDFDiodeClipper.hpp"
//%include "Analog/VAWDFPassiveLPF.hpp"
//%include "Analog/VAWDFSallenKey.hpp"
//%rename Analog::Filters::VirtualAnalogDiodeLadderFilter VADiodeLadderFilter;
//%include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
//%include "Analog/VirtualAnalogStateVariableFilter.hpp"
//%include "Analog/VoltageControlledFilter.hpp"
//%include "Analog/VoltageControlledOscillator.hpp"

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
