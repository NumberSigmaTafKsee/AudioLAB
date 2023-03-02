%module audiodsp_analog
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
#include "Analog/VADPWOscillators.hpp"

#include "Analog/VAImprovedMoogFilter.hpp"
#include "Analog/VAKrajeskiMoogFilter.hpp"
#include "Analog/VALadderFilter.hpp"
#include "Analog/VALadderFilter2.hpp"
#include "Analog/VAMicroTrackerMoogFilter.hpp"
#include "Analog/VAMoogCatFilter.hpp"
#include "Analog/VAMoogFilter.hpp"
#include "Analog/VAMoogFilter1.hpp"
#include "Analog/VAMoogFilter2.hpp"
#include "Analog/VAMoogFilter3.hpp"
#include "Analog/VAMoogFilter4.hpp"
#include "Analog/VAMoogFilterI.hpp"
#include "Analog/VAMoogFilterII.hpp"
#include "Analog/VAMoogFilters.hpp"
#include "Analog/VAMoogLadderFilter.hpp"
#include "Analog/VAMoogLadderFilters.hpp"
#include "Analog/VAMoogLadders.hpp"
#include "Analog/VAMoogLikeFilter.hpp"
#include "Analog/VAMoogNonLinearFilter.hpp"
#include "Analog/VAMoogNonLinearFilter2.hpp"
#include "Analog/VAMoogRKLadderFilter.hpp"
#include "Analog/VAMoogVCFFilter.hpp"
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

#include "Analog/VADinkyFilter.hpp"
#include "Analog/VADiodeLadderFilter2.hpp"
#include "Analog/VAGenSVF.hpp"
#include "Analog/VAMorphableFilter.hpp"
#include "Analog/VAMS20Filter.hpp"
#include "Analog/VAOBXDFilter.hpp"
#include "Analog/VARCFilter.hpp"
#include "Analog/VARKLadderFilter.hpp"
#include "Analog/VASlewLimiter.hpp"
#include "Analog/VASstFilters.hpp"
#include "Analog/VAVCS3DiodeFilter.hpp"
#include "Analog/VAVCS3Filter.hpp"
#include "Analog/VAXodFilters.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"

%include "SoundObject.hpp"

%include "Analog/VABlitSawOscillator.hpp"
%include "Analog/VABlitSquareOscillator.hpp"
%include "Analog/VABlitTriangleOscillator.hpp"
%include "Analog/VablitOscillators.hpp"
%include "Analog/VADPWSawOscillator.hpp"    
%include "Analog/VADPWPulseOscillator.hpp"
%include "Analog/VADPWTriangleOscillator.hpp"
%include "Analog/VAMinBLEP.hpp"
%ignore Analog::Oscillators::PolyBLEPOsc::blep;
%ignore Analog::Oscillators::PolyBLEPOsc::blamp;
%include "Analog/VAPolyBLEPOscillator.hpp"


// SVF
%include "Analog/VAAnalogSVF.hpp"
%include "Analog/VAGenSVF.hpp"
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
%rename Analog::Filters::SVF::AnalogSVF VASVF;
%include "Analog/VASVF.hpp"
%include "Analog/VASVFChamberlinFilter.hpp"
%rename Analog::Filters::SVF::StateVariableFilter VASVFFilter;
%include "Analog/VASVFFilter.hpp"
%include "Analog/VASVFSmoother.hpp"
%include "Analog/VASVSmoothFilter.hpp"



// Moog
%include "Analog/VAImprovedMoogFilter.hpp"
%include "Analog/VAKrajeskiMoogFilter.hpp"
%ignore Analog::Filters::LadderFilter2::TEMP;
%ignore Analog::Filters::LadderFilter2::THERMAL_VOLT;
%ignore Analog::Filters::LadderFilter2::OVER_TWO_THERMAL_VOLT;
%ignore Analog::Filters::LadderFilter2::NUMBER_OF_FILTERS;
%include "Analog/VALadderFilter.hpp"
%include "Analog/VALadderFilter2.hpp"
%include "Analog/VAMicroTrackerMoogFilter.hpp"
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
%ignore Analog::Moog::clip;
%rename Analog::Moog::MoogLike      VAMoogLikeFilter;
%rename Analog::Moog::ImprovedMoog  VAImprovedMoogFilter;
%rename Analog::Moog::KrajeskiMoog  VAKrajeskiMoogFilter;
%rename Analog::Moog::MicrotrackerMoog VAMicrotrackerMoogFilter;
%rename Analog::Moog::MoogFilter2   VAMoogLadderFilter2;
%rename Analog::Moog::MoogVCF       VAMoogLadderVCF;
%rename Analog::Moog::MoogFilter1   VAMoogLadderFilter1;
%include "Analog/VAMoogLadderFiltersBase.hpp"
%include "Analog/VAMoogLadderFilters.hpp"
%include "Analog/VAMoogLikeFilter.hpp"
%include "Analog/VAMoogNonLinearFilter.hpp"
%rename Analog::Filters::Moog::NonLinear2::MoogFilter VANonLinearMoogFilter;
%include "Analog/VAMoogNonLinearFilter2.hpp"
%include "Analog/VAMoogRKLadderFilter.hpp"
%rename Analog::Filters::Moog::MoogVCF::MoogVCF VAMoogVoltageControlledFilter;
%include "Analog/VAMoogVCFFilter.hpp"
%ignore Analog::Filters::Moog::StilsonMoog::gaintable;
%rename Analog::Filters::Moog::StilsonMoog::StilsonMoog VAStilsonMoogFilter;
%include "Analog/VAStilsonMoogFilter.hpp"
%rename Analog::Filters::Moog::StilsonMoogFilter2::StilsonMoog VAStilsonMoogFilter2;
%include "Analog/VAStilsonMoogFilter2.hpp"


// other filters
%include "Analog/VADinkyFilter.hpp"
%include "Analog/VADiodeLadderFilter2.hpp"
%ignore FirstOrderFilter;
%include "Analog/VAMorphableFilter.hpp"
%include "Analog/VAMS20Filter.hpp"
%include "Analog/VAOBXDFilter.hpp"
%include "Analog/VARCFilter.hpp"
%ignore Analog::Filters::RKLadderFilter::clip;
%ignore Analog::Filters::RKLadderFilter::crossfade;
%ignore Analog::Filters::RKLadderFilter::stepRK4;
%rename Analog::Filters::RKLadderFilter::LadderFilter RKLadderFilter;
%include "Analog/VARKLadderFilter.hpp"
%include "Analog/VASlewLimiter.hpp"
%include "Analog/VASstFilters.hpp"
%include "Analog/VATwoPoleEnvelopes.hpp"
%include "Analog/VAVCS3DiodeFilter.hpp"
%include "Analog/VAVCS3Filter.hpp"
%include "Analog/VAXodFilters.hpp"

%rename Analog::Filters::VirtualAnalogDiodeLadderFilter VADiodeLadderFilter;
%include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
%include "Analog/VirtualAnalogStateVariableFilter.hpp"
%include "Analog/VoltageControlledFilter.hpp"
%include "Analog/VoltageControlledOscillator.hpp"
*/
%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;
