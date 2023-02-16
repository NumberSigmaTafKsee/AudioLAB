%module Filters
%{

#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>



#include "Filters/FIR.hpp"
#include "Filters/FIRFilter.hpp"
#include "Filters/FIRFilters.hpp"
#include "Filters/FIRSincFilter.hpp"
#include "Filters/IIRAnalog.hpp"
#include "Filters/IIRAnalogFilter.hpp"
#include "Filters/IIRBesselFilterProcessor.hpp"
#include "Filters/IIRBiquadFilters.hpp"
#include "Filters/IIRButterworth.hpp"
#include "Filters/IIRButterworthFilterProcessor.hpp"
#include "Filters/IIRCauerFilter.hpp"
#include "Filters/IIRChebyshev2FilterProcessors.hpp"
#include "Filters/IIRChebyshevFilterProcessors.hpp"
#include "Filters/IIRChebyshevFilters.hpp"
#include "Filters/IIRDCBlock.hpp"
#include "Filters/IIRDCFilter.hpp"
#include "Filters/IIREllipticalFilterProcessor.hpp"
#include "Filters/IIRFilters.hpp"
#include "Filters/IIROptimalLFilterProcessor.hpp"
#include "Filters/IIRRbjFilterProcessor.hpp"
#include "Filters/IIRRBJFilters.hpp"
#include "Filters/IIRZolzerFilter.hpp"
%}


%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"

///////////////////////////////////////////////
// Filters
///////////////////////////////////////////////

%include "Filters/FIR.hpp"
%include "Filters/FIRFilter.hpp"
%include "Filters/FIRSincFilter.hpp"

%include "Filters/IIRFilters.hpp"
%include "Filters/IIRAnalog.hpp"
%include "Filters/IIRAnalogFilter.hpp"
%include "Filters/IIRChebyshevFilters.hpp"
%include "Filters/IIRZolzerFilter.hpp"
%include "Filters/IIRButterworth.hpp"

/* These are DSPfilters now */
%rename Filters::IIR::Bessel::LowPassFilter     IIRBesselLowPass;
%rename Filters::IIR::Bessel::HighPassFilter    IIRBesselHighPass;
%rename Filters::IIR::Bessel::BandPassFilter    IIRBesselBandPass;
%rename Filters::IIR::Bessel::BandStopFilter    IIRBesselBandStop;
%include "Filters/IIRBesselFilterProcessor.hpp"

%rename Filters::IIR::Butterworth::LowPassFilter    IIRButterworthLowPass;
%rename Filters::IIR::Butterworth::HighPassFilter   IIRButterworthHighPass;
%rename Filters::IIR::Butterworth::BandPassFilter   IIRButterworthBandPass;
%rename Filters::IIR::Butterworth::BandStopFilter   IIRButterworthBandStop;
%rename Filters::IIR::Butterworth::PeakFilter       IIRButterworthPeak;
%rename Filters::IIR::Butterworth::AllPassFilter    IIRButterworthAllPass;
%rename Filters::IIR::Butterworth::HighShelfFilter  IIRButterworthHighShelf;
%rename Filters::IIR::Butterworth::LowShelfFilter   IIRButterworthLowShelf;
%rename Filters::IIR::Butterworth::BandShelfFilter  IIRButterworthBandShelf;
%include "Filters/IIRButterworthFilterProcessor.hpp"


%rename Filters::IIR::ChebyshevII::LowPassFilter    IIRChebyshevIILowPass;
%rename Filters::IIR::ChebyshevII::HighPassFilter   IIRChebyshevIIHighPass;
%rename Filters::IIR::ChebyshevII::BandPassFilter   IIRChebyshevIIBandPass;
%rename Filters::IIR::ChebyshevII::BandStopFilter   IIRChebyshevIIBandStop;
%rename Filters::IIR::ChebyshevII::PeakFilter       IIRChebyshevIIPeak;
%rename Filters::IIR::ChebyshevII::AllPassFilter    IIRChebyshevIIAllPass;
%rename Filters::IIR::ChebyshevII::HighShelfFilter  IIRChebyshevIIHighShelf;
%rename Filters::IIR::ChebyshevII::LowShelfFilter   IIRChebyshevIILowShelf;
%rename Filters::IIR::ChebyshevII::BandShelfFilter  IIRChebyshevIIBandShelf;
%include "Filters/IIRChebyshev2FilterProcessors.hpp"

%rename Filters::IIR::ChebyshevI::LowPassFilter     IIRChebyshevILowPass;
%rename Filters::IIR::ChebyshevI::HighPassFilter    IIRChebyshevIHighPass;
%rename Filters::IIR::ChebyshevI::BandPassFilter    IIRChebyshevIBandPass;
%rename Filters::IIR::ChebyshevI::BandStopFilter    IIRChebyshevIBandStop;
%rename Filters::IIR::ChebyshevI::PeakFilter        IIRChebyshevIPeak;
%rename Filters::IIR::ChebyshevI::AllPassFilter     IIRChebyshevIAllPass;
%rename Filters::IIR::ChebyshevI::HighShelfFilter   IIRChebyshevIHighShelf;
%rename Filters::IIR::ChebyshevI::LowShelfFilter    IIRChebyshevILowShelf;
%rename Filters::IIR::ChebyshevI::BandShelfFilter   IIRChebyshevIBandShelf;
%include "Filters/IIRChebyshevFilterProcessors.hpp"

%rename Filters::IIR::Legendre::LowPassFilter   IIRLegendreLowPass;
%rename Filters::IIR::Legendre::HighPassFilter  IIRLegendreHighPass;
%rename Filters::IIR::Legendre::BandPassFilter  IIRLegendreBandPass;
%rename Filters::IIR::Legendre::BandStopFilter  IIRLegendreBandStop;
%rename Filters::IIR::Legendre::PeakFilter      IIRLegendrePeak;
%rename Filters::IIR::Legendre::AllPaaFilter    IIRLegendreAllPass;
%rename Filters::IIR::Legendre::HighShelfFilter IIRLegendreHighShelf;
%rename Filters::IIR::Legendre::LowShelfFilter  IIRLegendreLowShelf;
%rename Filters::IIR::Legendre::BandShelfFilter IIRLegendreBandShelf;
%include "Filters/IIROptimalLFilterProcessor.hpp"

%rename Filters::IIR::RBJ::LowPassFilter    IIRRBJLowPass;
%rename Filters::IIR::RBJ::HighPassFilter   IIRRBJHighPass;
%rename Filters::IIR::RBJ::BandPassFilter   IIRRBJBandPass;
%rename Filters::IIR::RBJ::BandStopFilter   IIRRBJBandStop;
%rename Filters::IIR::RBJ::PeakFilter       IIRRBJPeak;
%rename Filters::IIR::RBJ::AllPaaFilter     IIRRBJAllPass;
%rename Filters::IIR::RBJ::HighShelfFilter  IIRRBJHighShelf;
%rename Filters::IIR::RBJ::LowShelfFilter   IIRRBJLowShelf;
%rename Filters::IIR::RBJ::BandShelfFilter  IIRRBJBandShelf;
%include "Filters/IIRRbjFilterProcessor.hpp"


%rename Filters::IIR::RBJFilters::RBJLowPassFilter  IIRRBJLowPassBiquad;
%rename Filters::IIR::RBJFilters::RBJHighPassFilter IIRRBJHighPassBiquad;
%rename Filters::IIR::RBJFilters::RBJBandPassFilter IIRRBJBandPassBiquad;
%rename Filters::IIR::RBJFilters::RBJBandStopFilter IIRRBJBandStopBiquad;
%rename Filters::IIR::RBJFilters::RBJPeakFilter     IIRRBJPeakingBiquad;
%rename Filters::IIR::RBJFilters::RBJAllPassFilter  IIRRBJAllPassBiquad;
%rename Filters::IIR::RBJFilters::RBJLowShelfFilter IIRRBJLowShelfBiquad;
%rename Filters::IIR::RBJFilters::RBJHighShelfFilter IIRRBJHighShelfBiquad;
%rename Filters::IIR::RBJFilters::RBJBandShelfFilter IIRRBJBandShelfBiquad;
%include "Filters/IIRRBJFilters.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;

%template(complex_float_vector) std::vector<std::complex<float>>;
%template(complex_double_vector) std::vector<std::complex<double>>;
