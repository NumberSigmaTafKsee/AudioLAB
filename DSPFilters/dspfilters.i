%module dspfilters
%{
#include "IIRDspFilters.hpp"
#include "IIRBesselFilterProcessor.hpp"
#include "IIRButterworthFilterProcessor.hpp"
#include "IIRChebyshevFilterProcessors.hpp"
#include "IIRChebyshev2FilterProcessors.hpp"
#include "IIREllipticalFilterProcessor.hpp"
#include "IIROptimalLFilterProcessor.hpp"
#include "DspButterworthLowPassQ.hpp"
using namespace Filters;
%}

%include "std_vector.i"

%include "IIRDspFilters.hpp"

%include "DspBesselLowPass.hpp"
%include "DspBesselHighPass.hpp"
%include "DspBesselBandPass.hpp"
%include "DspBesselBandStop.hpp"

%include "DspButterworthLowPass.hpp"
%include "DspButterworthHighPass.hpp"
%include "DspButterworthBandPass.hpp"
%include "DspButterworthBandStop.hpp"
%include "DspButterworthLowShelf.hpp"
%include "DspButterworthHighShelf.hpp"
%include "DspButterworthBandShelf.hpp"

%include "DspChebyshev1LowPass.hpp"
%include "DspChebyshev1HighPass.hpp"
%include "DspChebyshev1BandPass.hpp"
%include "DspChebyshev1BandStop.hpp"
%include "DspChebyshev1BandShelf.hpp"
%include "DspChebyshev1LowShelf.hpp"
%include "DspChebyshev1HighShelf.hpp"

%include "DspChebyshev2LowPass.hpp"
%include "DspChebyshev2HighPass.hpp"
%include "DspChebyshev2BandPass.hpp"
%include "DspChebyshev2BandStop.hpp"
%include "DspChebyshev2BandShelf.hpp"
%include "DspChebyshev2LowShelf.hpp"
%include "DspChebyshev2HighShelf.hpp"

%include "DspEllipticLowPass.hpp"
%include "DspEllipticHighPass.hpp"
%include "DspEllipticBandPass.hpp"
%include "DspEllipticBandStop.hpp"

%include "DspLegendreLowPass.hpp"
%include "DspLegendreHighPass.hpp"
%include "DspLegendreBandPass.hpp"
%include "DspLegendreBandStop.hpp"

%include "DspButterworthLowPassQ.hpp"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;