%module kfrdsp
%{
#include "KfrDSP/KfrDsp.hpp"
#include "samples/kfr_sample.hpp"
#include "samples/kfr_sample_dsp.hpp"
#include "IIRFilters.hpp"
using namespace KfrDSP1;
using namespace Filters;


%}

typedef float SampleType;

%include "stdint.i"
%include "std_complex.i"
%include "std_vector.i"

%include "samples/kfr_sample.hpp"
%include "samples/kfr_sample_dsp.hpp"

//%include "KfrDSP/KfrAudio.hpp"
%include "KfrDSP/KfrCombFilters.hpp"
%include "KfrDSP/KfrDelayLine.hpp"
%include "KfrDSP/KfrNoise.hpp"
%include "KfrDSP/KfrUtils.hpp"
%include "KfrDSP/KfrRBJ.hpp"
%include "KfrDSP/KfrZolzer.hpp"
%include "KfrDSP/KfrBiquads.hpp"
%include "KfrDSP/KfrIIR.hpp"
%include "KfrDSP/KfrWindows.hpp"
%include "KfrDSP/KfrFIR.hpp"
%include "KfrDSP/KfrConvolution.hpp"
%include "KfrDSP/KfrDFT.hpp"
%include "KfrDSP/KfrFunctions.hpp"
%include "KfrDSP/KfrResample.hpp"

%include "IIRFilters.hpp"

%template(SampleVector)  kfr::univector<SampleType>;
%template(SampleMatrix)  kfr::univector2d<SampleType>;
%template(ComplexVector) kfr::univector<kfr::complex<SampleType>>;
%template(ComplexMatrix) kfr::univector2d<kfr::complex<SampleType>>;

%template(ConvolutionFilter) KfrDSP1::ConvolutionFilter<SampleType>;
%template(StereoConvolutionFilter) KfrDSP1::StereoConvolutionFilter<SampleType>;

%inline %{
float sampleRate=44100.0;
float invSampleRate=1.0/sampleRate;
Std::RandomMersenne noise;
%}