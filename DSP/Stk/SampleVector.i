%module sv
%{
#include "SampleVector.h"

using namespace Samples;
%}

%include "SampleVector.h"

%extend Samples::SampleVector {
    T __getitem__(size_t i) { return (*$self)[i-1]; }
    void __setitem__(size_t i, T v) { (*$self)[i-1] = v;}
}

%template(SampleVector) Samples::SampleVector<float>;