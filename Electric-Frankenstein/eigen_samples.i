%module sv
%{
#include "SampleVector.h"

using namespace Casino::TinyEigen;
%}

%include "SampleVector.h"

%extend Casino::TinyEigen::SampleVector {
    T __getitem__(size_t i) { return (*$self)[i-1]; }
    void __setitem__(size_t i, T v) { (*$self)[i-1] = v;}
}

%template(SampleVector) Casino::TinyEigen::SampleVector<float>;
