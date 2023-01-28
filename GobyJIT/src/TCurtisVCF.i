%module TCurtisVCF
%{
#include "TCurtisVCF.hpp"
%}

%include "TCurtisVCF.hpp"

%template( FloatCurtisVCF ) SoundWave::TCurtisVCF<float>;
%template( DoubleCurtisVCF ) SoundWave::TCurtisVCF<double>;