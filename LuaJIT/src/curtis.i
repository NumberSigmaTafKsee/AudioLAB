%module curtis
%{
#include "TCurtisVCF.hpp"
using namespace SoundWave;
%}

%include "TCurtisVCF.hpp"

%template (float_curtisvcf) SoundWave::TCurtisVCF<float>;