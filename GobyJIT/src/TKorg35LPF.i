%module TKorg35LPF
%{
#include "TKorg35LPF.hpp"
using namespace SoundWave;
%}

%include "TKorg35LPF.hpp"

%template (FloatTKorg35LPF) SoundWave::TKorg35LPF<float>;
%template (DoubleTKorg35LPF) SoundWave::TKorg35LPF<double>;