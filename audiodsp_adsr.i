%module audio_adsr
%{
#define DSPFLOATDOUBLE
#include "SoundObject.hpp"    
#include "audio_adsr.hpp"
%}

#define DSPFLOATDOUBLE
%include "SoundObject.hpp"    
%include "audio_adsr.hpp"