%module dr_libs 
%{
#include "dr_flac.h"
#include "dr_mp3.h"
#include "dr_wav.h"
%}
%include "stdint.i"
%include "dr_flac.h"
%include "dr_mp3.h"
%include "dr_wav.h"
