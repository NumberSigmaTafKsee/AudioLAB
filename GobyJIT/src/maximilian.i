%module maximilian
%{
#include <maximilian.h>
#include <libs/maxiAtoms.h>
#include <libs/maxiBark.h>
#include <libs/maxiClock.h>
#include <libs/maxiConvolve.h>
#include <libs/maxiFFT.h>
#include <libs/maxiGrains.h>
#include <libs/maxiMFCC.h>
#include <libs/maxiPolyBLEP.h>
#include <libs/maxiReverb.h>
#include <libs/maxiSynths.h>

%}

%rename maxiEnv::adsr(double,int) adsr_double_int;
%include <maximilian.h>

%include <libs/maxiAtoms.h>
%include <libs/maxiBark.h>
%include <libs/maxiClock.h>
%include <libs/maxiConvolve.h>
%include <libs/maxiFFT.h>
%include <libs/maxiGrains.h>
%include <libs/maxiMFCC.h>
%include <libs/maxiPolyBLEP.h>
%include <libs/maxiReverb.h>
%include <libs/maxiSynths.h>
