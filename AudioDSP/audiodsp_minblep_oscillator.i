%module va_minblep_oscillator
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VAMinBLEP.hpp"
%}

%include "SoundObject.hpp"
%include "FX/OnePole.hpp"
%include "Analog/VAMinBLEP.hpp"