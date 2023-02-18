%module audiodsp_blit2_oscillators
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VablitOscillators.hpp"
%}

%include "SoundObject.hpp"
%include "Analog/VablitOscillators.hpp"
