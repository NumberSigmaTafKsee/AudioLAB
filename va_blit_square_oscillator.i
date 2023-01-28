%module va_blit_square_oscillator
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VABlitSawOscillator.hpp"
#include "Analog/VABlitSquareOscillator.hpp"
%}

%include "SoundObject.hpp"
%include "Analog/VABlitSawOscillator.hpp"
%include "Analog/VABlitSquareOscillator.hpp"