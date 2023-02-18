%module audiodsp_blit_oscillators
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
#include "Analog/VABlitTriangleOscillator.hpp"
%}

%include "SoundObject.hpp"
%include "Analog/VABlitSawOscillator.hpp"
%include "Analog/VABlitSquareOscillator.hpp"
%include "Analog/VABlitTriangleOscillator.hpp"
