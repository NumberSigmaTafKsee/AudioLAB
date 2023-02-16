#pragma once

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>

#include "StdNoise.hpp"
// this will only work on Intel
#include "Undenormal.hpp"

extern float sampleRate;
extern float invSampleRate;
extern Default noise;

#include "kfr2_audio.hpp"
#include "kfr2_combfilters.hpp"
#include "kfr2_delaylines.hpp"
#include "kfr2_noise.hpp"
#include "kfr2_utils.hpp"

#include "kfr2_biquads.hpp"
#include "kfr2_rbj.hpp"
#include "kfr2_zolzer.hpp"
#include "kfr2_iir.hpp"
#include "kfr2_fir.hpp"

#include "kfr2_dft.hpp"
#include "kfr2_convolution.hpp"
#include "kfr2_functions.hpp"
#include "Kfr2_resample.hpp"
#include "kfr2_windows.hpp"
