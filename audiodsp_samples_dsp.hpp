#pragma once 

#include <complex>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cassert>

// CPU
//#include "audiodsp_samples.hpp"
#include <fftw3.h>

// GPU
// cuFFTW
// MatX
// Cuda Stl
// Cuda std::vector
// libcu++ = cuda stl
enum Direction {
    BACKWARD= FFTW_BACKWARD,
    FORWARD = FFTW_FORWARD,
};

#include "audiodsp_samples_dsp_windows.hpp"
#include "audiodsp_samples_dsp_c2cd_fft.hpp"
#include "audiodsp_samples_dsp_c2cd2d_fft.hpp"
#include "audiodsp_samples_dsp_c2cf_fft.hpp"
#include "audiodsp_samples_dsp_c2cf2d_fft.hpp"
#include "audiodsp_samples_dsp_c2rd_fft.hpp"
#include "audiodsp_samples_dsp_c2rd2d_fft.hpp"
#include "audiodsp_samples_dsp_c2rf_fft.hpp"
#include "audiodsp_samples_dsp_c2rf2d_fft.hpp"
#include "audiodsp_samples_dsp_r2cd_fft.hpp"
#include "audiodsp_samples_dsp_r2cd2d_fft.hpp"
#include "audiodsp_samples_dsp_r2cf_fft.hpp"
#include "audiodsp_samples_dsp_r2cf2d_fft.hpp"
#include "audiodsp_samples_dsp_real_fft.hpp"

