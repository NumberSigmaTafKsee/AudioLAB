// going to specialize everything
#pragma once

#include <complex>
#include <algorithm>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <chrono>
#include <cmath>
#include <cassert>

#include <ippcore.h>
#include <ipps.h>

#include "carlo_samples.hpp"

#include "carlo_ipp.hpp"
#include "carlo_ipparray.hpp"
#include "carlo_random.hpp"
#include "carlo_autocorr.hpp"
#include "carlo_xcorr.hpp"
#include "carlo_convolution.hpp"
#include "carlo_dct.hpp"
#include "carlo_dft.hpp"
#include "carlo_fft.hpp"
#include "carlo_firlms.hpp"
#include "carlo_firmr.hpp"
#include "carlo_firsr.hpp"
#include "carlo_hilbert.hpp"
#include "carlo_iir.hpp"
#include "carlo_resample.hpp"