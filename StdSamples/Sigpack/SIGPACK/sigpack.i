%module sigpack
%{
#define HAVE_FFTW
#define arma_inline inline
#define SP_VERSION_MAJOR 1
#define SP_VERSION_MINOR 2
#define SP_VERSION_PATCH 7
#include <cmath>
#include <armadillo>
#include "base/base.h"
#include "window/window.h"
#include "filter/filter.h"
#include "resampling/resampling.h"
#include "fftw/fftw.h"
#include "spectrum/spectrum.h"
#include "timing/timing.h"
#include "gplot/gplot.h"
#include "parser/parser.h"
#include "image/image.h"
#include "kalman/kalman.h"


using namespace sp;
using namespace fftw;
%}

#define arma_inline

%include "base/base.h"
%include "fftw/fftw.h"
%include "window/window.h"
%include "filter/filter.h"
%include "resampling/resampling.h"
%include "spectrum/spectrum.h"
%include "timing/timing.h"
%include "gplot/gplot.h"
%include "parser/parser.h"
%include "image/image.h"
%include "kalman/kalman.h"


