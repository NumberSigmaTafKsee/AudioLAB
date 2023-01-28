%module sigpack
%{
#define HAVE_FFTW
#include "base/base.h"
#include "window/window.h"
#include "filter/filter.h"
#include "resampling/resampling.h"
#include "spectrum/spectrum.h"
#include "timing/timing.h"
#include "gplot/gplot.h"
#include "parser/parser.h"
#include "fftw/fftw.h"
#include "image/image.h"
#include "kalman/kalman.h"
%}

%include "base/base.h"
%include "window/window.h"
%include "filter/filter.h"
%include "resampling/resampling.h"
%include "spectrum/spectrum.h"
%include "timing/timing.h"
%include "gplot/gplot.h"
%include "parser/parser.h"
%include "fftw/fftw.h"
%include "image/image.h"
%include "kalman/kalman.h"


