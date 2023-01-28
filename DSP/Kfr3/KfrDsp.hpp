#pragma once

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>

#include "Std/StdObject.h"
#include "Std/StdRandom.h"

//#include "audiosystem.h"
//#include "samples/kfr_sample.hpp"
//#include "samples/kfr_sample_dsp.hpp"


#include "Undenormal.hpp"
//#include "ADSR.h"
//#include "PolyBLEP.h"

//using namespace SoundAlchemy;

extern float sampleRate;
extern float invSampleRate;
extern Std::RandomMersenne noise;

#include "KfrAudio.hpp"
#include "KfrCombFilters.hpp"
#include "KfrDelayLines.hpp"
#include "KfrNoise.hpp"
#include "KfrUtils.hpp"

#include "KfrBiquads.hpp"
#include "KfrRBJ.hpp"
#include "KfrZolzer.hpp"
#include "KfrIIR.hpp"
#include "KfrFIR.hpp"

#include "KfrDFT.hpp"
#include "KfrConvolution.hpp"
#include "KfrFunctions.hpp"
#include "KfrResample.hpp"
#include "KfrWindows.hpp"
