%module audiodsp_diodes
%{
typedef float DspFloatType;
#include "SoundObject.hpp"
#include "FX/Diode.hpp"
#include "FX/DiodeClipper.hpp"
#include "FX/DiodeSim.hpp"
#include "FX/DiodeSimulator.hpp"
#include "Analog/VASlewLimiter.hpp"
#include "Analog/VAWDFDiodeClipper.hpp"
#include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
#include "Analog/VADiode.hpp"
#include "Analog/VADiodeClipper.hpp"
#include "Analog/VADiodeLadderFilter2.hpp"
#include "Analog/VADiodeSimulator.hpp"
#include "Analog/VAVCS3DiodeFilter.hpp"
#include "Analog/VAVCS3Filter.hpp"

%}

%include "SoundObject.hpp"
%include "FX/Diode.hpp"
%include "FX/DiodeClipper.hpp"
%include "FX/DiodeSim.hpp"
%include "FX/DiodeSimulator.hpp"
