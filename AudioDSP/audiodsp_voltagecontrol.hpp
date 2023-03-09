%module VoltageControl
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VCA.hpp"
#include "Analog/VCAProcessor.hpp"
#include "Analog/VCF.hpp"
#include "Analog/VCFProcessor.hpp"
#include "Analog/VCO.hpp"
#include "Analog/VCOProcessor.hpp"
#include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
#include "Analog/VirtualAnalogStateVariableFilter.hpp"
#include "Analog/VoltageControlledFilter.hpp"
//#include "Analog/VoltageControlledOscillator.hpp"

#include "Analog/VADiode.hpp"
#include "Analog/VADiodeClipper.hpp"
#include "Analog/VADiodeLadderFilter.cpp
#include "Analog/VADiodeLadderFilter.hpp"
#include "Analog/VADiodeSimulator.hpp"

#include "Analog/VAVoltageControlledFilter.hpp"
//#include "Analog/VAVoltageControlledOscillator.hpp"
/*
#include "Analog/VAWDFCompressor.hpp"
#include "Analog/VAWDFDiodeClipper.hpp"
#include "Analog/VAWDFPassiveLPF.hpp"
#include "Analog/VAWDFSallenKey.hpp"
*/

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"


//%ignore Analog::Filters::VoltageControlledFilter::clip;
//%rename Analog::Filters::VoltageControlledFilter::LadderFilter VAVoltageControlledFilter;
//%include "Analog/VAVoltageControlledFilter.hpp"
//%include "Analog/VAVoltageControlledOscillator.hpp"


/*
%include "Analog/VCA.hpp"
%include "Analog/VCAProcessor.hpp"
%include "Analog/VCF.hpp"
%include "Analog/VCFProcessor.hpp"
%include "Analog/VCO.hpp"
%include "Analog/VCOProcessor.hpp"
*/
//%rename Analog::Filters::VirtualAnalogDiodeLadderFilter VADiodeLadderFilter;
//%include "Analog/VirtualAnalogDiodeLadderFilter.hpp"
//%include "Analog/VirtualAnalogStateVariableFilter.hpp"
//%include "Analog/VoltageControlledFilter.hpp"
//%include "Analog/VoltageControlledOscillator.hpp"

//%include "Analog/VAWDFCompressor.hpp"
//%include "Analog/VAWDFDiodeClipper.hpp"
//%include "Analog/VAWDFPassiveLPF.hpp"
//%include "Analog/VAWDFSallenKey.hpp"
