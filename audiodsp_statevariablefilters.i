%module StateVariable
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

//#include "Analog/VAStateVariableCombFilter.hpp"
#include "Analog/VirtualAnalogStateVariableFilter.hpp"
#include "Analog/VAGenSVF.hpp"
#include "Analog/VAStateVariableFilter1.hpp"
#include "Analog/VAStateVariableFilter2.hpp"
#include "Analog/VAStateVariableFilters.hpp"
#include "Analog/VASVF.hpp"
#include "Analog/VASVFChamberlinFilter.hpp"
#include "Analog/VASVFFilter.hpp"
#include "Analog/VASVFSmoother.hpp"
#include "Analog/VASVSmoothFilter.hpp"
#include "Analog/VASVStateVariableFilter.hpp"

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"



%include "Analog/VirtualAnalogStateVariableFilter.hpp"

%include "Analog/VAGenSVF.hpp"
%rename Analog::Filters::StateVariableFilter2::StateVariableFilter VAStateVariableFilter2;
%include "Analog/VAStateVariableFilter2.hpp"

%rename Analog::Filters::SVF::AnalogSVF VASVF;
%include "Analog/VASVF.hpp"
%include "Analog/VASVFChamberlinFilter.hpp"
%rename Analog::Filters::SVF::StateVariableFilter VASVFFilter;
%include "Analog/VASVFFilter.hpp"
//%include "Analog/VASVStateVariableFilter.hpp"
//%include "Analog/VAStateVariableFilters.hpp"

