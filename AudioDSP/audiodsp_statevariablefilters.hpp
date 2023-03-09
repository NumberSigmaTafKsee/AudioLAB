

//#include "Analog/VAStateVariableCombFilter.hpp"

#include "Analog/VAStateVariableFilter1.hpp"
#include "Analog/VAStateVariableFilter2.hpp"
#include "Analog/VAStateVariableFilters.hpp"
#include "Analog/VASVF.hpp"
#include "Analog/VASVFChamberlinFilter.hpp"
#include "Analog/VASVFFilter.hpp"
#include "Analog/VASVFSmoother.hpp"
#include "Analog/VASVSmoothFilter.hpp"
#include "Analog/VASVStateVariableFilter.hpp"


/*
%include "Analog/VAGenSVF.hpp"
%include "Analog/VADiodeLadderFilter2.hpp"

///%include "Analog/VAStateVariableFilter1.hpp"
///%rename Analog::Filters::StateVariableFilter2::StateVariableFilter VAStateVariableFilter2;
///%include "Analog/VAStateVariableFilter2.hpp"

%rename Analog::Filters::SVF::AnalogSVF VASVF;
%include "Analog/VASVF.hpp"
%include "Analog/VASVFChamberlinFilter.hpp"
%rename Analog::Filters::SVF::StateVariableFilter VASVFFilter;
%include "Analog/VASVFFilter.hpp"
//%include "Analog/VASVStateVariableFilter.hpp"
//%include "Analog/VAStateVariableFilters.hpp"
//%include "Analog/VASVStateVariableFilter.hpp"
