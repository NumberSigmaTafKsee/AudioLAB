%module MoogFilters
%{
#include "SoundObject.hpp"

#include <cassert>
#include <random>
#include <chrono>
#include <complex>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "Analog/VAImprovedMoogFilter.hpp"
#include "Analog/VALadderFilter.hpp"
#include "Analog/VALadderFilter2.hpp"
#include "Analog/VAMicroTrackerMoogFilter.hpp"
#include "Analog/VAMoogCatFilter.hpp"
#include "Analog/VAMoogFilter.hpp"
#include "Analog/VAMoogFilter1.hpp"
#include "Analog/VAMoogFilter2.hpp"
#include "Analog/VAMoogFilter3.hpp"
#include "Analog/VAMoogFilter4.hpp"
#include "Analog/VAMoogFilterI.hpp"
#include "Analog/VAMoogFilterII.hpp"
#include "Analog/VAMoogFilters.hpp"
#include "Analog/VAMoogLadderFilters.hpp"
#include "Analog/VAMoogLadders.hpp"
#include "Analog/VAMoogLikeFilter.hpp"
#include "Analog/VAMoogNonLinearFilter.hpp"
#include "Analog/VAMoogNonLinearFilter2.hpp"
#include "Analog/VAMoogRKLadderFilter.hpp"
#include "Analog/VAMoogVCFFilter.hpp"
#include "Analog/VARKLadderFilter.hpp"
#include "Analog/VAKrajeskiMoogFilter.hpp"
#include "Analog/VAStilsonMoogFilter.hpp"
#include "Analog/VAStilsonMoogFilter2.hpp"

/*
#include "Analog/VAMoogHalfLadderFilter.cpp
#include "Analog/VAMoogHalfLadderFilter.hpp"
#include "Analog/VAMoogLadderFilter.cpp
#include "Analog/VAMoogLadderFilter.hpp"
*/

%}

%include "stdint.i"
%include "std_math.i"
%include "std_vector.i"
%include "std_list.i"
%include "std_map.i"
%include "lua_fnptr.i"

%include "SoundObject.hpp"


%include "Analog/VAImprovedMoogFilter.hpp"
%include "Analog/VAKrajeskiMoogFilter.hpp"
%include "Analog/VAMicroTrackerMoogFilter.hpp"

%ignore Analog::Filters::RKLadderFilter::clip;
%ignore Analog::Filters::RKLadderFilter::crossfade;
%ignore Analog::Filters::RKLadderFilter::stepRK4;
%rename Analog::Filters::RKLadderFilter::LadderFilter RKLadderFilter;
%include "Analog/VARKLadderFilter.hpp"

%ignore Analog::Filters::Moog::StilsonMoog::gaintable;
%rename Analog::Filters::Moog::StilsonMoog::StilsonMoog VAStilsonMoogFilter;
%include "Analog/VAStilsonMoogFilter.hpp"

%rename Analog::Filters::Moog::StilsonMoogFilter2::StilsonMoog VAStilsonMoogFilter2;
%include "Analog/VAStilsonMoogFilter2.hpp"

%ignore Analog::Filters::LadderFilter2::TEMP;
%ignore Analog::Filters::LadderFilter2::THERMAL_VOLT;
%ignore Analog::Filters::LadderFilter2::OVER_TWO_THERMAL_VOLT;
%ignore Analog::Filters::LadderFilter2::NUMBER_OF_FILTERS;
%include "Analog/VALadderFilter.hpp"
%include "Analog/VALadderFilter2.hpp"
%include "Analog/VAMoogCatFilter.hpp"

%include "Analog/VAMoogFilter.hpp"
%include "Analog/VAMoogFilter1.hpp"
%include "Analog/VAMoogFilter2.hpp"
%include "Analog/VAMoogFilter3.hpp"

%rename Analog::MoogFilters::MoogFilter4::MoogFilter VAMoogFilter4;
%include "Analog/VAMoogFilter4.hpp"
%include "Analog/VAMoogFilterI.hpp"
%include "Analog/VAMoogFilterII.hpp"
%include "Analog/VAMoogFilters.hpp"


//%include "Analog/VAMoogLikeFilter.hpp"
%include "Analog/VAMoogNonLinearFilter.hpp"

%rename Analog::Filters::Moog::NonLinear2::MoogFilter VANonLinearMoogFilter;
%include "Analog/VAMoogNonLinearFilter2.hpp"


%rename Analog::Filters::Moog::MoogVCF::MoogVCF VAMoogVoltageControlledFilter;
%include "Analog/VAMoogVCFFilter.hpp"

//%include "Analog/VAMoogHalfLadderFilter.hpp"
//%include "Analog/VAMoogLadder.hpp"

//%rename Analog::Filters::MoogLadder::MoogLadder MoogLadder1;
//%include "Analog/VAMoogLadderFilter.hpp"


