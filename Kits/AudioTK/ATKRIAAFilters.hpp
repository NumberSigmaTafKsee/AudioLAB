#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::EQ
{  
    struct RIAAFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::RIAACoefficients<DspFloatType>> * filter;

        RIAAFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::RIAACoefficients<DspFloatType>>(2);
        }
        ~RIAAFilter() {
            if(filter) delete filter;
        }
        enum {
            
        };
        void setPort(int port, DspFloatType v) {
            
        }
    };
    struct InverseRIAAFilter: public ATKFilter
    {
        ATK::IIRFilter<ATK::InverseRIAACoefficients<DspFloatType>> * filter;

        InverseRIAAFilter() : ATKFilter()
        {
            filter = new ATK::IIRFilter<ATK::InverseRIAACoefficients<DspFloatType>>(2);
        }
        ~InverseRIAAFilter() {
            if(filter) delete filter;
        }
        enum {
            
        };
        void setPort(int port, DspFloatType v) {
            
        }
    };
}    