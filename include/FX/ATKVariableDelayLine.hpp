#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Delays
{
    // Variable Delay Line = this is fixed..
    struct VDLDelay : public ATKFilter
    {
        ATK::VariableDelayLineFilter<DspFloatType> *filter;
        VDLDelay(size_t max_delay) :         
        ATKFilter()
        {
            filter = new ATK::VariableDelayLineFilter<DspFloatType>(max_delay);
            this->setFilter(filter);
        }
        ~VDLDelay() {
            if(filter) delete filter;
        }
        
        void setPort(size_t port, DspFloatType value)  {
        
        }
    };
    
}