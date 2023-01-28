#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Delays
{
    // Multiple Universal Fixed Delay Line
    struct MUFDLDelay : public ATKFilter
    {
        ATK::MultipleUniversalFixedDelayLineFilter<DspFloatType,2> *filter;
        MUFDLDelay(size_t max_delay) :         
        ATKFilter()
        {
            filter = new ATK::MultipleUniversalFixedDelayLineFilter<DspFloatType,2>(max_delay);
            this->setFilter(filter);
        }
        ~MUFDLDelay() {
            if(filter) delete filter;
        }
        enum {
            PORT_DELAY1,
            PORT_DELAY2,
            PORT_BLEND1,
            PORT_BLEND2,
            PORT_FEEDBACK,            
            PORT_FEEDFORWARD,            
        };
        void setPort(size_t port, DspFloatType value)  {
            switch(port)
            {
                case PORT_DELAY1: filter->set_delay(0,value); break;
                case PORT_DELAY2: filter->set_delay(1,value); break;
                case PORT_BLEND1: filter->set_blend(0,value); break;
                case PORT_BLEND2: filter->set_blend(1,value); break;
                case PORT_FEEDBACK: filter->set_feedback(0,1,value); break;                
                case PORT_FEEDFORWARD: filter->set_feedforward(0,1,value); break;                
            }
        }
    };
}