#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    struct GainLimiter : public ATKFilter
    {
        ATK::GainFilter<ATK::GainLimiterFilter<DspFloatType>> * filter;

        GainLimiter()
        {
            filter = new ATK::GainFilter<ATK::GainLimiterFilter<DspFloatType>>(2);
            this->setFilter(filter);
        }
        ~GainLimiter() {
            if(filter) delete filter;
        }
        enum {
            PORT_SOFTNESS,         
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_SOFTNESS: filter->set_softness(value); break;
            }
        }
    };
}