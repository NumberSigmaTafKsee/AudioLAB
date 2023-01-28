#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Dynamics
{
    struct GainExpander : public ATKFilter
    {
        ATK::GainFilter<ATK::GainExpanderFilter<DspFloatType>> * filter;

        GainExpander()
        {
            filter = new ATK::GainFilter<ATK::GainExpanderFilter<DspFloatType>>(2);
            this->setFilter(filter);
        }
        ~GainExpander() {
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