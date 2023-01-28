#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Dynamics
{
    struct GainColoredExpander : public ATKFilter
    {
        ATK::GainFilter<ATK::GainColoredExpanderFilter<DspFloatType>> * filter;
        
        GainColoredExpander()
        {
            filter = new ATK::GainFilter<ATK::GainColoredExpanderFilter<DspFloatType>>(2);
            this->setFilter(filter);
        }
        ~GainColoredExpander() {
            if(filter) delete filter;
        }
        enum {
            PORT_SOFTNESS,
            PORT_COLOR,
            PORT_QUALITY
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_SOFTNESS: filter->set_softness(value); break;
                case PORT_COLOR: filter->set_color(value); break;
                case PORT_QUALITY: filter->set_quality(value); break;
            }
        }
    };
}