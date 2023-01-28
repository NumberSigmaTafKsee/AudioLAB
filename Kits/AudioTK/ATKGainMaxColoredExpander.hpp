#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Dynamics
{
    struct GainMaxColoredExpander : public ATKFilter
    {
        ATK::GainFilter<ATK::GainMaxColoredExpanderFilter<DspFloatType>> * filter;

        GainMaxColoredExpander()
        {
            filter = new ATK::GainFilter<ATK::GainMaxColoredExpanderFilter<DspFloatType>>(2);
            this->setFilter(filter);
        }
        ~GainMaxColoredExpander() {
            if(filter) delete filter;
        }
        enum {
            PORT_SOFTNESS,
            PORT_COLOR,
            PORT_QUALITY,
            PORT_REDUCTION,
            PORT_REDUCTIONDB,            
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_SOFTNESS: filter->set_softness(value); break;
                case PORT_COLOR: filter->set_color(value); break;
                case PORT_QUALITY: filter->set_quality(value); break;
                case PORT_REDUCTION: filter->set_max_reduction(value); break;
                case PORT_REDUCTIONDB: filter->set_max_reduction_db(value); break;             
            }
        }
    };
}