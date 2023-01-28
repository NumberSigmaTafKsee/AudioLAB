#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Dynamics
{
    struct GainMaxExpander : public ATKFilter
    {
        using Filter = ATK::GainFilter<ATK::GainMaxExpanderFilter<DspFloatType>>;
        Filter * filter;

        GainMaxExpander()
        {
            filter = new Filter(2);
            this->setFilter(filter);
        }
        ~GainMaxExpander() {
            if(filter) delete filter;
        }
        enum {
            PORT_SOFTNESS,
            PORT_REDUCTION,
            PORT_REDUCTIONDB,            
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_SOFTNESS: filter->set_softness(value); break;
                case PORT_REDUCTION: filter->set_max_reduction(value); break;
                case PORT_REDUCTIONDB: filter->set_max_reduction_db(value); break;
            }
        }

    };
}