#pragma once

#include "ATK.hpp"

namespace Filters::AudioTK::Dynamics
{
    struct GainSwell : public ATKFilter
    {
        ATK::GainFilter<ATK::GainSwellFilter<DspFloatType>> * filter;

        GainSwell()
        {
            filter = new ATK::GainFilter<ATK::GainSwellFilter<DspFloatType>>(2);
            this->setFilter(filter);
        }
        ~GainSwell() {
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