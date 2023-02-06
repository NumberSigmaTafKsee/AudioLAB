#pragma once

#include "ATK.hpp"

namespace AudioTK
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

    struct MonoGainExpander : public ATKMonoFilter
    {
        ATK::GainFilter<ATK::GainExpanderFilter<DspFloatType>> * filter;

        MonoGainExpander() : ATKMonoFilter()
        {
            filter = new ATK::GainFilter<ATK::GainExpanderFilter<DspFloatType>>(1);
            this->setFilter(filter);
        }
        ~MonoGainExpander() {
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