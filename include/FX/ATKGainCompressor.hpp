#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    struct GainCompressor : public ATKFilter
    {
        ATK::GainFilter<ATK::GainCompressorFilter<DspFloatType>> * filter;

        GainCompressor()
        {
            filter = new ATK::GainFilter<ATK::GainCompressorFilter<DspFloatType>>(2);
            this->setFilter(filter);
        }
        ~GainCompressor() {
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