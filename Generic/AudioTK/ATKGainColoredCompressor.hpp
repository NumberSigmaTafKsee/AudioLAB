#pragma once

#include "ATK.hpp"

namespace AudioTK
{
    struct GainColoredCompressor : public ATKFilter
    {
        using Filter = ATK::GainFilter<ATK::GainColoredCompressorFilter<DspFloatType>>;
        Filter * filter;
        
        GainColoredCompressor()
        {
            filter = new Filter(2);
            this->setFilter(filter);
        }
        ~GainColoredCompressor() {
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