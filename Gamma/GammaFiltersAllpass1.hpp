#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Filters
{
    struct AllPass1 : public FilterProcessorPlugin<gam::AllPass1<DspFloatType>>
    {
        AllPass1() : FilterProcessorPlugin<gam::AllPass1<DspFloatType>>()
        {

        }
        
        enum {
            PORT_CUTOFF,
            PORT_CUTOFFF,
            PORT_RESET,
            PORT_HIGH,
            PORT_LOW,
            PORT_FREQ
        };
        int type = PORT_LOW;
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: this->freq(v); break;
                case PORT_CUTOFFF: this->freqF(v); break;
                //case PORT_RESET: this->reset(); break;
                case PORT_HIGH: type = PORT_HIGH; break;
                case PORT_LOW: type = PORT_LOW; break;
            }
        }
        DspFloatType getPort(int port) {
            if(port == PORT_FREQ) return this->freq();
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            if(type == PORT_LOW) return A*this->low(I);
            return A*this->high(I);
        }
    };
}