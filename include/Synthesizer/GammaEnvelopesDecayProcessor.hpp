#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Envelopes
{
    struct Decay : public GeneratorProcessor, public gam::Decay
    {
        Decay(DspFloatType decay, DspFloatType val=1) :
            GeneratorProcessor(),
            gam::Decay(decay,val)
        {
            
        }
        enum {
            PORT_DECAY,
            PORT_RESET,
            PORT_FINISH,
            PORT_VALUE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DECAY: this->decay(v); break;
                case PORT_RESET: this->reset(v); break;
                case PORT_FINISH: this->finish(v); break;
                case PORT_VALUE: this->value(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return (*this)();
        }
    };
}