#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Envelopes
{
    struct AD : public GeneratorProcessor, public gam::AD
    {
        AD(DspFloatType att, DspFloatType dec, DspFloatType amp=1, doubl crv=-4) 
        : GeneratorProcessor(),
            gam::AD(att,dec,amp,crv)
        {
            
        }
        enum {
            PORT_ATTACK,
            PORT_DECAY,
            PORT_AMP,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_ATTACK: this->attack(v); break;
                case PORT_DECAY: this->decay(v); break;
                case PORT_AMP: this->amp(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return (*this)();
        }
    };
}