#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Envelopes
{
    struct ADSR : public GeneratorProcessor, public gam::ADSR
    {
        ADSR(DspFloatType att, DspFloatType dec, DspFloatType sus, DspFloatType rel, DspFloatType amp=1, DspFloatType crv=-4) :
            GeneratorProcessor(),
            gam::ADSR(att,dec,rel,amp,crv)
        {
            
        }
        enum {
            PORT_ATTACK,
            PORT_DECAY,
            PORT_SUSTAIN,
            PORT_RELEASE,
            PORT_AMP,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_ATTACK: this->attack(v); break;
                case PORT_DECAY: this->decay(v); break;
                case PORT_SUSTAIN: this->sustain(v); break;
                case PORT_RELEASE: this->release(v); break;
                case PORT_AMP: this->amp(v);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)();
        }
    };
}