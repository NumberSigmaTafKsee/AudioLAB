#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::FX
{
    struct Pluck : public GeneratorProcessorPlugin<gam::Pluck>
    {
        Pluck() : GeneratorProcessorPlugin<gam::Pluck>()
        {

        }
        enum {
            PORT_RESET,
            PORT_FREQ,
            PORT_DECAY,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_RESET: this->reset(); break;
                case PORT_FREQ: this->freq(v); break;
                case PORT_DECAY: this->env.decay(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)(I);
        }
    };
}