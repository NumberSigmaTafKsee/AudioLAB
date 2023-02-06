#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Envelopes
{
    struct Gate : public GeneratorProcessor, public  gam::Gate
    {
        Gate(DspFloatType delay=0, DspFloatType thresh=0.801) :
            GeneratorProcessor(),
            gam::Gate(delay,thresh)
        {
            
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return (*this)();
        }
        enum {
            PORT_DELAY,
            PORT_DONE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DELAY: this->delay(v); break;
            }
        }
        DspFloatType getPort(int port) {
            if(port == PORT_DONE) return this->done(); 
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return (*this)(I);
        }
    };
}