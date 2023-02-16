#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::FX
{
    struct Quantizer : public FunctionProcessorPlugin<gam::Quantizer<DspFloatType>>
    {
        Quantizer() : FunctionProcessorPlugin<gam::Quantizer<DspFloatType>>()
        {

        }
        enum {
            PORT_FREQ,
            PORT_PERIOD,
            PORT_STEP,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->freq(v); break;
                case PORT_PERDIO: this->period(v); break;
                case PORT_STEP: this->step(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)(I);
        }
    };    
}