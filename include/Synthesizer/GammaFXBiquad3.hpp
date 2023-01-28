#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::FX
{
    struct Biquad3 : public FunctionProcessorPlugin<gam::Biquad3<DspFloatType>>
    {
        Biquad3() : FunctionProcessorPlugin<gam::Biquad3<DspFloatType>>()
        {

        }
        enum {
            PORT_FREQ,
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_FREQ: this->freq(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)(I);
        }
    };
}