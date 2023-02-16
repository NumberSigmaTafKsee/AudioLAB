#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::FX
{
    struct Burst : public FunctionProcessorPlugin<gam::Burst>
    {
        Burst() : FunctionProcessorPlugin<gam::Burst>()
        {

        }
        enum {
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_RESET) this->reset();
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)();
        }
    };
}