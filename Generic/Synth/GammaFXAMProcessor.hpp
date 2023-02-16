#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::FX
{
    struct AM : public FunctionProcessorPlugin<gam::AM<DspFloatType>>
    {
        AM() : FunctionProcessorPlugin<gam::AM<DspFloatType>>()
        {

        }
        enum {
            PORT_DEPTH
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DEPTH: this->depth(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)(I);
        }
    };
}