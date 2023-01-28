#pragma once

#include "DaisySP.hpp"
#include "DaisySP/Effects/wavefolder.h"

namespace DaisySP::FX
{
    struct WaveFolder : public FunctionProcessorPlugin<daisysp::Wavefolder>
    {
        WaveFolder() : FunctionProcessorPlugin<daisysp::Wavefolder>()
        {
            this->Init();
        }
        enum {
            PORT_GAIN,
            PORT_OFFSET,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_GAIN: this->SetGain(v); break;
                case PORT_OFFSET: this->SetOffset(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(I);
        }
    };
}