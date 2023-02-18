#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Noise
{
    struct WhiteNoise : public GeneratorProcessorPlugin<daisysp::WhiteNoise>
    {
        WhiteNoise() :GeneratorProcessorPlugin<daisysp::WhiteNoise>()
        {
            this->Init();
        }
        enum {
            PORT_AMP,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_AMP: this->SetAmp(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}
