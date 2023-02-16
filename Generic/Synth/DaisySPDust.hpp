#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Noise
{
    struct Dust : public GeneratorProcessorPlugin<daisysp::Dust>
    {
        Dust() : GeneratorProcessorPlugin<daisysp::Dust>()
        {

        }
        enum {
            PORT_DENSITY,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DENSITY: this->SetDensity(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=0) {
            return A*this->Process();
        }       
    };  
}