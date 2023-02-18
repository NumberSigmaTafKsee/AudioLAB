#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct FM2 : public GeneratorProcessorPlugin<daisysp::Fm2>
    {
        FM2() : GeneratorProcessorPlugin<daisysp::Fm2>()
        {

        }
        enum {
            PORT_FREQ,
            PORT_RATIO,
            PORT_INDEX,
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v)
        {
            switch(port) {
                case PORT_FREQ: this->SetFrequency(v); break;
                case PORT_RATIO: this->SetRatio(v); break;
                case PORT_INDEX: this->SetIndex(v); break;
                case PORT_RESET: this->Reset();
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {            
            return A*this->Process();
        }        
    };
}