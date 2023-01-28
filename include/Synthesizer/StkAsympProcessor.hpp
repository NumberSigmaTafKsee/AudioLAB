#pragma once

#include "StkHeaders.hpp"

namespace Stk::Envelopes
{
    struct ASymp : public GeneratorProcessorPlugin<stk::Asymp>
    {
        ASymp() : GeneratorProcessorPlugin<stk::Asymp>()
        {
            
        }
        enum {
            PORT_TAU,
            PORT_TIME,
            PORT_T60,
            PORT_VALUE,            
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_TAU: this->setTau(value); break;
                case PORT_TIME: this->setTime(value); break;
                case PORT_T60: this->setT60(value); break;
                case PORT_VALUE: this->setValue(value); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return this->tick();
        }
    };
}