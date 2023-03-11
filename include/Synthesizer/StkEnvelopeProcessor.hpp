#pragma once

#include "StkHeaders.hpp"

namespace Stk::Envelopes
{
    struct Envelope : public GeneratorProcessorPlugin<stk::Envelope>
    {
        Envelope() : GeneratorProcessorPlugin<stk::Envelope>()
        {
            
        }
        enum {
            PORT_RATE,
            PORT_TIME,
            PORT_TARGET,
            PORT_VALUE,
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_RATE: this->setRate(value); break;
                case PORT_TIME: this->setTime(value); break;
                case PORT_TARGET: this->setTarget(value); break;
                case PORT_VALUE: this->setValue(value); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return this->tick();
        }
    };
}