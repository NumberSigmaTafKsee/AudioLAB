
#pragma once
#include "StkHeaders.hpp"

namespace Stk::FX
{
    struct PitShift : public FunctionProcessorPlugin<stk::PitShift>
    {
        PitShift() : FunctionProcessorPlugin<stk::PitShift>()
        {

        }
        enum {
            PORT_CLEAR,
            PORT_SHIFT,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CLEAR: this->clear(); break;
                case PORT_SHIFT: this->setShift(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
    };
}