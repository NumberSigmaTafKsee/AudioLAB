#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct PoleZero : public FilterProcessorPlugin<stk::PoleZero>
    {
        PoleZero() : FilterProcessorPlugin<stk::PoleZero>()
        {

        }
        enum {
            PORT_B0,
            PORT_B1,
            PORT_A1,            
            PORT_ALLPASS,
            PORT_BLOCKZERO,
        };        
        void setPort(int port, DspFloatType v) {
            switch(port) {
            case PORT_B0: this->setB0(v); break;
            case PORT_B1: this->setB1(v); break;
            case PORT_A1: this->setA1(v); break;            
            case PORT_ALLPASS: this->setAllpass(v); break;
            case PORT_BLOCKZERO: this->setBlockZero(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
    };
}