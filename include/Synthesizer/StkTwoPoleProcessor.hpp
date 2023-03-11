#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct TwoPole : public FilterProcessorPlugin<stk::TwoPole>
    {
        TwoPole() : FilterProcessorPlugin<stk::TwoPole>()
        {

        }
        enum {
            PORT_B0,            
            PORT_A1,            
            PORT_A2,
            PORT_RESONANCE,            
        };        
        void setPort(int port, DspFloatType v) {
            switch(port) {
            case PORT_B0: this->setB0(v); break;            
            case PORT_A1: this->setA1(v); break;            
            case PORT_A2: this->setA2(v); break;                        
            default: printf("No port %d\n",port);
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b) {
            if(port == PORT_RESONANCE) this->setResonance(a,b,true);
            else printf("No port %d\n",port);
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
    };
}