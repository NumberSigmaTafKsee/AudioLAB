#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct TwoZero : public FilterProcessorPlugin<stk::TwoZero>
    {
        TwoZero() : FilterProcessorPlugin<stk::TwoZero>()
        {

        }
        enum {
            PORT_B0,
            PORT_B1,            
            PORT_B2,
            PORT_NOTCH,
            PORT_LASTOUT,       
        };        
        void setPort(int port, DspFloatType v) {
            switch(port) {
            case PORT_B0: this->setB0(v); break;
            case PORT_B1: this->setB1(v); break;
            case PORT_B2: this->setB2(v); break;                        
            default: printf("No port %d\n",port);
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b) {
            if(port == PORT_NOTCH) this->setNotch(a,b);
            else printf("No port %d\n",port);
        }
        DspFloatType getPort(int port) {
            if(port == PORT_LASTOUT) return this->lastOut();
            printf("No port %d\n",port);
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
    };
}