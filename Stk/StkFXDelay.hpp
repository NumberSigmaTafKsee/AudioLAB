#pragma once
#include "StkHeaders.hpp"

namespace Stk::FX
{
    // todo: ports
    struct Delay : public FunctionProcessorPlugin<stk::Delay>
    {
        Delay(unsigned long delay=0, unsigned long maxDelay=4096) : 
            FunctionProcessorPlugin<stk::Delay>()
        {
            this->setDelay(delay);
            this->setMaximumDelay(maxDelay);
        }
        enum {
            PORT_MAXDELAY,
            PORT_DELAY,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_MAXDELAY: this->setMaximumDelay(v); break;
                case PORT_DELAY: this->setDelay(v); break;
            }
        }
        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->tick(I);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
    };
}