#pragma once
#include "StkHeaders.hpp"

namespace Stk::FX
{
    struct Echo : public MonoFXProcessorPlugin<stk::Echo>
    {
        Echo(unsigned long max = (unsigned long)sampleRate) : 
            MonoFXProcessorPlugin<stk::Echo>()         
        {
            this->setMaximumDelay(max);
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
            return this->tick(I);
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output)
        {
            for(size_t i = 0; i < n; i++) output[i] = Tick(input[i]);
        }
    };
}