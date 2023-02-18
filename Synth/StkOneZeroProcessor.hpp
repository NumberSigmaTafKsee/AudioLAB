#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct OneZero : public FilterProcessorPlugin<stk::OneZero>
    {
        DspFloatType sampleRate;
        OneZero(DspFloatType sr=44100.0) : FilterProcessorPlugin<stk::OneZero>()
        {
            sampleRate = sr;
        }
        enum {
            PORT_B0,
            PORT_B1,
            PORT_ZERO,
            PORT_CUTOFF,
            PORT_HIGHPASS,
        };
        void setCutoff(DspFloatType f) {
            DspFloatType p = exp(-2*M_PI*f/sampleRate);            
            this->setZero(p);
        }
        void setHighpassCutoff(DspFloatType f) {
            DspFloatType p = -exp(-2*M_PI*f/sampleRate);
            this->setZero(p);
        }
        void setPort(int port, DspFloatType v) {
            switch(port) {
            case PORT_B0: this->setB0(v); break;
            case PORT_B1: this->setB1(v); break;
            case PORT_ZERO: this->setZero(v); break;
            case PORT_CUTOFF: setCutoff(v); break;
            case PORT_HIGHPASS: setHighpassCutoff(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
    };
}