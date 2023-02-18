#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct OnePole : public FilterProcessorPlugin<stk::OnePole>
    {
        DspFloatType sampleRate;
        OnePole(DspFloatType sr=44100.0) : FilterProcessorPlugin<stk::OnePole>()
        {
            sampleRate=sr;
        }
        enum {
            PORT_B0,
            PORT_A1,
            PORT_POLE,
            PORT_CUTOFF,
            PORT_HIGHPASS,
        };
        void setCutoff(DspFloatType f) {
            DspFloatType p = std::exp(-2*M_PI*f/sampleRate);
            this->setPole(p);            
        }
        void setHighpassCutoff(DspFloatType f) {
            DspFloatType p = -std::exp(-2*M_PI*f/sampleRate);
            this->setPole(p);            
        }
        void setPort(int port, DspFloatType v) {
            switch(port) {
            case PORT_B0: this->setB0(v); break;
            case PORT_A1: this->setA1(v); break;
            case PORT_POLE: this->setPole(v); break;
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