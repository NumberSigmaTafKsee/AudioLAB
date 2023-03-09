#pragma once

#include "DaisySP.hpp"
#include "Utility/samplehold.h"

namespace DaisySP::Util
{
    struct SampleHold : public FunctionProcessorPlugin<daisysp::SampleHold>
    {
        SampleHold() : FunctionProcessorPlugin<daisysp::SampleHold>()
        {
            
        }
        bool trigger = false;
        SampleHold::Mode mode    = MODE_SAMPLE_HOLD;
        enum {
            PORT_TRIGGER,
            PORT_MODE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRIGGER: trigger = (bool)v; break;
                case PORT_MODE: mode = (Mode)v; break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->Process(trigger,I,mode);
        }
    };
}