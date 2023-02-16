#pragma once

#include "DaisySP.hpp"

namespace DaisySP::FX
{
    struct SampleRateReducer : public FunctionProcessorPlugin<daisysp::SampleRateReducer>
    {
        SampleRateReducer() : FunctionProcessorPlugin<daisysp::SampleRateReducer>()
        {
            this->Init();
        }
        enum {
            PORT_FREQ,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->Process(I);
        }
    };
}