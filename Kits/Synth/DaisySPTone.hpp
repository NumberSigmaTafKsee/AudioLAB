#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    struct Tone : public FilterProcessorPlugin<daisysp::Tone>
    {
        Tone(DspFloatType sampleRate=44100.0) : FilterProcessorPlugin<daisysp::Tone>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
        };
        void setPort(int port, DspFloatType v) {
            float x = v;
            switch(port) {
                case PORT_FREQ: this->SetFreq(x); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            float x = I;
            return A*this->Process(x);
        }
    };
}