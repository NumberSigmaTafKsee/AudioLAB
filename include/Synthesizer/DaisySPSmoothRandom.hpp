#pragma once

#include "DaisySP.hpp"
#include "Utility/smooth_random.h"

namespace DaisySP::Util
{
    struct SmoothRandom : public GeneratorProcessorPlugin<daisysp::SmoothRandomGenerator>
    {
        SmoothRandom(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::SmoothRandomGenerator>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->Process();
        }
    };
}