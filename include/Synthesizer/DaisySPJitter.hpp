#pragma once

#include "DaisySP.hpp"
#include "Utility/jitter.h"

namespace DaisySP::Util
{
    struct Jitter : public GeneratorProcessorPlugin<daisysp::Jitter>
    {
        Jitter(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::Jitter>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_CPSMIN,
            PORT_CPSMAX,
            PORT_AMP,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CPSMIN: this->SetCpsMin(v); break;
                case PORT_CPSMAX: this->SetCpsMax(v); break;
                case PORT_AMP: this->SetAmp(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}