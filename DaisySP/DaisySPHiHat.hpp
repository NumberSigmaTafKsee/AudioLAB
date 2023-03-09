#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Drums
{
    struct HiHat : public GeneratorProcessorPlugin<daisysp::HiHat<>>
    {
        HiHat(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::HiHat<>>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_TRIG,
            PORT_SUSTAIN,
            PORT_ACCENT,
            PORT_FREQ,
            PORT_TONE,
            PORT_DECAY,
            PORT_NOISINESS,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRIG: this->Trig(); break;
                case PORT_SUSTAIN: this->SetSustain(v); break;
                case PORT_ACCENT: this->SetAccent(v); break;
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_TONE: this->SetTone(v); break;
                case PORT_DECAY: this->SetDecay(v); break;
                case PORT_NOISINESS: this->SetNoisiness(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(I);
        }
    };
}