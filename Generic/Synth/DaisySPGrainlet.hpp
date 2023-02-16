#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Noise
{
    struct Grainlet : public GeneratorProcessorPlugin<daisysp::GrainletOscillator>
    {
        Grainlet(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::GrainletOscillator>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_FORMANTFREQ,
            PORT_SHAPE,
            PORT_BLEED,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_FORMANTFREQ: this->SetFormantFreq(v); break;
                case PORT_SHAPE: this->SetShape(v); break;
                case PORT_BLEED: this->SetBleed(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}