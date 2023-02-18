#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Noise
{
    struct Particle : public GeneratorProcessorPlugin<daisysp::Particle>
    {
        Particle (DspFloatType sampleRate=44100.0): GeneratorProcessorPlugin<daisysp::Particle>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_RES,
            PORT_RANDOMFREQ,
            PORT_DENSITY,
            PORT_GAIN,
            PORT_SPREAD,
            PORT_SYNC,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_RES: this->SetResonance(v); break;
                case PORT_RANDOMFREQ: this->SetRandomFreq(v); break;
                case PORT_DENSITY: this->SetDensity(v); break;
                case PORT_GAIN: this->SetGain(v); break;
                case PORT_SPREAD: this->SetSpread(v); break;
                case PORT_SYNC: this->SetSync(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}