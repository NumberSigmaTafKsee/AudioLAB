#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct ZOscillator : public GeneratorProcessorPlugin<daisysp::ZOscillator>
    {
        ZOscillator(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::ZOscillator>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_FORMANTFREQ,
            PORT_SHAPE,
            PORT_MODE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_FORMANTFREQ: this->SetFormantFreq(v); break;
                case PORT_SHAPE: this->SetShape(v); break;
                case PORT_MODE: this->SetMode(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}    
