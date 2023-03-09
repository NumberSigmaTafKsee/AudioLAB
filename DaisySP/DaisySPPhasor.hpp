#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Control
{
    struct Phasor : public GeneratorProcessorPlugin<daisysp::Phasor>
    {
        Phasor(DspFloatType freq, DspFloatType phase=0.0,DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::Phasor>()
        {
            this->Init(sampleRate,freq,phase);
        }
        enum {
            PORT_FREQ,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}
