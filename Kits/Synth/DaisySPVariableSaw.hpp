#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct VariableSaw : public GeneratorProcessorPlugin<daisysp::VariableSawOscillator>
    {
        VariableSaw(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::VariableSawOscillator>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_PW,
            PORT_WAVESHAPE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_PW: this->SetPW(v); break;
                case PORT_WAVESHAPE: this->SetWaveshape(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}