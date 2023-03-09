#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct OscillatorBank : public GeneratorProcessorPlugin<daisysp::OscillatorBank>
    {
        OscillatorBank(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::OscillatorBank>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_AMP,
            PORT_GAIN,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_GAIN: this->SetGain(v); break;
            }
        }
        void setPort2(int port, DspFloatType x, DspFloatType y) {
            switch(port) {
                case PORT_AMP: this->SetSingleAmp(x,y);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}