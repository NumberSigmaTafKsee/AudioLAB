#pragma once

#include "DaisySP.hpp"

namespace DaisySP::PhysicalModels
{
    struct KarplusString : public GeneratorProcessorPlugin<daisysp::String>
    {
        KarplusString(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::String>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_NONLINEAR,
            PORT_BRIGHT,
            PORT_DAMP
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_NONLINEAR: this->SetNonLinearity(v); break;
                case PORT_BRIGHT: this->SetBrightness(v); break;
                case PORT_DAMP: this->SetDamping(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(I);
        }
    };
}