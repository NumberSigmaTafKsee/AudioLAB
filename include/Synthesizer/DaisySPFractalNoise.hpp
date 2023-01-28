#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Noise
{
    template<int order=3>
    struct FractalNoise : public GeneratorProcessorPlugin<daisysp::FractalRandomGenerator<DspFloatType,order>>
    {
        FractalNoise(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::FractalRandomGenerator<DspFloatType,order>>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_COLOR,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_COLOR: this->SetColor(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}