#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Noise
{
    struct ClockedNoise : public GeneratorProcessorPlugin<daisysp::ClockedNoise>
    {
        ClockedNoise(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::ClockedNoise>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_SYNC,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_SYNC: this->Sync(); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {            
            return A*this->Process();
        }        
    };
}