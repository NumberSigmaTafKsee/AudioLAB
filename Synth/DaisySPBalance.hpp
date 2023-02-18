#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Dynamics
{
    struct Balance : public FunctionProcessorPlugin<daisysp::Balance>
    {
        Balance(DspFloatType sampleRate=44100.0) : FunctionProcessorPlugin<daisysp::Balance>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_CUTOFF,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: this->SetCutoff(v); break;
            }
        }        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(I,X);
        }
        void ProcessBlock(size_t n, float * in, float * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
    };
}