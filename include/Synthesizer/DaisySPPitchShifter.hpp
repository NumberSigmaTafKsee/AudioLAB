#pragma once

#include "DaisySP.hpp"

namespace DaisySP::FX
{
    struct PitchShifter : public MonoFXProcessorPlugin<daisysp::PitchShifter>
    {
        PitchShifter(DspFloatType sampleRate=44100.0) : MonoFXProcessorPlugin<daisysp::PitchShifter>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_TRANSPOSE,
            PORT_DELSIZE,
            PORT_FUN,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRANSPOSE: this->SetTransposition(v); break;
                case PORT_DELSIZE: this->SetDelSize(v); break;
                case PORT_FUN: this->SetFun(v); break;                
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            float x = I;
            return A*this->Process(x);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            for(size_t i = 0; i < n; i++) _out[i] = Tick(_in[i]);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}