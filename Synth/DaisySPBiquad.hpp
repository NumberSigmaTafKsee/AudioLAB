#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    struct Biquad : public FilterProcessorPlugin<daisysp::Biquad>
    {
        Biquad(DspFloatType sampleRate=44100.0) : FilterProcessorPlugin<daisysp::Biquad>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_CUTOFF,
            PORT_RES,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: this->SetCutoff(v); break;
                case PORT_RES: this->SetRes(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {            
            return A*this->Process(I);
        }
        void ProcessBlock(size_t n, float * in, float * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
    };
}