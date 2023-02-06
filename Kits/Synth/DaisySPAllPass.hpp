#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    struct AllPass : public FilterProcessorPlugin<daisysp::Allpass>
    {
        std::vector<float> buffer;
        AllPass(DspFloatType sampleRate=44100.0) : FilterProcessorPlugin<daisysp::Allpass>()
        {
            buffer.resize(1024);
            this->Init(sampleRate,buffer.data(),buffer.size());
        }
        enum {
            PORT_FREQ,
            PORT_REVTIME
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_REVTIME: this->SetRevTime(v); break;
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