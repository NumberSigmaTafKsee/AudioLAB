#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    struct Comb : public FilterProcessorPlugin<daisysp::Comb>
    {
        std::vector<float> buffer;
        Comb(DspFloatType sampleRate=44100.0) : FilterProcessorPlugin<daisysp::Comb>()
        {
            buffer.resize(1024);
            this->Init(sampleRate,buffer.data(),1024);
        }
        enum {
            PORT_PERIOD,
            PORT_FREQ,
            PORT_REVTIME,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_PERIOD: this->SetPeriod(v); break;
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