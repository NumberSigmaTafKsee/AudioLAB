#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    struct Mode : public FilterProcessorPlugin<daisysp::Mode>
    {
        Mode(DspFloatType sampleRate=44100.0) : FilterProcessorPlugin<daisysp::Mode>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_CLEAR,
            PORT_FREQ,
            PORT_Q,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CLEAR: this->Clear(); break;
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_Q: this->SetQ(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A* this->Process(I);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            for(size_t i = 0; i < n; i++) _out[i] = Tick(_in[i]);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}