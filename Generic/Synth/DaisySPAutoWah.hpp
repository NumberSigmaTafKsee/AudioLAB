#pragma once

#include "DaisySP.hpp"

namespace DaisySP::FX
{
    struct AutoWah : public MonoFXProcessorPlugin<daisysp::Autowah>
    {
        AutoWah(DspFloatType sampleRate=44100.0) : MonoFXProcessorPlugin<daisysp::Autowah>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_WAH,
            PORT_DRYWET,
            PORT_LEVEL
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_WAH: this->SetWah(v); break;
                case PORT_DRYWET: this->SetDryWet(v); break;
                case PORT_LEVEL: this->SetLevel(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(I);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            for(size_t i = 0; i < n; i++) _out[i] = Tick(_in[i]);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}