#pragma once

#include "DaisySP.hpp"

namespace DaisySP::FX
{
    struct BitCrush : public MonoFXProcessorPlugin<daisysp::Bitcrush>
    {
        BitCrush(DspFloatType sampleRate=44100.0) : MonoFXProcessorPlugin<daisysp::Bitcrush>( )
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_DEPTH,
            PORT_RATE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_DEPTH: this->SetBitDepth(v); break;
                case PORT_RATE: this->SetCrushRate(v); break;
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