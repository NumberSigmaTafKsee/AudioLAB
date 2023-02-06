#pragma once

#include "DaisySP.hpp"

namespace DaisySP::FX
{
    struct Decimator : public MonoFXProcessorPlugin<daisysp::Decimator>
    {
        Decimator() :  MonoFXProcessorPlugin<daisysp::Decimator>()
        {
            this->Init();
        }
        enum {
            PORT_DOWNSAMPLE_FACTOR,
            PORT_BITCRUSH_FACTOR,
            PORT_BITS_TO_CRUSH,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DOWNSAMPLE_FACTOR: this->SetDownsampleFactor(v); break;
                case PORT_BITCRUSH_FACTOR: this->SetBitcrushFactor(v); break;
                case PORT_BITS_TO_CRUSH: this->SetBitsToCrush(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return this->Process(I);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            for(size_t i = 0; i < n; i++) _out[i] = Tick(_in[i]);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}