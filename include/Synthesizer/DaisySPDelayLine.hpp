#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Util
{
    template<int size=1024>
    struct DelayLine : public MonoFXProcessorPlugin<daisysp::DelayLine<DspFloatType,size>>
    {
        DelayLine(size_t delay) : MonoFXProcessorPlugin<daisysp::DelayLine<DspFloatType,size>>()
        {
            this->Init();
            this->SetDelay(delay);
        }
        enum {
            PORT_DELAY,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DELAY: this->SetDelay(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=0) {
            DspFloatType x = this->Read();
            this->Write( X*I + Y*x);
            return A*x;
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            for(size_t i = 0; i < n; i++) _out[i] = Tick(_in[i]);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}