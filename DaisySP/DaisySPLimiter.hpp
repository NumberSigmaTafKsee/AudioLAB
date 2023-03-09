#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Dynamics
{   
    struct Limiter : public MonoFXProcessorPlugin<daisysp::Limiter>
    {
        float pre_gain = 1.0;
        Limiter() : MonoFXProcessorPlugin<daisysp::Limiter>()
        {
            this->Init();
        }
        enum {
            PORT_PREGAIN,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_PREGAIN: pre_gain = v; break;
            }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> temp(n);
            for(size_t i = 0; i < n; i++) temp[i] = in[i];            
            dynamic_cast<daisysp::Limiter*>(this)->ProcessBlock(temp.data(),n,pre_gain);
            for(size_t i = 0; i < n; i++) out[i] = temp[i];
        }
    };
}