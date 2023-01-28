#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Filters
{
    template<int max_size,int max_block>
    struct FIR : public FilterProcessorPlugin<daisysp::FIRFilterImplGeneric<max_size,max_block>>
    {
        FIR(const float * ir, size_t len, bool rev) : FilterProcessorPlugin<daisysp::FIRFilterImplGeneric<max_size,max_block>>()
        {
            this->SetIR(ir,len,rev);
        }
        enum {
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_RESET: this->Reset(); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(I);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            this->ProcessBlock(_in.data(),_out.data(),n);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}