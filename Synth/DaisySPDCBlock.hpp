#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Util
{
    struct DCBlock : public FilterProcessorPlugin<daisysp::DcBlock>
    {
        DCBlock(DspFloatType sampleRate=44100.0) : FilterProcessorPlugin<daisysp::DcBlock>()
        {
            this->Init(sampleRate);
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