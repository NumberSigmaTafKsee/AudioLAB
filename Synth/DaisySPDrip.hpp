#pragma once

#include "DaisySP.hpp"

namespace DaisySP::PhysicalModels
{
    struct Drip : public GeneratorProcessorPlugin<daisysp::Drip>
    {
        Drip(DspFloatType deattack,DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::Drip>()
        {
            this->Init(sampleRate,deattack);
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=0) {
            return A*this->Process((bool)I);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            std::vector<float> _in(n),_out(n);
            for(size_t i = 0; i < n; i++) _in[i] = in[i];
            for(size_t i = 0; i < n; i++) _out[i] = Tick(_in[i]);
            for(size_t i = 0; i < n; i++) out[i] = _out[i];
        }
    };
}