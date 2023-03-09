#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct Cubic : public MonoFXProcessorPlugin<stk::Cubic>
    {
        Cubic() : MonoFXProcessorPlugin<stk::Cubic>()
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->tick(I);
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
        }
    };
}