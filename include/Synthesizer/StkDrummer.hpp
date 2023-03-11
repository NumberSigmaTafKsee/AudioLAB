#pragma once
#include "StkHeaders.hpp"

namespace Stk::Drummer
{
     // todo: ports
    struct Drummer : public GeneratorProcessorPlugin<stk::Drummer>
    {
        Drummer() : GeneratorProcessorPlugin<stk::Drummer>()
        {

        }

        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick();
        }
    };
}