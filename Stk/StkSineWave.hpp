#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    struct SineWave : public GeneratorProcessorPlugin<stk::SineWave>
    {
        SineWave() : GeneratorProcessorPlugin<stk::SineWave>()
        {

        }

        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick();
        }
    };
}