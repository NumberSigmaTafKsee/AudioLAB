#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    struct Simple : public GeneratorProcessorPlugin<stk::Simple>
    {
        Simple() : GeneratorProcessorPlugin<stk::Simple>()
        {

        }
        enum {
            PORT_FREQ,
            PORT_KEYON,
            PORT_KEYOFF,
            PORT_NOTEON,
            PORT_NOTEOFF,
            PORT_CC,
        };
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick(I);
        }
    };
}