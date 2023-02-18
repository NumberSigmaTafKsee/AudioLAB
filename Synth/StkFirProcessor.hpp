#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct Fir : public FilterProcessorPlugin<stk::Fir>
    {
        Fir(std::vector<stk::StkFloat> &coeffs) : FilterProcessorPlugin<stk::Fir>()
        {
            this->setCoefficients(coeffs,true);
        }

        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->tick(I);
        }
    };
}