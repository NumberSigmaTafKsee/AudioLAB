#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
    // todo: ports
    struct Iir : public FilterProcessorPlugin<stk::Iir>
    {
        Iir(std::vector<stk::StkFloat> & a, std::vector<stk::StkFloat> & b) : FilterProcessorPlugin<stk::Iir>()
        {
            this->setCoefficients(a,b,true);
        }        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->tick(I);
        }
    };
}