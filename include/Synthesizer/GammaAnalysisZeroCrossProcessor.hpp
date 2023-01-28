#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Analysis
{
    // ports
    // none

    // A = amp
    // X = none
    // Y = none
    struct ZeroCross : public FunctionProcessor, public gam::ZeroCross<DspFloatType>
    {
        ZeroCross(DspFloatType prev = 0.0)
        : FunctionProcessor(),
            gam::ZeroCross<DspFloatType>(prev)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {                    
            return (*this)(I);            
        }
    };
}