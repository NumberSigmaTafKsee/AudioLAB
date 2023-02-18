#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Analysis
{
    struct ZeroCrossRate : public FunctionProcessor, public gam::ZeroCrossRate<DspFloatType>
    {
        ZeroCrossRate(int winsize=256)
        : FunctionProcessor(),
            gam::ZeroCrossRate<DspFloatType>(winsize)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {                        
            DspFloatType r = (*this)(I);
            return r;
        }
    };
}