#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"
#include <functional>

namespace Gamma::Analysis
{
    
    struct Inspector : public FunctionProcessorPlugin<gam::Inspector>
    {                
        std::function<void (Inspector &, unsigned)> callback = [](Inspector&,unsigned){};

        Inspector(int winsize=256)
        : FunctionProcessorPlugin<gam::Inspector>()
        {
            
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {                        
            unsigned r = (*this)(I);
            if(r != 0) callback(*this,r);
            return I;
        }
    };
};