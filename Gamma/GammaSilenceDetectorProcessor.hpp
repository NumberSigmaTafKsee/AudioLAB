#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Analysis
{
    // ports
    // count
    // reset
    // done

    // X != 0       count(X)
    // Y >= 1.0     reset()

    struct SilenceDetect : public FunctionProcessorPlugin<gam::SilenceDetect>
    {
        std::function<void (SilenceDetect&)> callback = [](SilenceDetect & s) {}

        SilenceDetect(unsigned count = 1000) 
        : FunctionProcessorPlugin<gam::SilenceDetect>,
          gam::SilenceDetect(count)
        {

        }

        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {            
            if(fabs(X) > 0) this->count(X);
            if(fabs(Y) >= 1.0 || this->done()) 
            {                   
                this->reset();
            }
            if((*this)(I,A)) callback(*this);
            return A*I;
        }    
    };
}