#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Analysis
{
    // ports
    // inv

    // A = amp
    // X = hi
    // Y = lo

    struct Threshold : public FunctionProcessor, public gam::Threshold<DspFloatType>
    {        
        Threshold(DspFloatType thresh, DspFloatType f = 10.0) 
        : FunctionProcessor(),gam::Threshold<DspFloatType>(thresh,f)
        {

        }
        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1.0, DspFloatType Y=-1.0)
        {                        
            DspFloatType x = (*this)(I);
            if(x == 0) return -A*I;
            return A*I;
        }
    };

    struct ModThreshold : public FunctionProcessor, public gam::Threshold<DspFloatType>
    {        
        ModThreshold(DspFloatType thresh, DspFloatType f = 10.0) 
        : FunctionProcessor(),
            gam::Threshold<DspFloatType>(thresh,f)
        {

        }        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1.0, DspFloatType Y=-1.0)
        {                        
            DspFloatType x = (*this)(I,X,Y);            
            return A*x;
        }
    };

    struct InvThreshold : public FunctionProcessor, public gam::Threshold<DspFloatType>
    {        
        InvThreshold(DspFloatType thresh, DspFloatType f = 10.0) 
        : FunctionProcessor(),
            gam::Threshold<DspFloatType>(thresh,f)
        {

        }
        
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1.0, DspFloatType Y=-1.0)
        {                        
            DspFloatType x = this->inv(I);            
            if(x == 1) return -A*I;
            return A*I;
        }
    };
}