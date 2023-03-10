#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Analysis
{
    struct MaxAbs : public FunctionProcessor,public gam::MaxAbs<double>
    {    
        MaxAbs(int winsize=256)
        : FunctionProcessor(),
            gam::MaxAbs<double>(winsize)
        {

        }

        double Tick(double I, double A=1, double X=0, double Y=0)
        {                        
            return A*(*this)(I);
        }
    };
}