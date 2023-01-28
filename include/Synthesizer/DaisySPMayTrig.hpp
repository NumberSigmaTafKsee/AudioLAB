#pragma once

#include "DaisySP.hpp"
#include "Utility/maytrig.h"

namespace DaisySP::Util
{
    struct MayTrig : public FunctionProcessorPlugin<daisysp::Maytrig>
    {
        std::function<void (MayTrig * p)> callback;
        MayTrig() : FunctionProcessorPlugin<daisysp::Maytrig>()
        {
            callback = [](MayTrig * ptr){};
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            if( (bool)this->Process(I) == true ) callback(this);
            return A*I;
        }
    };
}