#pragma once
#include "Amplifiers.hpp"

namespace FX::Distortion
{
    struct FoldAmp : public AmplifierProcessor
    {
        FoldAmp() : AmplifierProcessor()
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*Fold(I);
        }
    };

       
}