#pragma once
#include "ClipFunctions.hpp"

namespace FX::Distortion::Clip
{
    struct SerpentCurve : public AmplifierProcessor
    {
        DspFloatType gain = 1;
        SerpentCurve(DspFloatType g = 1.0) : AmplifierProcessor()
        {
            gain = pow(10.0,g/20.0);
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*serpent_curve(I,gain);
        }
    };
}