#pragma once
#include "ClippingFunctions.hpp"

namespace FX::WaveShaping
{
    struct ATanSoftClipper : public AmplifierProcessor
    {
        DspFloatType gain = 1.0;
        ATanSoftClipper(DspFloatType g = 1.0) : AmplifierProcessor()
        {
            gain = pow(10.0,g/20.0);
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*arctanSoftClipping(gain*I);
        }
        /*
        Eigen::VectorXf Vectorize(const Eigen::VectorXf& I) {
            Eigen::VectorXf r = atanSoftClipping(gain*I);
        }
        */
    };
}