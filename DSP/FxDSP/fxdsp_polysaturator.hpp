#pragma once

#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "fxdsp.hpp"

namespace FXDSP
{
    template<typename T>
    struct PolySaturator
    {
        T a;
        T b;
        T n;

        PolySaturator(T n)
        {
            setN(n);
        }
        void setN(T n)
        {
            saturator->a = std::pow(1./n, 1./n);
            saturator->b = (n + 1) / n;
            saturator->n = n;
        }
        void ProcessBlock(size_t n, T * in_buffer, T * out_buffer)
        {
            T buf[n_samples];
            VectorScalarMultiply(buf, (T*)in_buffer, saturator->a, n_samples);
            VectorAbs(buf, buf, n_samples);
            VectorPower(buf, buf, saturator->n, n_samples);
            VectorScalarAdd(buf, buf, -saturator->b, n_samples);
            VectorNegate(buf, buf, n_samples);
            VectorVectorMultiply(out_buffer, in_buffer, buf, n_samples);
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            return -(std::pow(std::fabs(saturator->a * in_sample), saturator->n) - saturator->b) * in_sample;
        }
    };
}