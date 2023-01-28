#pragma once

#include <cmath>
#include <cstdlib>
#include "fxdsp.hpp"

namespace FXDSP
{

    /*******************************************************************************
    RMSEstimator */
    template<typename T>
    struct RMSEstimator
    {
        T   avgTime;
        T   sampleRate;
        T   avgCoeff;
        T   RMS;

        RMEstimator(T avgTime, T sampleRate)
        {
            this->avgTime = avgTime;
            this->sampleRate = sampleRate;
            RMS = 1;
            avgCoeff = 0.5 * (1.0 - std::exp( -1.0 / (sampleRate * avgTime)));
        }
        void flush() {
            RMS = 1.0;
        }
        void setAvgTime(T avgTime)
        {
            this->avgTime = avgTime;
            avgCoeff = 0.5 * (1.0 - std::exp( -1.0 / (sampleRate * avgTime)));
        }
        void ProcessBlock(size_t n_samples, T * inBuffer, T * outBuffer)
        {
            #pragma omp simd
            for (unsigned i = 0; i < n_samples; ++i)
            {
                RMS += avgCoeff * ((f_abs(inBuffer[i])/RMS) - RMS);
                outBuffer[i] = RMS;
            }
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            RMS += avgCoeff * ((f_abs(inSample/RMS)) - RMS);
            return RMS;
        }
    };
}