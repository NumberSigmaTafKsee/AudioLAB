#pragma once
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include "fxdsp.hpp"
#include "fxdsp_rbj.hpp"

namespace FXDSP
{
    /* OnePoleFilter ********************************************************/
    template<typename T>
    struct OnePole
    {
        T a0;
        T b1;
        T y1;
        T cutoff;
        T sampleRate;
        Filter_t type;

        OnePole(T cutoff, T sampleRate, Filter_t type)
        {
            a0 = 1;
            b1 = 0;
            y1 = 0;
            type = type;
            sampleRate = sampleRate;
            setCutoff(cutoff);    
        }
        void setType(Filter_t type)
        {
            if (type == LOWPASS || type == HIGHPASS)
            {
                type = type;
                setCutoff(cutoff);                
            }
            else throw std::runtime_error("Onepole can only be low or highpass");            
        }
        void setCutoff(T cutoff)
        {
            this->cutoff = cutoff;
            if (type == LOWPASS)
            {
                b1 = expf(-2.0 * M_PI * (cutoff / sampleRate));
                a0 = 1.0 - b1;
            }
            else
            {
                b1 = -expf(-2.0 * M_PI * (0.5 - (cutoff / sampleRate)));
                a0 = 1.0 + b1;
            }
        }
        void setSampleRate(T sampleRate)
        {
            this->sampleRate = sampleRate;
            setCutoff(cutoff);
        }
        void setCoefficients(T * beta, T * alpha)
        {
            b1 = *beta;
            a0 = *alpha;
        }
        void flush() {
            y1 = (T)0.0;
        }
        void ProcessBlock(size_t n, T * inBuffer, T * outBuffer)
        {
            #pragma omp simd
            for (unsigned i = 0; i < n; ++i)
            {
                outBuffer[i] = y1 = inBuffer[i] * a0 + y1 * b1;
            }
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            return  y1 = inSample * a0 + y1 * b1;
        }
        T alpha() { return a0; }
        T beta()  { return b0; }

    };
}