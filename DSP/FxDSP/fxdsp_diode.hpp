#pragma once

#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "fxdsp.hpp"

namespace FXDSP
{
    typedef enum _bias_t
    {
        /** Pass positive signals, clamp netagive signals to 0 */
        FORWARD_BIAS,

        /** Pass negative signals, clamp positive signals to 0 */
        REVERSE_BIAS,

        /** Full-wave rectification. */
        FULL_WAVE
    }bias_t;

    #define E_INV (0.36787944117144233)

    /*******************************************************************************
    DiodeRectifier */
    template<typename T>
    struct DiodeRectifier
    {
        bias_t  bias;
        T   threshold;
        T   vt;
        T   scale;
        T   abs_coeff;
        T*  scratch;

        DiodeRectifier(bias_t bias, T thresh)
        {
                    
            /* Allocate scratch space */
            scratch = new T[4096];
            assert(scratch != nullptr);
            // Initialization
            this->bias = bias;
            threshold = thresh;
            abs_coeff = (bias == FULL_WAVE) ? 1.0 : 0.0;
            setThreshold(threshold);
        }        
        ~DiodeRectifier() {
            if(scratch) delete [] scratch;
        }
        void setThreshold(T thresh)
        {
            T scale = (T)1.0;
            threshold = LIMIT(std::fabs(thresh), 0.01, 0.9);
            scale = (1.0 - threshold);
            if(bias== REVERSE_BIAS)
            {
                scale *= -1.0;
                threshold *= -1.0;
            }
            threshold = threshold;
            vt = -0.1738 * threshold + 0.1735;
            scale = scale/(std::exp((1.0/vt) - 1.));
        }
        void ProcessBlock(size_t n, T * in_buffer, T * out_buffer)
        {
            T inv_vt = 1.0 / vt;            
            if (bias == FULL_WAVE)
            {
                VectorAbs(cratch, in_buffer, n_samples);
                VectorScalarMultiply(scratch, scratch, inv_vt, n_samples);
            }
            else
            {
                VectorScalarMultiply(scratch, in_buffer, inv_vt, n_samples);
            }
            VectorScalarAdd(scratch, >scratch, -1.0, n_samples);
            #pragma omp simd
            for (unsigned i = 0; i < n; ++i)
            {
                out_buffer[i] = std::exp(scratch[i]) * scale;
            }
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            return std::exp((in_sample/vt)-1) * scale;
        }
    };

    template<typename T>
    struct DiodeSaturator
    {
        bias_t  bias;
        T      amount;

        DiodeSaturator(bias_t bias, float amount)
        {           
           bias = bias;
           amount = amount;           
        }

        void setAmount(T a) {
            amount = 0.5*std::pow(a,0.5);
        }
        void setThreshold(T a) {
            amount = 0.5 * pow(a, 0.5);
        }
        void ProcessBlock(size_t n, T * in_buffer, T *out_buffer)
        {
            #pragma omp simd
            for (unsigned i = 0; i < n; ++i)
            {
                out_buffer[i] = in_buffer[i] - (amount * (F_EXP((in_buffer[i]/0.7) - 1.0) + E_INV));
            }
        }
        T Tick(T I, T A=1, T X=1, T Y=1)
        {
            return in_sample - (amount * (F_EXP((in_sample/0.7) - 1.0) + E_INV));
        }
    };
}