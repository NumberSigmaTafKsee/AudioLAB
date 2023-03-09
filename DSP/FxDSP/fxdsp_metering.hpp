#pragma once

#include <cmath>
#include <cfloat>
#include "fxdsp.hpp"

/* pow(10, (-12./20.)) */
#define K12_REF (0.25118864315095801309)

/* pow(10, (-12./20.)); */
#define K14_REF (0.19952623149688797355)

/* pow(10, (-20./20.)); */
#define K20_REF (0.1)


static const float ref[] = {1.0, (float)K12_REF, (float)K14_REF, (float)K20_REF};
static const double refD[] = {1.0, K12_REF, K14_REF, K20_REF};

namespace FXDSP
{
    typedef enum
    {
        FULL_SCALE,
        K_12,
        K_14,
        K_20
    } MeterScale;

    template<typename T>
    T phase_correlation(T* left, T* right, unsigned n_samples)
    {
        T product = (T)0.0;
        T lsq = (T)0.0;
        T rsq = (T)0.0;
        T denom = (T)0.0;

        #pragma omp simd aligned(left,right)
        for (unsigned i = 0; i < n_samples; ++i)
        {
            T left_sample = left[i];
            T right_sample = right[i];

            product += left_sample * right_sample;
            lsq += left_sample * left_sample;
            rsq += right_sample * right_sample;
        }

        denom = lsq * rsq;

        if (denom > 0.0)
        {
            return product / std::sqrt(denom);
        }
        else
        {
            return 1.0;
        }
    }

    template<typename T>
    T balance(T* left, T* right, unsigned n_samples)
    {
        T r = (T)0.0;
        T l = (T)0.0;
        T rbuf[n_samples];
        T lbuf[n_samples];
        VectorPower(rbuf, right, 2.0, n_samples);
        VectorPower(lbuf, left, 2.0, n_samples);
        r = VectorSum(rbuf, n_samples);
        l = VectorSum(lbuf, n_samples);
        return  (r - l) / ((r + l) + FLT_MIN);
    }

    template<typename T>
    T vu_peak(T* signal, unsigned n_samples, MeterScale scale)
    {
        T peak = VectorMax(signal, n_samples);
        return 20.0 * std::log10(peak / ref[scale]);
    }
}
