#pragma once

#include "fxdsp.hpp"

namespace FXDSP
{
    template<typename T>
    void StereoToMono(T*        dest,
                const T*  left,
                const T*  right,
                unsigned      length)
    {
        T scale = SQRT_TWO_OVER_TWO;
        VectorVectorSumScale(dest, left, right, &scale, length);
    }

    
    template<typename T>
    void MonoToStereo(T*         left,
                T*         right,
                const T*   mono,
                unsigned       length)
    {
        T scale = SQRT_TWO_OVER_TWO;
        VectorScalarMultiply(left, mono, scale, length);
        CopyBuffer(right, left, length);
    }   
}