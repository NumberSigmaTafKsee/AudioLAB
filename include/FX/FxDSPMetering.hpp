#pragma once

namespace DSPFX
{

    struct FxMetering
    {

        enum
        {
            FULL_SCALE,
            K_12,
            K_14,
            K_20
        };

        static DspFloatType phase_correlation(DspFloatType * left, DspFloatType * right, size_t n) {
            #ifdef DSPFLOATDOUBLE
            return ::phase_correlationD(left,right,n);
            #else
            return ::phase_correlation(left,right,n);
            #endif
        }
        static DspFloatType balance(DspFloatType * left, DspFloatType * right, size_t n) {
            #ifdef DSPFLOATDOUBLE
            return ::balanceD(left,right,n);
            #else
            return ::balance(left,right,n);
            #endif
        }
        static DspFloatType vu_peak(DspFloatType * signal, size_t n, int scale) {
            #ifdef DSPFLOATDOUBLE
            return ::vu_peakD(signal,n,(MeterScale)scale);
            #else
            return ::vu_peak(signal,n,(MeterScale)scale);
            #endif
        }
    };
}
