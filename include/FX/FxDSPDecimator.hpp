#pragma once

namespace DSPFX
{

    struct FxDecimator : public MonoFXProcessor
    {
        #ifdef DSPFLOATDOUBLE
        DecimatorD * deci;
        #else
        Decimator * deci;
        #endif

        /** Resampling Factor constants */
        enum factor
        {
            /** 2x resampling */
            X2 = 0,

            /** 4x resampling */
            X4,

            /** 8x resampling */
            X8,

            /** 16x resampling */
            /*X16,*/

            /** number of resampling factors */
            N_FACTORS
        };
        FxDecimator(factor fac) : MonoFXProcessor() {
            #ifdef DSPFLOATDOUBLE
            deci = DecimatorInitD((ResampleFactor_t)fac);
            #else
            deci = DecimatorInit((ResampleFactor_t)fac);
            #endif
            assert(deci != NULL);
        }
        ~FxDecimator() {
            #ifdef DSPFLOATDOUBLE
            if(deci) DecimatorFreeD(deci);
            #else
            if(deci) DecimatorFree(deci);
            #endif
        }
        void flush() {
            #ifdef DSPFLOATDOUBLE
            DecimatorFlushD(deci);
            #else
            DecimatorFlush(deci);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            Error_t err = DecimatorProcessD(deci,out,in,n);
            #else
            Error_t err = DecimatorProcess(deci,out,in,n);
            #endif
        }
    };

}

