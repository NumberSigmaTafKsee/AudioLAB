#pragma once

namespace DSPFX
{

    struct FxRMSEstimator
    {
        #ifdef DSPFLOATDOUBLE
        RMSEstimatorD *rms;
        #else
        RMSEstimator *rms;
        #endif

        FxRMSEstimator(DspFloatType avg, DspFloatType sr=44100)
        {
            #ifdef DSPFLOATDOUBLE
            rms = RMSEstimatorInitD(avg,sr);
            #else
            rms = RMSEstimatorInit(avg,sr);
            #endif
            assert(rms != NULL);
        }
        ~FxRMSEstimator()
        {
            #ifdef DSPFLOATDOUBLE
            if(rms) RMSEstimatorFreeD(rms);
            #else
            if(rms) RMSEstimatorFree(rms);
            #endif
        }
        void flush() { 
            #ifdef DSPFLOATDOUBLE
            RMSEstimatorFlushD(rms); 
            #else
            RMSEstimatorFlush(rms); 
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            RMSEstimatorProcessD(rms,out,in,n);
            #else
            RMSEstimatorProcess(rms,out,in,n);
            #endif
        }
        // this is never modulated
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            #ifdef DSPFLOATDOUBLE
            return RMSEstimatorTickD(rms,I);
            #else
            return RMSEstimatorTick(rms,I);
            #endif
        }
    };
}
