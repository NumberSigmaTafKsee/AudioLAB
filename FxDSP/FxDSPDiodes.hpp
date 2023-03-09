#pragma once

namespace DSPFX
{

    struct FxDiodeRectifier : public MonoFXProcessor
    {
        enum _bias_t
        {
        /** Pass positive signals, clamp netagive signals to 0 */
        FORWARD_BIAS,

        /** Pass negative signals, clamp positive signals to 0 */
        REVERSE_BIAS,

        /** Full-wave rectification. */
        FULL_WAVE
        };

        #ifdef DSPFLOATDOUBLE
        DiodeRectifierD * rect;
        #else
        DiodeRectifier * rect;
        #endif

        FxDiodeRectifier( bias_t bias, DspFloatType thresh) : MonoFXProcessor()
        {
            #ifdef DSPFLOATDOUBLE
            rect = DiodeRectifierInitD((bias_t)bias,thresh);
            #else
            rect = DiodeRectifierInit((bias_t)bias,thresh);
            #endif

            assert(rect != NULL);
        }
        ~FxDiodeRectifier() {
            #ifdef DSPFLOATDOUBLE
            if(rect) DiodeRectifierFreeD(rect);
            #else
            if(rect) DiodeRectifierFree(rect);
            #endif
        }
        void setThreshold(DspFloatType t) {
            #ifdef DSPFLOATDOUBLE
            Error_t err = DiodeRectifierSetThresholdD(rect,t);
            #else
            Error_t err = DiodeRectifierSetThreshold(rect,t);
            #endif
        }
        enum {
            PORT_THRESHOLD,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_THRESHOLD: setThreshold(v); break;
                default: printf("No Port %d\n",port);
            }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            Error_t err = DiodeRectifierProcessD(rect,out,in,n);
            #else
            Error_t err = DiodeRectifierProcess(rect,out,in,n);
            #endif
        }
        DspFloatType Tick(DspFloatType in) {
            #ifdef DSPFLOATDOUBLE
            return DiodeRectifierTickD(rect,in);
            #else
            return DiodeRectifierTick(rect,in);
            #endif
        }
    };

    struct FxDiodeSaturator : public MonoFXProcessor
    {
        enum _bias_t
        {
        /** Pass positive signals, clamp netagive signals to 0 */
        FORWARD_BIAS,

        /** Pass negative signals, clamp positive signals to 0 */
        REVERSE_BIAS,

        /** Full-wave rectification. */
        FULL_WAVE
        };
        #ifdef DSPFLOATDOUBLE
        DiodeSaturatorD * rect;
        #else
        DiodeSaturator * rect;
        #endif


        FxDiodeSaturator( bias_t bias, DspFloatType amt) : MonoFXProcessor()
        {
            #ifdef DSPFLOATDOUBLE
            rect = DiodeSaturatorInitD((bias_t)bias,amt);
            #else
            rect = DiodeSaturatorInit((bias_t)bias,amt);
            #endif
            assert(rect != NULL);
        }
        ~FxDiodeSaturator() {
            #ifdef DSPFLOATDOUBLE
            if(rect) DiodeSaturatorFreeD(rect);
            #else
            if(rect) DiodeSaturatorFree(rect);
            #endif
        }
        void setThreshold(DspFloatType t) {
            #ifdef DSPFLOATDOUBLE
            Error_t err = DiodeSaturatorSetAmountD(rect,t);
            #else
            Error_t err = DiodeSaturatorSetAmount(rect,t);
            #endif
        }
        enum {
            PORT_THRESHOLD,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_THRESHOLD: setThreshold(v); break;
                default: printf("No Port %d\n",port);
            }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            Error_t err = DiodeSaturatorProcessD(rect,out,in,n);
            #else
            Error_t err = DiodeSaturatorProcess(rect,out,in,n);
            #endif
        }
        DspFloatType Tick(DspFloatType in, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            #ifdef DSPFLOATDOUBLE
            return DiodeSaturatorTickD(rect,in);
            #else
            return DiodeSaturatorTick(rect,in);
            #endif
        }        
    };
}
