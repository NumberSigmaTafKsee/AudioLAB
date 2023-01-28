#pragma once

namespace DSPFX
{

    struct FxOnePoleFilter : FilterProcessor
    {
        #ifdef DSPFLOATDOUBLE
        OnePoleD * filter;
        #else
        OnePole * filter;
        #endif

        enum Type
        {
        /** Lowpass */
        LOWPASS,

        /** Highpass */
        HIGHPASS,

        /** Bandpass */
        BANDPASS,

        /** Allpass */
        ALLPASS,

        /** Notch */
        NOTCH,

        /** Peaking */
        PEAK,

        /** Low Shelf */
        LOW_SHELF,

        /** High Shelf */
        HIGH_SHELF,

        /** Number of Filter types */
        N_FILTER_TYPES
        };
        DspFloatType cutoff;
        // only low and highpass are valie
        FxOnePoleFilter(DspFloatType cut, int type=LOWPASS, DspFloatType sr=44100.0)
        : FilterProcessor()
        {
            cutoff = cut;
            #ifdef DSPFLOATDOUBLE
            filter = OnePoleInitD(cut,sr,(Filter_t)type);
            #else
            filter = OnePoleInit(cut,sr,(Filter_t)type);
            #endif
            assert(filter != NULL);
        }
        ~FxOnePoleFilter() {
            #ifdef DSPFLOATDOUBLE
            if(filter) OnePoleFreeD(filter);
            #else
            if(filter) OnePoleFree(filter);
            #endif
        }
        void flush() { 
            #ifdef DSPFLOATDOUBLE
            OnePoleFlushD(filter); 
            #else
            OnePoleFlush(filter); 
            #endif
        }
        void setType(int type) {
            #ifdef DSPFLOATDOUBLE
            OnePoleSetTypeD(filter,(Filter_t)type);
            #else
            OnePoleSetType(filter,(Filter_t)type);
            #endif
        }
        void setCutoff(DspFloatType c) {
            cutoff = c;
            #ifdef DSPFLOATDOUBLE
            OnePoleSetCutoffD(filter,c);
            #else
            OnePoleSetCutoff(filter,c);
            #endif
        }
        void setCoefficients(DspFloatType beta, DspFloatType alpha) {
            #ifdef DSPFLOATDOUBLE
            OnePoleSetCoefficientsD(filter,&beta, &alpha);
            #else
            OnePoleSetCoefficients(filter,&beta, &alpha);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #ifdef DSPFLOATDOUBLE
            OnePoleProcessD(filter,out,in,n);
            #else
            OnePoleProcess(filter,out,in,n);
            #endif
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType c = cutoff;            
            setCutoff(cutoff * X * Y);
            #ifdef DSPFLOATDOUBLE
            DspFloatType out = OnePoleTickD(filter,I);
            #else
            DspFloatType out = OnePoleTick(filter,I);
            #endif
            setCutoff(c);
            return out * A;
        }
        enum {
            PORT_CUTOFF,
            PORT_COEFFS,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;                
                default: printf("no port %d\n",port);
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b) {
            switch(port) {
                case PORT_COEFFS: setCoefficients(a,b); break;
                default: printf("no port %d\n",port);
            }
        }
        
    };
}
