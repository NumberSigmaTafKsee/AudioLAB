#pragma once

namespace DSPFX
{

    struct FxLRFilter
    {
        #ifdef DSPFLOATDOUBLE
        LRFilterD * filter;
        #else
        LRFilter * filter;
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
        DspFloatType cutoff, Q;
        int type = LOWPASS;
        FxLRFilter(int type, DspFloatType cut, DspFloatType q, DspFloatType sr = 44100) {
            #ifdef DSPFLOATDOUBLE
            filter = LRFilterInitD((Filter_t)type,cut,q,sr);
            #else
            filter = LRFilterInit((Filter_t)type,cut,q,sr);
            #endif
            assert(filter != NULL);
            cutoff = cut;
            Q = q;
            this->type = type;
        }
        ~FxLRFilter() {
            #ifdef DSPFLOATDOUBLE
            if(filter) LRFilterFreeD(filter);
            #else
            if(filter) LRFilterFree(filter);
            #endif
        }
        void flush() { 
            #ifdef DSPFLOATDOUBLE
            LRFilterFlushD(filter); 
            #else
            LRFilterFlush(filter);    
            #endif
        }
        void setParams(int type, DspFloatType cut, DspFloatType q) {
            this->type = type;
            cutoff = cut;
            Q = q;
            #ifdef DSPFLOATDOUBLE
            LRFilterSetParamsD(filter,(Filter_t)type,cutoff,Q);
            #else
            LRFilterSetParams(filter,(Filter_t)type,cutoff,Q);
            #endif
        }
        void setCutoff(DspFloatType cut) {
            cutoff = cut;
            #ifdef DSPFLOATDOUBLE
            LRFilterSetParamsD(filter,(Filter_t)type,cutoff,Q);
            #else
            LRFilterSetParams(filter,(Filter_t)type,cutoff,Q);
            #endif
        }
        void setQ(DspFloatType q) {
            Q = q;
            #ifdef DSPFLOATDOUBLE
            LRFilterSetParamsD(filter,(Filter_t)type,cutoff,Q);
            #else
            LRFilterSetParams(filter,(Filter_t)type,cutoff,Q);
            #endif
        }
        void setType(int type) {
            this->type = type;
            #ifdef DSPFLOATDOUBLE
            LRFilterSetParamsD(filter,(Filter_t)type,cutoff,Q);
            #else
            LRFilterSetParams(filter,(Filter_t)type,cutoff,Q);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            LRFilterProcessD(filter,out,in,n);
            #else
            LRFilterProcess(filter,out,in,n);
            #endif
        }
        enum {
            PORT_CUTOFF,
            PORT_Q,
            PORT_TYPE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
                case PORT_TYPE: setType(v); break;
                default: printf("No Port %d\n",port);
            }
        }
    };
}
