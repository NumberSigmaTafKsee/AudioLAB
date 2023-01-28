#pragma once

namespace DSPFX
{

    struct FxMultiBandFilter 
    {
        #ifdef DSPFLOATDOUBLE
        MultibandFilterD * filter;
        #else
        MultibandFilter * filter;
        #endif

        FxMultiBandFilter(DspFloatType low, DspFloatType high, DspFloatType sr = 44100.0)        
        {
            #ifdef DSPFLOATDOUBLE
            filter = MultibandFilterInitD(low,high,sr);
            #else
            filter = MultibandFilterInit(low,high,sr);
            #endif
            assert(filter != NULL);
        }
        ~FxMultiBandFilter() {
            #ifdef DSPFLOATDOUBLE
            if(filter) MultibandFilterFreeD(filter);
            #else
            if(filter) MultibandFilterFree(filter);
            #endif
        }
        void flush() { 
            #ifdef DSPFLOATDOUBLE
            MultibandFilterFlushD(filter); 
            #else
            MultibandFilterFlush(filter); 
            #endif
        }
        void setHighCutoff(DspFloatType c) {
            #ifdef DSPFLOATDOUBLE
            MultibandFilterSetHighCutoffD(filter,c);
            #else
            MultibandFilterSetHighCutoff(filter,c);
            #endif
        }
        void setLowCutoff(DspFloatType c) {
            #ifdef DSPFLOATDOUBLE
            MultibandFilterSetLowCutoffD(filter,c);
            #else
            MultibandFilterSetLowCutoff(filter,c);
            #endif
        }
        void update(DspFloatType low, DspFloatType high) {
            #ifdef DSPFLOATDOUBLE
            MultibandFilterUpdateD(filter,low,high);
            #else
            MultibandFilterUpdate(filter,low,high);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * low, DspFloatType * mid, DspFloatType * high) {
            #ifdef DSPFLOATDOUBLE
            MultibandFilterProcessD(filter,low,mid,high,in,n);
            #else
            MultibandFilterProcess(filter,low,mid,high,in,n);
            #endif
        }
        enum {
            PORT_HIGHCUT,
            PORT_LOWCUT,
            PORT_UPDATE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_HIGHCUT: setHighCutoff(v); break;
                case PORT_LOWCUT: setLowCutoff(v); break;
                default: printf("No Port %d\n",port);
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b) {
            switch(port) {
                case PORT_HIGHCUT: update(a,b); break;                
                default: printf("No Port %d\n",port);
            }
        }
    };

}
