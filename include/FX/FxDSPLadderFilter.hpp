#pragma once

namespace DSPFX
{

    struct FxLadderFilter 
    {
        #ifdef DSPFLOATDOUBLE
        LadderFilterD * filter;
        #else
        LadderFilter * filter;
        #endif

        FxLadderFilter(DspFloatType sr) {
            #ifdef DSPFLOATDOUBLE
            filter = LadderFilterInitD(sr);
            #else
            filter = LadderFilterInit(sr);
            #endif
            assert(filter != NULL);
        }
        ~FxLadderFilter() {
            #ifdef DSPFLOATDOUBLE
            if(filter) LadderFilterFreeD(filter);
            #else
            if(filter) LadderFilterFree(filter);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            LadderFilterProcessD(filter,out,in,n);
            #else
            LadderFilterProcess(filter,out,in,n);
            #endif
        }
        void setCutoff(DspFloatType c) {
            #ifdef DSPFLOATDOUBLE
            LadderFilterSetCutoffD(filter,c);
            #else
            LadderFilterSetCutoff(filter,c);
            #endif
        }
        void setResonance(DspFloatType q) {
            #ifdef DSPFLOATDOUBLE
            LadderFilterSetResonanceD(filter,q);
            #else
            LadderFilterSetResonance(filter,q);
            #endif
        }
        void setTemperature(DspFloatType t) {
            #ifdef DSPFLOATDOUBLE
            LadderFilterSetTemperatureD(filter,t);
            #else
            LadderFilterSetTemperature(filter,t);
            #endif
        }
        enum {
            PORT_CUTOFF,
            PORT_RESONANCE,
            PORT_TEMP,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_RESONANCE: setResonance(v); break;
                case PORT_TEMP: setTemperature(v); break;
                default: printf("No Port %d\n",port);
            }
        }
    };
}
