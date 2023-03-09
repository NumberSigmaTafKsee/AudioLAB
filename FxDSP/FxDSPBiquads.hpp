#pragma once

#include "FxDSP.hpp"
#include <cassert>

namespace DSPFX
{
    struct FxBiquadFilter : public MonoFXProcessor
    {
        #ifdef DSPFLOATDOUBLE
        BiquadFilterD * filter;
        #else
        BiquadFilter * filter;
        #endif

        FxBiquadFilter(const DspFloatType * Bc, const DspFloatType * Ac) : MonoFXProcessor()
        {
            #ifdef DSPFLOATDOUBLE
            filter = BiquadFilterInitD(Bc,Ac);
            #else
            filter = BiquadFilterInit(Bc,Ac);
            #endif
            assert(filter != NULL);
        }    
        ~FxBiquadFilter() {
            
            #ifdef DSPFLOATDOUBLE
            if(filter) BiquadFilterFreeD(filter);
            #else
            if(filter) BiquadFilterFree(filter);
            #endif
        }
        void flush() {
            #ifdef DSPFLOATDOUBLE
            BiquadFilterFlushD(filter);
            #else
            BiquadFilterFlush(filter);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #ifdef DSPFLOATDOUBLE
            BiquadFilterProcessD(filter,out,in,n);
            #else
            BiquadFilterProcess(filter,out,in,n);
            #endif
        }
        DspFloatType Tick(DspFloatType In, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            #ifdef DSPFLOATDOUBLE
            return BiquadFilterTickD(filter,In);
            #else
            return BiquadFilterTick(filter,In);
            #endif
        }
        void updateCoefficients(const DspFloatType * B, const DspFloatType * A) {
            #ifdef DSPFLOATDOUBLE
            BiquadFilterUpdateKernelD(filter,B,A);
            #else
            BiquadFilterUpdateKernel(filter,B,A);
            #endif
        }
        enum {
            PORT_COEFF,
        };
        void setPortV(int port, const std::vector<DspFloatType> & c) {
            if(port == PORT_COEFF) updateCoefficients(&c[0],&c[3]);
            else printf("No Port %d\n",port);
        }
    };

    struct FxRBJFilter : public MonoFXProcessor
    {
        
        #ifdef DSPFLOATDOUBLE
        RBJFilterD * filter;
        #else
        RBJFilter * filter;
        #endif

        DspFloatType Cutoff,Q;
        
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

        FxRBJFilter(int type,DspFloatType cutoff, DspFloatType sampleRate) : MonoFXProcessor()
        {
            Cutoff = cutoff;
            Q      = 0.5;
            #ifdef DSPFLOATDOUBLE
            filter = RBJFilterInitD((Filter_t)type,cutoff,sampleRate);
            #else
            filter = RBJFilterInit((Filter_t)type,cutoff,sampleRate);
            #endif
            assert(filter != NULL);
        }    
        ~FxRBJFilter() {
            #ifdef DSPFLOATDOUBLE
            if(filter) RBJFilterFreeD(filter);
            #else
            if(filter) RBJFilterFree(filter);
            #endif
        }
        void setType(Type type) {
            #ifdef DSPFLOATDOUBLE
            RBJFilterSetTypeD(filter,(Filter_t)type);
            #else
            RBJFilterSetType(filter,(Filter_t)type);
            #endif
        }
        void setCutoff(DspFloatType cut) {
            Cutoff = cut;
            #ifdef DSPFLOATDOUBLE
            RBJFilterSetCutoffD(filter,cut);
            #else
            RBJFilterSetCutoff(filter,cut);
            #endif

        }
        void setQ(DspFloatType q) {
            Q = q;
            #ifdef DSPFLOATDOUBLE
            RBJFilterSetQD(filter,Q);
            #else
            RBJFilterSetQ(filter,Q);
            #endif
        }
        void flush() {
            #ifdef DSPFLOATDOUBLE
            RBJFilterFlushD(filter);
            #else
            RBJFilterFlush(filter);
            #endif
        }
        enum {
            PORT_TYPE,
            PORT_CUTOFF,
            PORT_Q,
            PORT_FLUSH,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TYPE: setType((Type)v); break;
                case PORT_CUTOFF: setCutoff(v); break;
                case PORT_Q: setQ(v); break;
                case PORT_FLUSH: flush(); break;
            }
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #ifdef DSPFLOATDOUBLE
            RBJFilterProcessD(filter,out,in,n);
            #else
            RBJFilterProcess(filter,out,in,n);
            #endif
        }          
    };
}
