#pragma once

namespace DSPFX
{

    struct FxPolySaturator : public FunctionProcessor
    {
        #ifdef DSPFLOATDOUBLE
        PolySaturatorD * sat;
        #else
        PolySaturator * sat;
        #endif

        DspFloatType N;

        FxPolySaturator(DspFloatType n) 
        : FunctionProcessor()
        {
            N = n;
            #ifdef DSPFLOATDOUBLE
            sat = PolySaturatorInitD(n);
            #else
            sat = PolySaturatorInit(n);
            #endif
            assert(sat != NULL);
        }
        ~FxPolySaturator() {
            #ifdef DSPFLOATDOUBLE
            if(sat) PolySaturatorFreeD(sat);
            #else
            if(sat) PolySaturatorFree(sat);
            #endif
        }
        void setN(DspFloatType n) {
            N = n;
            #ifdef DSPFLOATDOUBLE
            PolySaturatorSetND(sat,n);
            #else
            PolySaturatorSetN(sat,n);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            #ifdef DSPFLOATDOUBLE
            PolySaturatorProcessD(sat,out,in,n);
            #else
            PolySaturatorProcess(sat,out,in,n);
            #endif
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType t = N;
            setN(N*X*Y);
            #ifdef DSPFLOATDOUBLE
            DspFloatType r = PolySaturatorTickD(sat,I);
            #else
            DspFloatType r = PolySaturatorTick(sat,I);
            #endif
            setN(t);
            return A*r;
        }
        enum {
            PORT_N,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_N: setN(v); break;
                default: printf("no port %d\n",port);
            }
            
        }        
    };
}
