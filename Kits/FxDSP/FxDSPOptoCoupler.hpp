#pragma once

namespace DSPFX
{
    struct FxOptoCoupler
    {
        enum 
        {
            /** Light-Dependent-Resistor output. Based
             on Vactrol VTL series. datasheet:
            http://pdf.datasheetcatalog.com/datasheet/perkinelmer/VT500.pdf
            Midpoint Delay values:
            Turn-on delay:   ~10ms
            Turn-off delay:  ~133ms
            */
            OPTO_LDR,

            /** TODO: Add Phototransistor/voltage output opto model*/
            OPTO_PHOTOTRANSISTOR
        };

        #ifdef DSPFLOATDOUBLE
        OptoD * oc;
        #else
        Opto * oc;
        #endif

        DspFloatType delay;

        FxOptoCoupler(int type, DspFloatType delay, DspFloatType sr=44100.0)
        {
            this->delay = delay;
            #ifdef DSPFLOATDOUBLE
            oc = OptoInitD((Opto_t)type,delay,sr);
            #else
            oc = OptoInit((Opto_t)type,delay,sr);
            #endif
            assert(oc != NULL);
        }
        ~FxOptoCoupler() {
            #ifdef DSPFLOATDOUBLE
            OptoFreeD(oc);
            #else
            OptoFree(oc);
            #endif
        }
        void setDelay(DspFloatType d) {
            delay = d;
            #ifdef DSPFLOATDOUBLE
            OptoSetDelayD(oc,d);
            #else
            OptoSetDelay(oc,d);
            #endif
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #ifdef DSPFLOATDOUBLE
            OptoProcessD(oc,out,in,n);
            #else
            OptoProcess(oc,out,in,n);
            #endif
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType d = delay;
            setDelay(d*X*Y);
            #ifdef DSPFLOATDOUBLE
            DspFloatType r = OptoTickD(oc,I);
            #else
            DspFloatType r = OptoTick(oc,I);
            #endif
            setDelay(d);
            return A*r;
        }
        enum {
            PORT_DELAY,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DELAY: setDelay(v); break;
                default: printf("no port %d\n",port);
            }
        }
    };
}
