#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::Delay
{
    struct Delay1 : public FunctionProcessor, public gam::Delay1<DspFloatType>
    {
        Delay1(DspFloatType value=0) 
        : FunctionProcessor(),
          gam::Delay1<DspFloatType>(value)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return (*this)(I);
        }
    };

    struct Delay2 : public FunctionProcessor, public gam::Delay2<DspFloatType>
    {
        Delay2(DspFloatType value=0) : 
        FunctionProcessor(),
        gam::Delay2<DspFloatType>(value)
        {

        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return (*this)(I);
        }
    };

    struct Delay : public FunctionProcessor, public gam::Delay<DspFloatType>
    {
        Delay(DspFloatType value=0) : FunctionProcessor(),
        gam::Delay<DspFloatType>(value)
        {

        }
        enum {
            PORT_DELAY,
            PORT_DELAYSAMPLES,
            PORT_DELAYSAMPLESR,
            PORT_DELAYUNIT,
            PORT_FREQ,
            PORT_IPLTYPE,
            PORT_MAXDELAY
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DELAY: this->delay(v); break;
                case PORT_DELAYSAMPLES: this->delaySamples(v); break;
                case PORT_DELAYSAMPLESR: this->delaySamplesR(v); break;
                case PORT_DELAYUNIT: this->delayUnit(v); break;
                case PORT_FREQ: this->freq(v); break;
                case PORT_IPLTYPE: this->ipolType((gam::ipl::Type)v); break;
                case PORT_MAXDELAY: this->maxDelay(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {            
            return A*(*this)(I);
        }
    };

    struct Multitap : public FunctionProcessor, public gam::Multitap<DspFloatType>
    {
        Multitap(DspFloatType delay, unsigned numTaps) : 
        FunctionProcessor(),
        gam::Multitap<DspFloatType>(delay,numTaps)
        {

        }
        enum {
            PORT_DELAY,
            PORT_FREQ,
            PORT_TAPS,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TAPS: this->taps(v); break;
            }
        }
        void setPort2(int port, DspFloatType x, DspFloatType y) {
            switch(port) {
                case PORT_DELAY: this->delay(x,y); break;
                case PORT_FREQ: this->freq(x,y); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return A*(*this)(I);
        }
    };
}
