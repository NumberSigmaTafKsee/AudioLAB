#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Envelopes
{

    struct Seg : public GeneratorProcessor, public gam::Seg
    {
        Seg(DspFloatType len=0.5, DspFloatType start=1, DspFloatType end=0, DspFloatType phase=0) :
            GeneratorProcessor(),
            gam::Seg(len,start,end,phase)
        {
            
        }
        enum {
            PORT_FREQ,
            PORT_LENGTH,
            PORT_PERIOD,
            PORT_PHASE,
            PORT_RESET,
            PORT_DONE,
            PORT_VAL,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->freq(v); break;
                case PORT_LENGTH: this->length(v); break;
                case PORT_PERIOD: this->period(v); break;
                case PORT_PHASE: this->phase(v); break;
                case PORT_RESET: this->reset(); break;
            }
        }
        DspFloatType getPort(int port) {
            switch(port) {
                case PORT_DONE: return this->done();
                case PORT_VAL: return this->val();
            }
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return (*this)();
        }
    };

    struct SegExp : public GeneratorProcessor, public gam::SegExp
    {
        SegExp(DspFloatType len=0.5, DspFloatType crv=-3, DspFloatType start=1, DspFloatType end=0) :
            GeneratorProcessor(),
            gam::Seg(len,crv,start,end)
        {
            
        }
        enum {
            PORT_DONE,
            PORT_CURVE,
            PORT_PERIOD,
            PORT_RESET,            
            PORT_SET,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_PERIOD: this->period(v); break;
                case PORT_RESET: this->reset(); break;
                case PORT_CURVE: this->curve(v); break;
            }
        }
        void setPort2(int port, DspFloatType x, DspFloatType y) {
            switch(port) {
                case PORT_SET: this->set(x,y); break;
            }
        }
        DspFloatType getPort(int port) {
            switch(port) {
                case PORT_DONE: return this->done();                
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return (*this)();
        }
    };
}