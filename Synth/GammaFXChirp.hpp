#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::FX
{
    struct Chirp : public FunctionProcessorPlugin<gam::Chirp<DspFloatType>>
    {
        Chirp() : FunctionProcessorPlugin<gam::Chirp<DspFloatType>>()
        {

        }
        enum {
            PORT_DECAY,
            PORT_FREQ,
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_DECAY: this->decay(v); break;
                case PORT_RESET: this->reset(v); break;
            }
        }
        void setPort2(int port, DspFloatType x, DspFloatType y) {
            switch(port) {
                case PORT_FREQ: this->freq(x,y);
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*(*this)();
        }
    };
}