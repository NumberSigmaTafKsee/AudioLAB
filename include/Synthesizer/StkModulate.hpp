#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    struct Modulate : public GeneratorProcessorPlugin<stk::Modulate>
    {
        Modulate() : GeneratorProcessorPlugin<stk::Modulate>()
        {

        }
        enum {
            PORT_RESET,
            PORT_VIBRATORATE,
            PORT_VIBRATOGAIN,
            PORT_RANDOMRATE,
            PORT_RANDOMGAIN,
            PORT_LASTOUT,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_RESET: this->reset(); break;
                case PORT_VIBRATORATE: this->setVibratoRate(v); break;
                case PORT_VIBRATOGAIN: this->setVibratoGain(v); break;
                case PORT_RANDOMRATE: this->setRandomRate(v); break;
                case PORT_RANDOMGAIN: this->setRandomGain(v); break;                
            }
        }
        DspFloatType getPort(int port) {
            if(port == PORT_LASTOUT) return this->lastOut();
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick();
        }
    };
}