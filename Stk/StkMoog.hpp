#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    struct Moog : public GeneratorProcessorPlugin<stk::Moog>
    {
        Moog() : GeneratorProcessorPlugin<stk::Moog>()
        {

        }
        enum {
            PORT_FREQ,
            PORT_NOTEON,
            PORT_MODSPEED,
            PORT_MODDEPTH,
            PORT_CC,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: this->setFrequency(v); break;
                case PORT_MODSPEED: this->setModulationSpeed(v); break;
                case PORT_MODDEPTH: this->setModulationDepth(v); break;
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b) {
            switch(port) {
                case PORT_NOTEON: this->noteOn(a,b); break;
                case PORT_CC: this->controlChange(a,b); break;
            }
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick();
        }
    };
}