#pragma once

#include "StkHeaders.hpp"

namespace Stk 
{
     // todo: ports
    struct FMVoices : public GeneratorProcessorPlugin<stk::FMVoices>
    {
        FMVoices() : GeneratorProcessorPlugin<stk::FMVoices>()
        {

        }
        enum {
            PORT_FREQ,
            PORT_RATIO,
            PORT_GAIN,
            PORT_MODSPEED,
            PORT_MODDEPTH,
            PORT_CONTROL1,
            PORT_CONTROL2,
            PORT_KEYON,
            PORT_KEYOFF,
            PORT_NOTEOFF,
            PORT_CONTROLCHANGE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_KEYON: this->keyOn(); break;
                case PORT_KEYOFF: this->keyOff(); break;
                case PORT_NOTEOFF: this->noteOff(v); break;
                case PORT_FREQ: this->setFrequency(v); break;
                case PORT_MODSPEED: this->setModulationSpeed(v); break;
                case PORT_MODDEPTH: this->setModulationDepth(v); break;
                case PORT_CONTROL1: this->setControl1(v); break;
                case PORT_CONTROL2: this->setControl2(v); break;
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b) {
            switch(port) {
                case PORT_RATIO: this->setRatio(a,b); break;
                case PORT_GAIN: this->setGain(a,b); break;
                case PORT_CONTROLCHANGE: this->controlChange(a,b); break;
            }
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->tick();
        }
    };
}