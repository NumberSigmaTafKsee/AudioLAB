#pragma once

#include "StkHeaders.hpp"

namespace Stk::Envelopes
{
    struct ADSR : public GeneratorProcessorPlugin<stk::ADSR>
    {
        ADSR() 
        : GeneratorProcessorPlugin<stk::ADSR>()
        {

        }
        ADSR(DspFloatType a, DspFloatType d, DspFloatType s, DspFloatType r)
        : GeneratorProcessorPlugin<stk::ADSR>()          
        {
            this->setAllTimes(a,d,s,r);
        }

        enum {
            PORT_ATTACKRATE,
            PORT_DECAYRATE,
            PORT_SUSTAIN,
            PORT_RELEASERATE,
            PORT_ATTACKTIME,
            PORT_DECAYTIME,
            PORT_RELEASETIME,
            PORT_TARGET,
            PORT_VALUE,
        };
        void setPort(int port, DspFloatType value) {
            switch(port) {
                case PORT_SUSTAIN: this->setSustainLevel(value); break;
                case PORT_ATTACKRATE: this->setAttackRate(value); break;
                case PORT_DECAYRATE: this->setDecayRate(value); break;
                case PORT_RELEASERATE: this->setReleaseRate(value); break;
                case PORT_ATTACKTIME: this->setAttackTime(value); break;
                case PORT_DECAYTIME: this->setDecayTime(value); break;
                case PORT_RELEASETIME: this->setReleaseTime(value); break;                                
                case PORT_TARGET: this->setTarget(value); break;
                case PORT_VALUE: this->setValue(value); break;
            }
        }
        void noteOn() {
            this->keyOn();
        }
        void noteOff() {
            this->keyOff();
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return this->tick();
        }
    };
}