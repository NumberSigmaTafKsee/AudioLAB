#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Control
{
    struct ADSR : public GeneratorProcessorPlugin<daisysp::Adsr>
    {
        bool gate = false;
        ADSR(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::Adsr>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_RETRIGGER,            
            PORT_ATTACKTIME,
            PORT_DECAYTIME,
            PORT_RELEASETIME,
            PORT_SUSTAINLEVEL,
            PORT_GATE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_RETRIGGER: this->Retrigger((bool)v); break;
                case PORT_ATTACKTIME: this->SetAttackTime(v); break;
                case PORT_DECAYTIME: this->SetDecayTime(v); break;
                case PORT_RELEASETIME: this->SetReleaseTime(v); break;
                case PORT_SUSTAINLEVEL: this->SetSustainLevel(v); break;
                case PORT_GATE: gate = (bool)v; break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return this->Process(gate);
        }
    };

}