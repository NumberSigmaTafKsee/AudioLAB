#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Control
{
    struct ADenv : public GeneratorProcessorPlugin<daisysp::AdEnv>
    {
        ADenv(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::AdEnv>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_TRIGGER,
            PORT_ATTACK,
            PORT_DECAY,
            PORT_CURVE,
            PORT_MIN,
            PORT_MAX,            
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRIGGER: if(v == 1.0) this->Trigger(); break;
                case PORT_ATTACK: this->SetTime(daisysp::AdEnvSegment::ADENV_SEG_ATTACK,v); break;
                case PORT_DECAY: this->SetTime(daisysp::AdEnvSegment::ADENV_SEG_DECAY,v); break;
                case PORT_CURVE: this->SetCurve(v); break;
                case PORT_MIN: this->SetMin(v); break;
                case PORT_MAX: this->SetMax(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return this->Process();
        }
    };
}