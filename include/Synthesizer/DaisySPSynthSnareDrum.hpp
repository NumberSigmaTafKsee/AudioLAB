#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Drums
{
    struct SynthSnareDrum : public GeneratorProcessorPlugin<daisysp::SyntheticSnareDrum>
    {
        SynthSnareDrum(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::SyntheticSnareDrum>()
        {
            this->Init(sampleRate);
        }

        enum {
            PORT_TRIG,
            PORT_TRIGGER,
            PORT_SUSTAIN,
            PORT_ACCENT,
            PORT_FREQ,        
            PORT_DECAY,           
            PORT_FMAMT,
            PORT_SNAPPY,
        };
        bool trigger = false;
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRIG: this->Trig(); break;
                case PORT_TRIGGER: trigger = (bool)v; break;
                case PORT_SUSTAIN: this->SetSustain(v); break;
                case PORT_ACCENT: this->SetAccent(v); break;
                case PORT_FREQ: this->SetFreq(v); break;                
                case PORT_DECAY: this->SetDecay(v); break;                
                case PORT_FMAMT: this->SetFmAmount(v); break;
                case PORT_SNAPPY: this->SetSnappy(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(trigger);
        }
    };
}