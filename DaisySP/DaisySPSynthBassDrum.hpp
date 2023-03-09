#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Drums
{
    struct SynthBassDrum : public GeneratorProcessorPlugin<daisysp::SyntheticBassDrum>
    {
        SynthBassDrum(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::SyntheticBassDrum>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_TRIG,
            PORT_TRIGGER,
            PORT_SUSTAIN,
            PORT_ACCENT,
            PORT_FREQ,
            PORT_TONE,
            PORT_DECAY,
            PORT_DIRT,
            PORT_FMENVAMT,
            PORT_FMENVDCY,
        };
        bool trigger = false;
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRIG: this->Trig(); break;
                case PORT_TRIGGER: trigger = (bool)v; break;
                case PORT_SUSTAIN: this->SetSustain(v); break;
                case PORT_ACCENT: this->SetAccent(v); break;
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_TONE: this->SetTone(v); break;
                case PORT_DECAY: this->SetDecay(v); break;
                case PORT_DIRT: this->SetDirtiness(v); break;
                case PORT_FMENVAMT: this->SetFmEnvelopeAmount(v); break;
                case PORT_FMENVDCY: this->SetFmEnvelopeDecay(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process(trigger);
        }
    };
}