#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Drums
{
    struct AnalogBassDrum : public GeneratorProcessorPlugin<daisysp::AnalogBassDrum>
    {        
        AnalogBassDrum(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::AnalogBassDrum>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_TRIG,            
            PORT_SUSTAIN,
            PORT_ACCENT,
            PORT_FREQ,
            PORT_TONE,
            PORT_DECAY,
            PORT_ATTACKFM,
            PORT_SELFFM,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRIG: this->Trig(); break;                
                case PORT_SUSTAIN: this->SetSustain(v); break;
                case PORT_ACCENT: this->SetAccent(v); break;
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_TONE: this->SetTone(v); break;
                case PORT_DECAY: this->SetDecay(v); break;
                case PORT_ATTACKFM: this->SetAttackFmAmount(v); break;
                case PORT_SELFFM: this->SetSelfFmAmount(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1 )
        {
            return A*this->Process(I);
        }        
    };
}