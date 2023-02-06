#pragma once

#include "DaisySP.hpp"

namespace DaisySP::PhysicalModels
{
    struct ModalVoice : public GeneratorProcessorPlugin<daisysp::ModalVoice>
    {
        ModalVoice(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::ModalVoice>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_TRIG,
            PORT_SUSTAIN,
            PORT_FREQ,
            PORT_ACCENT,
            PORT_STRUCTURE,
            PORT_BRIGHTNESS,
            PORT_DAMP,
            PORT_AUX
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_TRIG: this->Trig(); break;
                case PORT_SUSTAIN: this->SetSustain(v); break;
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_ACCENT: this->SetAccent(v); break;
                case PORT_STRUCTURE: this->SetStructure(v); break;
                case PORT_BRIGHTNESS: this->SetBrightness(v); break;
                case PORT_DAMP: this->SetDamping(v); break;                
            }
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A* this->Process(I);
        }
    };
}