#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct VOSIM : public GeneratorProcessorPlugin<daisysp::VosimOscillator>
    {
        VOSIM(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::VosimOscillator>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FREQ,
            PORT_FORM1FREQ,
            PORT_FORM2FREQ,
            PORT_SHAPE
        };
        void setPort(int port, DspFloatType v) {
            switch(port)
            {
                case PORT_FREQ: this->SetFreq(v); break;
                case PORT_FORM1FREQ: this->SetForm1Freq(v); break;
                case PORT_FORM2FREQ: this->SetForm2Freq(v); break;
                case PORT_SHAPE: this->SetShape(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}