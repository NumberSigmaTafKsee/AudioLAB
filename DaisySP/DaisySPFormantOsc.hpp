#pragma once

#include "DaisySP.hpp"

namespace DaisySP::Synthesis
{
    struct FormantOsc : public GeneratorProcessorPlugin<daisysp::FormantOscillator>
    {
        FormantOsc(DspFloatType sampleRate=44100.0) : GeneratorProcessorPlugin<daisysp::FormantOscillator>()
        {
            this->Init(sampleRate);
        }
        enum {
            PORT_FORMANT_FREQ,
            PORT_CARRIER_FREQ,
            PORT_PHASE_SHIFT,
        };
        void setPort(int port, DspFloatType v)
        {
            switch(port) {
            case PORT_FORMANT_FREQ: this->SetFormantFreq(v); break;
            case PORT_CARRIER_FREQ: this->SetCarrierFreq(v); break;
            case PORT_PHASE_SHIFT: this->SetPhaseShift(v); break;
            }
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1) {
            return A*this->Process();
        }
    };
}